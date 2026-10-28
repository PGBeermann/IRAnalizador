// IRAnalizador — servidor principal (Node.js + Express)
import express from 'express';
import helmet from 'helmet';
import compression from 'compression';
import multer from 'multer';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import { initChem, createContext, rdkitVersion } from './src/chem.js';
import { simulate } from './src/simulator.js';
import { explainStructure, explainPeaks, explainMissing } from './src/explain.js';
import {
  parseSpectrumFile, normalizeSpectrum, pickPeaks, assignPeaks, similarity, diagnose, phaseFromState
} from './src/analysis.js';
import { searchSpectra, getNistSpectrum, resolveCompound } from './src/external.js';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const PORT = Number(process.env.PORT) || 3000;
const app = express();

app.set('trust proxy', 1);
app.disable('x-powered-by');
// CSP desactivada: JSME (GWT) genera e inyecta scripts dinámicamente.
app.use(helmet({ contentSecurityPolicy: false, crossOriginEmbedderPolicy: false }));
app.use(compression());
app.use(express.json({ limit: '2mb' }));

// ---------------------------------------------------------------- estáticos
app.use('/vendor/plotly', express.static(path.join(__dirname, 'node_modules/plotly.js-dist-min'), { maxAge: '7d' }));
app.use('/vendor/jsme', express.static(path.join(__dirname, 'node_modules/jsme-editor'), { maxAge: '7d' }));
app.use(express.static(path.join(__dirname, 'public'), { maxAge: 0, extensions: ['html'] }));

// ---------------------------------------------------------------- límite simple de peticiones
const hits = new Map();
function rateLimit(max, windowMs = 60e3) {
  return (req, res, next) => {
    const k = `${req.ip}:${req.path}`;
    const now = Date.now();
    const e = hits.get(k) || { n: 0, t: now + windowMs };
    if (now > e.t) { e.n = 0; e.t = now + windowMs; }
    e.n++;
    hits.set(k, e);
    if (hits.size > 5000) hits.clear();
    if (e.n > max) return res.status(429).json({ error: 'Demasiadas solicitudes. Espere un minuto e inténtelo de nuevo.' });
    next();
  };
}

const wrap = (fn) => (req, res, next) => Promise.resolve(fn(req, res, next)).catch(next);
const PHASES = new Set(['gas', 'condensada']);

function moleculeInfo(ctx, highlight = []) {
  return {
    smiles: ctx.smiles,
    formula: ctx.formula,
    mw: Math.round(ctx.mw * 100) / 100,
    exactMass: Math.round(ctx.exactMass * 10000) / 10000,
    ihd: ctx.ihd,
    inchi: ctx.inchi,
    inchikey: ctx.inchikey,
    svg: ctx.svg(highlight)
  };
}

function runSimulation(structure, phase) {
  const ctx = createContext(structure);
  try {
    const sim = simulate(ctx, phase);
    return { ctx, sim, explanation: explainStructure(ctx, sim) };
  } catch (e) {
    ctx.dispose();
    throw e;
  }
}

/** Análisis común de un espectro experimental. */
function analyzeRaw(raw, { structure, phase, yMode = 'auto' }) {
  const spec = normalizeSpectrum(raw, yMode);
  const ph = PHASES.has(phase) ? phase : phaseFromState(raw.meta?.state);
  const peaks = pickPeaks(spec);
  const out = {
    meta: raw.meta,
    phase: ph,
    spectrum: spec,
    diagnosis: diagnose(peaks)
  };
  if (structure) {
    const { ctx, sim, explanation } = runSimulation(structure, ph);
    try {
      const assigned = assignPeaks(peaks, sim.bands, ph, spec.range);
      out.molecule = moleculeInfo(ctx);
      out.simulation = { spectrum: sim.spectrum, bands: sim.bands };
      out.groups = explanation.groups;
      out.explanation = explanation;
      out.peaks = assigned.peaks;
      out.peakText = explainPeaks(assigned, { phase: ph });
      out.missing = explainMissing(assigned.missing, ph);
      out.similarity = similarity(sim.spectrum, spec);
    } finally {
      ctx.dispose();
    }
  } else {
    const assigned = assignPeaks(peaks, null, ph, spec.range);
    out.peaks = assigned.peaks;
    out.peakText = explainPeaks(assigned, { phase: ph });
    out.missing = [];
  }
  return out;
}

// ---------------------------------------------------------------- API
app.get('/api/health', (req, res) => {
  res.json({ ok: true, rdkit: rdkitVersion(), uptime: Math.round(process.uptime()) });
});

app.post('/api/simulate', rateLimit(120), wrap(async (req, res) => {
  const { structure, phase = 'condensada' } = req.body || {};
  const { ctx, sim, explanation } = runSimulation(structure, PHASES.has(phase) ? phase : 'condensada');
  try {
    res.json({
      molecule: moleculeInfo(ctx),
      phase: sim.phase,
      spectrum: sim.spectrum,
      bands: sim.bands,
      groups: explanation.groups,
      explanation
    });
  } finally {
    ctx.dispose();
  }
}));

app.post('/api/depict', rateLimit(300), wrap(async (req, res) => {
  const { structure, atoms = [] } = req.body || {};
  const ctx = createContext(structure);
  try {
    res.json({ svg: ctx.svg(Array.isArray(atoms) ? atoms.filter(Number.isInteger) : []) });
  } finally {
    ctx.dispose();
  }
}));

app.get('/api/search', rateLimit(30), wrap(async (req, res) => {
  const q = String(req.query.q || '').slice(0, 300);
  const type = String(req.query.type || 'auto');
  const r = await searchSpectra(q, type);
  if (!r.compound) return res.status(404).json({ error: `No se encontró "${q}" en PubChem. Pruebe con el nombre en inglés, el número CAS o un SMILES.` });
  res.json(r);
}));

app.get('/api/resolve', rateLimit(60), wrap(async (req, res) => {
  const r = await resolveCompound(String(req.query.q || '').slice(0, 300));
  if (!r) return res.status(404).json({ error: 'Compuesto no encontrado en PubChem.' });
  res.json(r);
}));

app.post('/api/analyze', rateLimit(40), wrap(async (req, res) => {
  const { source = {}, structure, phase, yMode } = req.body || {};
  let raw;
  if (source.type === 'nist') {
    raw = await getNistSpectrum(String(source.id), Number(source.index) || 0);
    raw.meta.database = 'NIST Chemistry WebBook (SRD 69)';
    raw.meta.url = `https://webbook.nist.gov/cgi/cbook.cgi?ID=${source.id}&Type=IR-SPEC&Index=${Number(source.index) || 0}`;
  } else if (Array.isArray(source.x) && Array.isArray(source.y)) {
    raw = { x: source.x.map(Number), y: source.y.map(Number), xunits: source.xunits || '1/CM', yunits: source.yunits || '', meta: source.meta || {} };
  } else {
    return res.status(400).json({ error: 'Falta el espectro a analizar.' });
  }
  res.json(analyzeRaw(raw, { structure, phase, yMode }));
}));

const upload = multer({ storage: multer.memoryStorage(), limits: { fileSize: 5 * 1024 * 1024, files: 1 } });
app.post('/api/upload', rateLimit(20), upload.single('file'), wrap(async (req, res) => {
  if (!req.file) return res.status(400).json({ error: 'No se recibió ningún archivo.' });
  const text = req.file.buffer.toString('utf8');
  const raw = parseSpectrumFile(text, req.file.originalname);
  raw.meta = { ...raw.meta, database: 'Archivo del usuario', filename: req.file.originalname };
  res.json(analyzeRaw(raw, { structure: req.body.structure || null, phase: req.body.phase, yMode: req.body.yMode || 'auto' }));
}));

app.use('/api', (req, res) => res.status(404).json({ error: 'Ruta no encontrada.' }));

// ---------------------------------------------------------------- errores
app.use((err, req, res, _next) => {
  const status = err.status || (err.name === 'TimeoutError' ? 504 : 400);
  if (status >= 500 || !err.expose) console.error('[error]', req.method, req.path, err.message);
  const msg = err.name === 'TimeoutError'
    ? 'La base de datos externa tardó demasiado en responder. Inténtelo de nuevo.'
    : err.code === 'LIMIT_FILE_SIZE' ? 'El archivo supera 5 MB.' : err.message || 'Error inesperado.';
  res.status(status).json({ error: msg });
});

await initChem();
app.listen(PORT, () => {
  console.log(`IRAnalizador escuchando en http://0.0.0.0:${PORT} (RDKit ${rdkitVersion()})`);
});
