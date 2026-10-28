// Acceso a bases de datos externas: PubChem (identificación) y NIST WebBook (espectros IR).
import { parseJcamp } from './analysis.js';

const UA = 'IRAnalizador/1.0 (academic IR spectroscopy app)';
const PUBCHEM = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';
const NIST = 'https://webbook.nist.gov/cgi/cbook.cgi';

// ---------------------------------------------------------- caché simple con caducidad
class TTLCache {
  constructor(max = 300, ttlMs = 12 * 3600e3) { this.max = max; this.ttl = ttlMs; this.map = new Map(); }
  get(k) {
    const e = this.map.get(k);
    if (!e) return undefined;
    if (Date.now() > e.t) { this.map.delete(k); return undefined; }
    this.map.delete(k); this.map.set(k, e);
    return e.v;
  }
  set(k, v) {
    this.map.set(k, { v, t: Date.now() + this.ttl });
    if (this.map.size > this.max) this.map.delete(this.map.keys().next().value);
  }
}
const cache = new TTLCache();

async function fetchText(url, opts = {}, timeout = 20000) {
  const res = await fetch(url, { ...opts, headers: { 'User-Agent': UA, ...(opts.headers || {}) }, signal: AbortSignal.timeout(timeout) });
  const text = await res.text();
  return { ok: res.ok, status: res.status, text };
}

// Cola para NIST: una petición cada NIST_GAP ms como mínimo (NIST limita la frecuencia).
const NIST_GAP = Number(process.env.NIST_MIN_INTERVAL_MS) || 1500;
let nistChain = Promise.resolve();
let nistLast = 0;
function nistFetch(url) {
  const run = nistChain.then(async () => {
    const wait = nistLast + NIST_GAP - Date.now();
    if (wait > 0) await new Promise((r) => setTimeout(r, wait));
    try {
      return await fetchText(url);
    } finally {
      nistLast = Date.now();
    }
  });
  nistChain = run.catch(() => {});
  return run;
}

function decodeHtml(s) {
  return s.replace(/<[^>]+>/g, '').replace(/&amp;/g, '&').replace(/&lt;/g, '<').replace(/&gt;/g, '>')
    .replace(/&#(\d+);/g, (_, n) => String.fromCharCode(Number(n))).replace(/&nbsp;/g, ' ').replace(/\s+/g, ' ').trim();
}

// ---------------------------------------------------------- PubChem
const PROPS = 'Title,IUPACName,MolecularFormula,MolecularWeight,SMILES,ConnectivitySMILES,InChI,InChIKey';

export function classifyQuery(q) {
  const s = q.trim();
  if (/^\d{2,7}-\d{2}-\d$/.test(s)) return 'cas';
  if (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(s)) return 'inchikey';
  if (/^InChI=/.test(s)) return 'inchi';
  if (/^\d+$/.test(s)) return 'cid';
  if (/[=#()[\]@+\\/]/.test(s) && !/\s/.test(s)) return 'smiles';
  if (/^(?:Cl|Br|[BCNOPSFIcnops]|\d)+$/.test(s)) return 'smiles';
  return 'name';
}

async function pubchemProps(ns, value) {
  const url = ns === 'cid'
    ? `${PUBCHEM}/compound/cid/${encodeURIComponent(value)}/property/${PROPS}/JSON`
    : `${PUBCHEM}/compound/${ns}/property/${PROPS}/JSON`;
  const opts = ns === 'cid' ? {} : {
    method: 'POST',
    headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
    body: new URLSearchParams({ [ns]: value }).toString()
  };
  const r = await fetchText(url, opts);
  if (!r.ok) return null;
  const j = JSON.parse(r.text);
  const p = j?.PropertyTable?.Properties?.[0];
  return p || null;
}

async function pubchemCAS(cid) {
  const r = await fetchText(`${PUBCHEM}/compound/cid/${cid}/synonyms/JSON`);
  if (!r.ok) return { cas: [], synonyms: [] };
  const syn = JSON.parse(r.text)?.InformationList?.Information?.[0]?.Synonym || [];
  const cas = syn.filter((s) => /^\d{2,7}-\d{2}-\d$/.test(s) && validCAS(s));
  return { cas: [...new Set(cas)].slice(0, 5), synonyms: syn.slice(0, 12) };
}

function validCAS(cas) {
  const digits = cas.replace(/-/g, '');
  const check = Number(digits.at(-1));
  const body = digits.slice(0, -1).split('').reverse();
  const sum = body.reduce((s, d, i) => s + Number(d) * (i + 1), 0);
  return sum % 10 === check;
}

/** Resuelve un nombre, CAS, SMILES, InChI/InChIKey o CID en PubChem. */
export async function resolveCompound(query, type = 'auto') {
  const q = String(query || '').trim();
  if (!q) throw new Error('Escriba un nombre, número CAS o SMILES.');
  const key = `pc:${type}:${q}`;
  const hit = cache.get(key);
  if (hit) return hit;
  let kind = type === 'auto' ? classifyQuery(q) : type;
  const order = kind === 'cas' ? ['name'] : kind === 'smiles' ? ['smiles', 'name'] : kind === 'name' ? ['name', 'smiles'] : [kind];
  let props = null;
  for (const ns of order) {
    try { props = await pubchemProps(ns, q); } catch { props = null; }
    if (props) { kind = ns; break; }
  }
  if (!props) return null;
  const { cas, synonyms } = await pubchemCAS(props.CID).catch(() => ({ cas: [], synonyms: [] }));
  if (kind === 'name' && /^\d{2,7}-\d{2}-\d$/.test(q) && !cas.includes(q)) cas.unshift(q);
  const out = {
    cid: props.CID,
    name: props.Title || props.IUPACName || q,
    iupac: props.IUPACName || '',
    formula: props.MolecularFormula || '',
    mw: Number(props.MolecularWeight) || null,
    smiles: props.SMILES || props.IsomericSMILES || props.ConnectivitySMILES || props.CanonicalSMILES || '',
    inchi: props.InChI || '',
    inchikey: props.InChIKey || '',
    cas,
    synonyms,
    pubchemUrl: `https://pubchem.ncbi.nlm.nih.gov/compound/${props.CID}`
  };
  cache.set(key, out);
  return out;
}

// ---------------------------------------------------------- NIST WebBook
export function nistIdFromCAS(cas) {
  return 'C' + cas.replace(/-/g, '');
}

async function nistJcamp(id, index) {
  const key = `nist:${id}:${index}`;
  const hit = cache.get(key);
  if (hit !== undefined) return hit;
  const r = await nistFetch(`${NIST}?JCAMP=${encodeURIComponent(id)}&Index=${index}&Type=IR`);
  if (!r.ok) {
    // error temporal (p. ej. limitación de NIST): no se guarda en caché
    const e = new Error(r.status === 429
      ? 'NIST WebBook está limitando las consultas en este momento. Espere un minuto e inténtelo de nuevo.'
      : `NIST respondió con el código ${r.status}. Inténtelo de nuevo en unos segundos.`);
    e.status = 502;
    throw e;
  }
  const text = /##TITLE/.test(r.text) && !/Spectrum not found/i.test(r.text) ? r.text : null;
  cache.set(key, text);
  return text;
}

/**
 * Lista los espectros IR disponibles en NIST con UNA sola petición
 * (página del compuesto con Mask=80 = datos IR). El JCAMP se descarga
 * sólo cuando el usuario elige un espectro.
 */
export async function listNistSpectra(id) {
  const key = `nistlist:${id}`;
  const hit = cache.get(key);
  if (hit) return hit;
  const r = await nistFetch(`${NIST}?ID=${encodeURIComponent(id)}&Mask=80&Units=SI`);
  if (!r.ok) {
    const e = new Error(r.status === 429
      ? 'NIST WebBook está limitando las consultas en este momento. Espere un minuto e inténtelo de nuevo.'
      : `NIST respondió con el código ${r.status}.`);
    e.status = 502;
    throw e;
  }
  const html = r.text;
  const list = [];
  if (/id="IR-SPEC"|Infrared Spectrum/i.test(html) && /JCAMP=|IR-SPEC/i.test(html)) {
    const seen = new Set();
    // enlaces a otros espectros: ...Type=IR-SPEC&amp;Index=N...>texto</a>
    for (const m of html.matchAll(/<a[^>]+Type=IR-SPEC&(?:amp;)?Index=(\d+)[^>]*>([\s\S]*?)<\/a>/gi)) {
      const index = Number(m[1]);
      const label = decodeHtml(m[2]);
      if (seen.has(index) || !label || /^(top|infrared spectrum|ir spectrum)$/i.test(label)) continue;
      seen.add(index);
      list.push({ id, index, title: label, state: label });
    }
    // espectro mostrado en la propia página (si no figura ya)
    const shown = html.match(/JCAMP=[A-Z]\d+&(?:amp;)?Index=(\d+)/);
    if (shown && !seen.has(Number(shown[1]))) {
      const h3 = [...html.matchAll(/<h3>([^<]+)<\/h3>/gi)].map((x) => decodeHtml(x[1])).find((t) => /spectrum/i.test(t));
      list.unshift({ id, index: Number(shown[1]), title: h3 || 'Espectro IR', state: h3 || '' });
    }
    list.sort((a, b) => a.index - b.index);
    // Los espectros cuantitativos (NIST/EPA, fase gas) aparecen como una tabla de resoluciones:
    // se conserva sólo el de resolución más cercana a 1 cm⁻¹.
    const quant = list.filter((s) => /^\d+(\.\d+)?$/.test(s.title));
    if (quant.length) {
      const best = quant.reduce((a, b) => (Math.abs(Number(b.title) - 1) < Math.abs(Number(a.title) - 1) ? b : a));
      best.title = `Gas · espectro cuantitativo NIST/EPA, resolución ${best.title} cm⁻¹`;
      best.state = 'gas';
      for (const q of quant) if (q !== best) list.splice(list.indexOf(q), 1);
    }
    for (const s of list) {
      if (/^gas$/i.test(s.title)) s.title = 'Fase gas';
      s.url = `https://webbook.nist.gov/cgi/cbook.cgi?ID=${id}&Type=IR-SPEC&Index=${s.index}`;
      s.phase = /gas|vapor/i.test(s.state) ? 'gas' : /solution|liquid|solid|film|mull|kbr|nujol|neat/i.test(s.state) ? 'condensada' : '';
    }
  }
  cache.set(key, list);
  return list;
}

/** Busca el ID de NIST a partir del InChI (si no hay CAS útil). */
async function nistIdFromInchi(inchi) {
  if (!inchi) return null;
  const r = await nistFetch(`${NIST}?InChI=${encodeURIComponent(inchi)}&Units=SI`);
  if (!r.ok) return null;
  const ids = [...r.text.matchAll(/[?&]ID=([A-Z]\d{3,})/g)].map((m) => m[1]);
  if (!ids.length) return null;
  const freq = {};
  ids.forEach((i) => { freq[i] = (freq[i] || 0) + 1; });
  return Object.entries(freq).sort((a, b) => b[1] - a[1])[0][0];
}

/** Búsqueda completa: compuesto (PubChem) + espectros IR (NIST). */
export async function searchSpectra(query, type = 'auto') {
  const compound = await resolveCompound(query, type);
  if (!compound) return { compound: null, spectra: [], nistId: null };
  let spectra = [];
  let nistId = null;
  let nistError = null;
  for (const cas of compound.cas.slice(0, 3)) {
    const id = nistIdFromCAS(cas);
    try {
      spectra = await listNistSpectra(id);
    } catch (e) {
      nistError = e.message;
      continue;
    }
    if (spectra.length) { nistId = id; break; }
  }
  if (!spectra.length && !nistError) {
    const id = await nistIdFromInchi(compound.inchi).catch(() => null);
    if (id) {
      spectra = await listNistSpectra(id).catch((e) => { nistError = e.message; return []; });
      if (spectra.length) nistId = id;
    }
  }
  return {
    compound,
    spectra,
    nistId,
    nistError: spectra.length ? null : nistError,
    nistUrl: nistId ? `https://webbook.nist.gov/cgi/cbook.cgi?ID=${nistId}&Type=IR-SPEC&Index=0` : null
  };
}

/** Descarga y lee un espectro concreto de NIST. */
export async function getNistSpectrum(id, index = 0) {
  if (!/^[A-Z]\d{3,}$/.test(id)) throw new Error('Identificador NIST no válido.');
  const text = await nistJcamp(id, Number(index) || 0);
  if (!text) throw new Error('NIST no tiene ese espectro.');
  return parseJcamp(text);
}
