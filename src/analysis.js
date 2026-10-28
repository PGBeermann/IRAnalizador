// Procesamiento de espectros experimentales: lectura, normalización,
// detección de picos y asignación de bandas.
import { convert } from 'jcampconverter';
import { CORRELATION, INT_LABEL } from './groups.js';
import { intensityClass } from './simulator.js';

// ------------------------------------------------------------ LECTURA
function headerFields(text) {
  const out = {};
  const head = text.split(/##XYDATA|##PEAK TABLE|##XYPOINTS/i)[0];
  for (const m of head.matchAll(/^##\$?([^=]+)=(.*)$/gm)) {
    out[m[1].trim().toUpperCase()] = m[2].trim();
  }
  return out;
}

/** Lee un archivo JCAMP-DX (IR). */
export function parseJcamp(text) {
  const h = headerFields(text);
  if (/spectrum not found/i.test(h.TITLE || '')) throw new Error('Espectro no encontrado.');
  const res = convert(text, { keepRecordsRegExp: /.*/ });
  const entry = (res.flatten || []).find((e) => e.spectra && e.spectra.length && e.spectra[0].data);
  if (!entry) throw new Error('El archivo JCAMP-DX no contiene datos espectrales legibles.');
  const d = entry.spectra[0].data;
  return {
    x: Array.from(d.x),
    y: Array.from(d.y),
    xunits: h.XUNITS || '1/CM',
    yunits: h.YUNITS || '',
    meta: {
      title: h.TITLE || '',
      cas: h['CAS REGISTRY NO'] || '',
      formula: (h.MOLFORM || '').replace(/\s+/g, ''),
      state: h.STATE || '',
      origin: h.ORIGIN || h.OWNER || '',
      source: h['NIST SOURCE'] || h['SOURCE REFERENCE'] || '',
      instrument: h['SPECTROMETER/DATA SYSTEM'] || '',
      resolution: h.RESOLUTION || '',
      date: h.DATE || '',
      pathLength: h['PATH LENGTH'] || '',
      sampling: h['SAMPLING PROCEDURE'] || ''
    }
  };
}

/** Lee CSV/TXT de dos columnas (x, y). */
export function parseTable(text) {
  const x = [];
  const y = [];
  for (const line of text.split(/\r?\n/)) {
    const t = line.trim();
    if (!t || /^[#;%a-df-z"]/i.test(t)) continue;
    const parts = t.split(/[\s,;\t]+/).filter(Boolean);
    if (parts.length < 2) continue;
    const a = parseFloat(parts[0].replace(',', '.'));
    const b2 = parseFloat(parts[1].replace(',', '.'));
    if (Number.isFinite(a) && Number.isFinite(b2)) {
      x.push(a);
      y.push(b2);
    }
  }
  if (x.length < 20) throw new Error('No se encontraron suficientes pares de datos (x, y) en el archivo.');
  return { x, y, xunits: '1/CM', yunits: '', meta: { title: '', state: '' } };
}

export function parseSpectrumFile(text, filename = '') {
  if (/##\s*TITLE/i.test(text) || /\.(jdx|dx|jcamp)$/i.test(filename)) return parseJcamp(text);
  return parseTable(text);
}

export function phaseFromState(state = '') {
  return /gas|vapor|vapour/i.test(state) ? 'gas' : 'condensada';
}

// ------------------------------------------------------------ NORMALIZACIÓN
/**
 * Lleva el espectro a una malla uniforme (2 cm⁻¹, 400–4000) en absorbancia
 * normalizada (0–1) y %T equivalente.
 * yMode: 'auto' | 'T' | 'T%' | 'A'
 */
export function normalizeSpectrum(raw, yMode = 'auto') {
  let x = raw.x.slice();
  let y = raw.y.slice();
  const xu = (raw.xunits || '').toUpperCase();
  if (/MICRO/.test(xu)) x = x.map((v) => 1e4 / v);
  else if (/NANO/.test(xu)) x = x.map((v) => 1e7 / v);

  // ordenar por x creciente
  const idx = x.map((_, i) => i).filter((i) => Number.isFinite(x[i]) && Number.isFinite(y[i]));
  idx.sort((a, b) => x[a] - x[b]);
  x = idx.map((i) => x[i]);
  y = idx.map((i) => y[i]);
  if (x.length < 10) throw new Error('El espectro tiene muy pocos puntos.');

  const yu = (raw.yunits || '').toUpperCase();
  let mode = yMode;
  if (mode === 'auto') {
    if (/TRANSM/.test(yu)) mode = 'T';
    else if (/ABSORB|MOL|BASE 10/.test(yu)) mode = 'A';
    else {
      const sorted = y.slice().sort((a, b) => a - b);
      const med = sorted[Math.floor(sorted.length / 2)];
      const max = sorted[sorted.length - 1];
      mode = max > 1.5 && max <= 110 ? 'T%' : med > 0.5 && max <= 1.2 ? 'T' : 'A';
    }
  }
  if (mode === 'T') {
    const max = Math.max(...y);
    if (max > 1.5) mode = 'T%';
  }

  let A;
  if (mode === 'T' || mode === 'T%') {
    const f = mode === 'T%' ? 0.01 : 1;
    A = y.map((v) => -Math.log10(Math.min(1.5, Math.max(1e-4, v * f))));
  } else {
    A = y.slice();
  }

  const lo = Math.max(400, Math.ceil(x[0] / 2) * 2);
  const hi = Math.min(4000, Math.floor(x[x.length - 1] / 2) * 2);
  if (hi - lo < 200) throw new Error('El intervalo espectral es demasiado corto o no está en cm⁻¹.');
  const gx = [];
  const gA = [];
  let j = 0;
  for (let v = lo; v <= hi; v += 2) {
    while (j < x.length - 2 && x[j + 1] < v) j++;
    const x0 = x[j];
    const x1 = x[j + 1];
    const t = x1 === x0 ? 0 : (v - x0) / (x1 - x0);
    gx.push(v);
    gA.push(A[j] + (A[j + 1] - A[j]) * Math.min(1, Math.max(0, t)));
  }
  // línea base: percentil 3 % → 0 ; máximo → 1
  const sorted = gA.slice().sort((a, b) => a - b);
  const base = sorted[Math.floor(sorted.length * 0.03)];
  let max = -Infinity;
  for (const v of gA) if (v - base > max) max = v - base;
  const An = gA.map((v) => Math.round(Math.max(0, (v - base) / (max || 1)) * 1e4) / 1e4);
  const T = An.map((v) => Math.round(100 * Math.pow(10, -v) * 100) / 100);
  return { x: gx, A: An, T, yMode: mode, range: [lo, hi] };
}

// ------------------------------------------------------------ PICOS
function smooth(y, half = 3) {
  // Savitzky–Golay cuadrático de 7 puntos
  const c = [-2, 3, 6, 7, 6, 3, -2];
  const norm = 21;
  const out = y.slice();
  for (let i = half; i < y.length - half; i++) {
    let s = 0;
    for (let k = -half; k <= half; k++) s += c[k + half] * y[i + k];
    out[i] = s / norm;
  }
  return out;
}

export function pickPeaks(spec, opts = {}) {
  const { x } = spec;
  const y = smooth(spec.A);
  const minProm = opts.minProminence ?? 0.035;
  const minHeight = opts.minHeight ?? 0.04;
  const win = 3;
  const peaks = [];
  for (let i = win; i < y.length - win; i++) {
    let isMax = true;
    for (let k = 1; k <= win; k++) {
      if (y[i] < y[i - k] || y[i] < y[i + k]) { isMax = false; break; }
    }
    if (!isMax || y[i] < minHeight) continue;
    if (y[i] === y[i - 1]) continue;
    // prominencia
    let l = i;
    let lmin = y[i];
    while (l > 0 && y[l - 1] <= y[i]) { l--; if (y[l] < lmin) lmin = y[l]; }
    let r = i;
    let rmin = y[i];
    while (r < y.length - 1 && y[r + 1] <= y[i]) { r++; if (y[r] < rmin) rmin = y[r]; }
    const prom = y[i] - Math.max(lmin, rmin);
    if (prom < minProm) continue;
    // anchura a media prominencia
    const half = y[i] - prom / 2;
    let a = i;
    while (a > 0 && y[a] > half) a--;
    let bb = i;
    while (bb < y.length - 1 && y[bb] > half) bb++;
    const width = x[bb] - x[a];
    peaks.push({ x: x[i], height: Math.round(spec.A[i] * 1000) / 1000, prominence: Math.round(prom * 1000) / 1000, width });
  }
  // fusionar picos de una banda saturada (meseta con transmitancia ≈ 0)
  peaks.sort((p, q) => p.x - q.x);
  for (let k = peaks.length - 1; k > 0; k--) {
    const a = peaks[k - 1];
    const c = peaks[k];
    if (c.x - a.x <= 30 && a.height >= 0.95 && c.height >= 0.95) {
      a.x = Math.round((a.x + c.x) / 4) * 2;
      a.width = Math.max(a.width, c.width) + (c.x - a.x);
      a.prominence = Math.max(a.prominence, c.prominence);
      a.saturated = true;
      peaks.splice(k, 1);
    }
  }
  // limitar a los más prominentes
  peaks.sort((p, q) => q.prominence - p.prominence);
  const kept = peaks.slice(0, opts.max ?? 45);
  kept.sort((p, q) => q.x - p.x);
  for (const p of kept) {
    p.intensity = intensityClass(p.height);
    p.intensityText = INT_LABEL[p.intensity];
    p.shape = p.width > 150 ? 'ancha' : p.width > 45 ? 'media' : 'aguda';
  }
  return kept;
}

// ------------------------------------------------------------ ASIGNACIÓN
function scorePair(p, bd, phase) {
  const [lo, hi] = bd.range;
  const tol = phase === 'gas' ? 25 : 15;
  if (p.x < lo - tol || p.x > hi + tol) return 0;
  const sigma = Math.max(14, (hi - lo) / 2);
  const sPos = Math.exp(-0.5 * ((p.x - bd.center) / sigma) ** 2);
  const he = Math.min(1, bd.height ?? 0.5);
  const sInt = Math.max(0.25, 1 - 0.7 * Math.abs(he - p.height));
  let sShape = 1;
  if (bd.fwhm > 150 && p.width < 50) sShape = 0.55;
  if (bd.fwhm < 40 && p.width > 180) sShape = 0.6;
  const sOver = bd.overtone ? 0.6 : 1;
  return sPos * sInt * sShape * sOver;
}

function correlationFor(p) {
  return CORRELATION
    .filter((r) => p.x >= Math.min(r.lo, r.hi) - 10 && p.x <= Math.max(r.lo, r.hi) + 10)
    .map((r) => {
      let s = 1;
      if (/ancha/.test(r.shape) && p.shape === 'aguda') s *= 0.5;
      if (r.shape === 'aguda' && p.shape === 'ancha') s *= 0.5;
      if ((r.label === 'mf' || r.label === 'f') && p.intensity === 'md') s *= 0.5;
      if ((r.label === 'md' || r.label === 'd') && p.height > 0.3) s *= 0.2;
      s /= 1 + Math.abs(r.hi - r.lo) / 400; // preferir intervalos específicos
      return { ...r, score: s };
    })
    .sort((a, b) => b.score - a.score)
    .slice(0, 3);
}

/**
 * Asigna picos experimentales a bandas esperadas (si hay estructura)
 * o a la tabla de correlación general.
 */
export function assignPeaks(peaks, bands, phase = 'condensada', range = [400, 4000]) {
  const result = peaks.map((p) => ({ ...p, assignment: null, alternatives: [], note: '' }));
  if (bands && bands.length) {
    const pairs = [];
    result.forEach((p, i) => bands.forEach((bd, j) => {
      const s = scorePair(p, bd, phase);
      if (s > 0.05) pairs.push({ i, j, s });
    }));
    pairs.sort((a, b) => b.s - a.s);
    const usedBand = new Map();
    for (const { i, j, s } of pairs) {
      const p = result[i];
      if (p.assignment) continue;
      const already = usedBand.get(j);
      // una banda puede explicar un segundo pico (desdoblamiento, ramas P/R) si la puntuación es buena
      if (already !== undefined && (s < 0.35 || already >= 2)) continue;
      p.assignment = { ...bands[j], score: Math.round(s * 100) / 100, shared: already !== undefined };
      usedBand.set(j, (already || 0) + 1);
    }
    for (const p of result) {
      const alts = pairs.filter((q) => result[q.i] === p && (!p.assignment || bands[q.j].id !== p.assignment.id)).slice(0, 2);
      p.alternatives = alts.map((q) => ({ center: bands[q.j].center, vib: bands[q.j].vib, groupName: bands[q.j].groupName, score: Math.round(q.s * 100) / 100 }));
    }
  }
  // picos sin asignar
  const strong = result.filter((p) => p.height > 0.5);
  for (const p of result) {
    if (p.saturated) p.satNote = 'Banda saturada (transmitancia ≈ 0 %): la muestra era demasiado concentrada en esta zona y la posición exacta del máximo es aproximada.';
    if (p.assignment) {
      if (p.assignment.shared) p.note = phase === 'gas'
        ? 'Componente adicional de la misma vibración: en fase gaseosa las bandas muestran ramas rotovibracionales P y R a ambos lados del centro.'
        : 'Componente adicional de la misma vibración (desdoblamiento por acoplamiento, resonancia de Fermi o diferentes conformaciones).';
      continue;
    }
    if ((p.x >= 2300 && p.x <= 2380) || (p.x >= 660 && p.x <= 675 && p.shape === 'aguda' && p.height < 0.3)) {
      p.note = 'Probable CO₂ atmosférico no compensado en el espectro de fondo (artefacto, no pertenece a la muestra).';
      p.artifact = true;
      continue;
    }
    const ot = strong.find((s) => s !== p && Math.abs(2 * s.x - p.x) < 40);
    if (ot && p.height < 0.2) {
      p.note = `Posible sobretono (2 × ${ot.x} cm⁻¹) o banda de combinación: aparece cerca del doble de la frecuencia de una banda fundamental intensa, con intensidad baja.`;
      continue;
    }
    if (phase === 'gas') {
      const near = result.find((q) => q !== p && q.assignment && Math.abs(q.x - p.x) <= 40 && q.height > p.height);
      if (near) {
        p.note = `Rama rotovibracional (P/R) de la banda a ${near.x} cm⁻¹: en gas, la rotación molecular desdobla cada vibración en ramas a ambos lados del centro.`;
        p.related = near.x;
        continue;
      }
    }
    p.correlation = correlationFor(p);
    if (p.x < 1500) {
      p.note = 'Banda de la región de huella dactilar (1500–400 cm⁻¹): corresponde a modos esqueléticos acoplados (C–C, C–O, C–N y flexiones) propios de la molécula concreta; no se asigna a un solo enlace, pero sirve para identificar el compuesto por comparación.';
    } else if (bands && bands.length) {
      p.note = 'No coincide con ninguna banda esperada para la estructura. Puede deberse a una impureza (p. ej., agua ~3400 y ~1640 cm⁻¹), a una banda de combinación o a que la estructura no corresponde al espectro.';
    }
  }
  // bandas esperadas no observadas
  const assignedIds = new Set(result.filter((p) => p.assignment).map((p) => p.assignment.id));
  const missing = (bands || []).filter((bd) => !assignedIds.has(bd.id) && !bd.overtone && (bd.height ?? 0) >= 0.2 &&
    bd.center >= range[0] + 10 && bd.center <= range[1] - 10);
  return { peaks: result, missing };
}

/** Similitud coseno entre dos espectros de absorbancia (con suavizado para tolerar pequeños desplazamientos). */
export function similarity(simSpec, expSpec) {
  const m = new Map(simSpec.x.map((v, i) => [v, simSpec.A[i]]));
  const a = [];
  const b2 = [];
  for (let i = 0; i < expSpec.x.length; i++) {
    const v = m.get(expSpec.x[i]);
    if (v === undefined) continue;
    a.push(v);
    b2.push(expSpec.A[i]);
  }
  const blur = (arr) => {
    const out = new Array(arr.length).fill(0);
    const w = 10;
    for (let i = 0; i < arr.length; i++) {
      let s = 0;
      let n = 0;
      for (let k = -w; k <= w; k++) {
        const j = i + k;
        if (j < 0 || j >= arr.length) continue;
        const g = Math.exp(-0.5 * (k / 5) ** 2);
        s += arr[j] * g;
        n += g;
      }
      out[i] = s / n;
    }
    return out;
  };
  const A1 = blur(a);
  const B1 = blur(b2);
  let dot = 0;
  let na = 0;
  let nb = 0;
  for (let i = 0; i < A1.length; i++) {
    dot += A1[i] * B1[i];
    na += A1[i] ** 2;
    nb += B1[i] ** 2;
  }
  return na && nb ? Math.round((dot / Math.sqrt(na * nb)) * 1000) / 10 : 0;
}

/** Hallazgos generales (sin estructura): qué grupos sugiere el espectro. */
export function diagnose(peaks) {
  const f = [];
  const inR = (lo, hi, minH = 0, shape) => peaks.filter((p) => p.x >= lo && p.x <= hi && p.height >= minH && (!shape || shape.includes(p.shape)));
  const co = inR(1650, 1870, 0.45);
  const chSp3 = inR(2840, 2999, 0.08);
  const chSp2 = inR(3001, 3110, 0.03);
  const ohBroad = inR(3200, 3600, 0.12, ['ancha', 'media']).filter((p) => p.width > 100);
  const ohFree = inR(3580, 3680, 0.05, ['aguda', 'media']);
  const acid = inR(2500, 3300, 0.1, ['ancha']).filter((p) => p.width > 250);
  const nh = inR(3300, 3500, 0.08, ['aguda', 'media']);
  const ald = inR(2690, 2780, 0.05);
  const triple = inR(2100, 2280, 0.08);
  const alkyneCH = inR(3270, 3330, 0.2, ['aguda']);
  const arom = inR(1580, 1610, 0.08).length && inR(1480, 1520, 0.08).length;
  const nitroA = inR(1500, 1570, 0.5);
  const nitroS = inR(1330, 1390, 0.4);
  const coStrong = inR(1000, 1300, 0.5);

  if (chSp3.length) f.push({ ok: true, text: `Bandas C–H por debajo de 3000 cm⁻¹ (${chSp3.map((p) => p.x).join(', ')}): hay carbonos saturados (sp³).` });
  if (chSp2.length) f.push({ ok: true, text: `Bandas C–H por encima de 3000 cm⁻¹ (${chSp2.map((p) => p.x).join(', ')}): hay C–H sobre carbonos insaturados (aromáticos o alquenos).` });
  if (acid.length && co.some((p) => p.x <= 1725)) f.push({ ok: true, text: 'Banda O–H muy ancha (3300–2500 cm⁻¹) junto con C=O ~1700 cm⁻¹: patrón típico de ÁCIDO CARBOXÍLICO.' });
  else if (ohBroad.length) f.push({ ok: true, text: `Banda ancha a ${ohBroad[0].x} cm⁻¹: O–H asociado por puentes de hidrógeno (alcohol o fenol) o N–H.` });
  if (ohFree.length) f.push({ ok: true, text: `Banda aguda a ${ohFree[0].x} cm⁻¹: O–H libre (alcohol/fenol en gas o disolución diluida).` });
  if (!ohBroad.length && !acid.length && !ohFree.length) f.push({ ok: false, text: 'No hay banda ancha entre 3600 y 3200 cm⁻¹: probablemente no hay grupos O–H.' });
  if (nh.length >= 2 && !ohBroad.length) f.push({ ok: true, text: `Dos bandas entre 3500 y 3300 cm⁻¹ (${nh.map((p) => p.x).join(', ')}): posible NH₂ (amina o amida primaria).` });
  if (co.length) f.push({ ok: true, text: `Banda intensa a ${co.map((p) => p.x).join(', ')} cm⁻¹: hay un grupo CARBONILO (C=O).` + (co.length >= 2 && Math.abs(co[0].x - co[1].x) > 40 && Math.abs(co[0].x - co[1].x) < 80 && co[0].x > 1740 ? ' Dos bandas separadas ~60 cm⁻¹ por encima de 1740: posible anhídrido.' : '') });
  else f.push({ ok: false, text: 'No hay banda intensa entre 1870 y 1650 cm⁻¹: no hay grupo carbonilo.' });
  if (ald.length && co.length) f.push({ ok: true, text: `Banda a ${ald[0].x} cm⁻¹ junto al C=O: doblete de Fermi típico de ALDEHÍDO.` });
  if (co.some((p) => p.x >= 1730 && p.x <= 1760) && coStrong.length) f.push({ ok: true, text: 'C=O a 1760–1730 cm⁻¹ con bandas C–O intensas entre 1300 y 1000 cm⁻¹: compatible con ÉSTER.' });
  if (alkyneCH.length) f.push({ ok: true, text: `Banda aguda e intensa a ${alkyneCH[0].x} cm⁻¹: ≡C–H de alquino terminal.` });
  if (triple.length) f.push({ ok: true, text: `Banda en la región de triples enlaces (${triple.map((p) => p.x).join(', ')} cm⁻¹): C≡N, C≡C o sistema acumulado (N=C=O, N₃).` });
  if (arom) f.push({ ok: true, text: 'Bandas cerca de 1600 y 1500 cm⁻¹: compatibles con un anillo AROMÁTICO (confirmar con C–H > 3000 y flexiones 900–675 cm⁻¹).' });
  if (nitroA.length && nitroS.length) f.push({ ok: true, text: `Dos bandas muy intensas a ${nitroA[0].x} y ${nitroS[0].x} cm⁻¹: compatibles con un grupo NITRO.` });
  return f;
}
