// Simulación de espectros IR por frecuencias de grupo.
import { RULES, INT_LABEL } from './groups.js';

export const GRID = { min: 400, max: 4000, step: 2 };

export function makeGrid(min = GRID.min, max = GRID.max, step = GRID.step) {
  const x = [];
  for (let v = min; v <= max + 1e-9; v += step) x.push(Math.round(v * 100) / 100);
  return x;
}

function applyPhase(band, phase) {
  if (phase !== 'gas') return { ...band };
  const g = band.gas;
  if (!g) return { ...band, fwhm: band.fwhm > 100 ? band.fwhm : band.fwhm * 0.85 };
  const delta = (g.center ?? band.center) - band.center;
  return {
    ...band,
    center: g.center ?? band.center,
    fwhm: g.fwhm ?? band.fwhm * 0.85,
    eps: g.eps ?? band.eps,
    vib: g.vib ?? band.vib,
    range: [band.range[0] + delta, band.range[1] + delta],
    why: band.why + ' ' + (g.why || `En fase gaseosa esta banda aparece cerca de ${g.center} cm⁻¹ (sin interacciones intermoleculares).`)
  };
}

/** Perfil pseudo-Voigt normalizado a 1 en el centro. */
function profile(x, c, fwhm) {
  const t = (x - c) / (fwhm / 2);
  const t2 = t * t;
  return 0.5 / (1 + t2) + 0.5 * Math.exp(-Math.LN2 * t2);
}

/** Detecta grupos funcionales y devuelve grupos con bandas agregadas. */
export function detectGroups(ctx, phase = 'condensada') {
  const groups = [];
  let bandId = 0;
  for (const rule of RULES) {
    let inst = [];
    try {
      inst = rule.detect(ctx) || [];
    } catch (e) {
      console.warn(`[regla ${rule.id}]`, e.message);
    }
    if (!inst.length) continue;
    const map = new Map();
    const atoms = new Set();
    for (const it of inst) {
      it.atoms.forEach((a) => atoms.add(a));
      for (const raw of it.bands) {
        const bd = applyPhase(raw, phase);
        if (bd.center < GRID.min - 50 || bd.center > GRID.max + 50) continue;
        const key = `${bd.vib}|${Math.round(bd.center)}`;
        const prev = map.get(key);
        if (prev) {
          prev.eps += bd.eps;
          prev.n += 1;
          it.atoms.forEach((a) => prev.atoms.add(a));
        } else {
          map.set(key, { ...bd, n: 1, atoms: new Set(it.atoms) });
        }
      }
    }
    const variants = [...new Set(inst.map((i) => i.variant).filter(Boolean))];
    const count = inst.reduce((s, i) => s + (i.count || 1), 0);
    groups.push({
      id: rule.id,
      name: rule.name,
      category: rule.category,
      theory: rule.theory,
      count,
      variants,
      atoms: [...atoms],
      bands: [...map.values()].map((bd) => ({
        id: `b${++bandId}`,
        groupId: rule.id,
        groupName: rule.name,
        center: Math.round(bd.center),
        fwhm: bd.fwhm,
        eps: bd.eps,
        label: bd.label,
        labelText: INT_LABEL[bd.label] || bd.label,
        vib: bd.vib,
        why: bd.why,
        range: [Math.round(Math.min(...bd.range)), Math.round(Math.max(...bd.range))],
        n: bd.n,
        overtone: bd.overtone,
        atoms: [...bd.atoms]
      }))
    });
  }
  return groups;
}

/** Construye el espectro simulado (absorbancia normalizada y %T). */
export function buildSpectrum(bands, x = makeGrid()) {
  const A = new Float64Array(x.length);
  for (const bd of bands) {
    const reach = bd.fwhm * 6;
    for (let i = 0; i < x.length; i++) {
      const d = x[i] - bd.center;
      if (d < -reach || d > reach) continue;
      A[i] += bd.eps * profile(x[i], bd.center, bd.fwhm);
    }
  }
  let max = 0;
  for (const v of A) if (v > max) max = v;
  const scale = max > 0 ? 1 / max : 1;
  const An = Array.from(A, (v) => Math.round(v * scale * 1e4) / 1e4);
  const T = An.map((v) => Math.round(100 * Math.pow(10, -v) * 100) / 100);
  // altura relativa de cada banda aislada (para clasificar intensidades)
  for (const bd of bands) bd.height = Math.round(bd.eps * scale * 1000) / 1000;
  return { x, A: An, T };
}

export function intensityClass(h) {
  if (h >= 0.65) return 'f';
  if (h >= 0.3) return 'm';
  if (h >= 0.1) return 'd';
  return 'md';
}

/** Simulación completa a partir de un contexto químico. */
export function simulate(ctx, phase = 'condensada') {
  const groups = detectGroups(ctx, phase);
  const bands = groups.flatMap((g) => g.bands);
  const spectrum = buildSpectrum(bands);
  bands.sort((a, b) => b.center - a.center);
  return { groups, bands, spectrum, phase };
}
