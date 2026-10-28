// Generación de explicaciones detalladas en español.
import { INT_LABEL } from './groups.js';

export const REGIONS = [
  { id: 'xh', name: 'Estiramientos X–H', lo: 2500, hi: 4000, text: 'Región de estiramientos de enlaces con hidrógeno (O–H, N–H, C–H). Como el H es muy ligero, la masa reducida μ es pequeña y la frecuencia ν = (1/2πc)·√(k/μ) es alta. La línea de 3000 cm⁻¹ separa los C–H de carbonos insaturados (sp², sp: por encima) de los saturados (sp³: por debajo).' },
  { id: 'triple', name: 'Triples enlaces y acumulados', lo: 2000, hi: 2500, text: 'Región de C≡N, C≡C y sistemas X=C=Y. Los triples enlaces tienen constantes de fuerza muy grandes. Es una zona casi vacía, por lo que cualquier banda aquí es muy diagnóstica (cuidado con el CO₂ atmosférico, ~2350 cm⁻¹).' },
  { id: 'double', name: 'Dobles enlaces', lo: 1500, hi: 2000, text: 'Región de C=O, C=C, C=N y flexiones N–H. El carbonilo (1870–1630 cm⁻¹) suele ser la banda más intensa del espectro porque el enlace es muy polar y su vibración produce un gran cambio de momento dipolar.' },
  { id: 'finger', name: 'Huella dactilar', lo: 400, hi: 1500, text: 'Región de huella dactilar: estiramientos de enlaces sencillos (C–O, C–N, C–C, C–X) y numerosas flexiones que se acoplan entre sí. El patrón es único para cada molécula; algunas bandas intensas (C–O, flexiones C–H fuera del plano, C–Cl) sí son asignables.' }
];

function ihdText(ihd) {
  if (!Number.isInteger(ihd)) return `IHD = ${ihd}: valor no entero, típico de especies cargadas o sales; interprételo con cautela.`;
  if (ihd === 0) return 'IHD = 0: molécula saturada y acíclica. No puede contener C=O, C=C, C≡N ni anillos; sus grupos funcionales sólo pueden ser O–H, N–H, C–O, C–N, C–X, etc.';
  let t = `IHD = ${ihd}: la molécula tiene ${ihd} insaturación${ihd > 1 ? 'es' : ''} (cada doble enlace o anillo cuenta 1; cada triple enlace, 2).`;
  if (ihd >= 4) t += ' Un IHD ≥ 4 es compatible con un anillo bencénico (1 anillo + 3 dobles enlaces).';
  return t;
}

function fmtBand(bd) {
  return `${bd.center} cm⁻¹ (${bd.range[0]}–${bd.range[1]}; ${INT_LABEL[bd.label] || bd.label})`;
}

/** Explicación de una estructura simulada. */
export function explainStructure(ctx, sim) {
  const { groups, bands } = sim;
  const main = groups.filter((g) => g.id !== 'alkyl');
  const summary = [];
  summary.push(`Fórmula molecular ${ctx.formula}, masa molar ${ctx.mw.toFixed(2)} g/mol. ${ihdText(ctx.ihd)}`);
  if (main.length) {
    summary.push(`Grupos funcionales reconocidos en la estructura: ${main.map((g) => g.name + (g.variants.length ? ` (${g.variants.join('; ')})` : '')).join(', ')}.`);
  } else {
    summary.push('La estructura sólo contiene enlaces C–C y C–H saturados (hidrocarburo saturado): el espectro estará dominado por los estiramientos C–H (3000–2850 cm⁻¹) y las flexiones CH₂/CH₃ (1470–1370 cm⁻¹).');
  }
  const strongest = bands.filter((bd) => !bd.overtone).sort((a, b) => b.height - a.height).slice(0, 3);
  if (strongest.length) {
    summary.push(`Bandas más intensas previstas: ${strongest.map((bd) => `${bd.center} cm⁻¹ (${bd.vib}, ${bd.groupName})`).join('; ')}.`);
  }
  summary.push(sim.phase === 'gas'
    ? 'Simulación en fase gaseosa: sin puentes de hidrógeno (O–H/N–H agudos y a mayor frecuencia) y C=O desplazados unos 15 cm⁻¹ hacia arriba.'
    : 'Simulación en fase condensada (líquido puro o sólido): se consideran los puentes de hidrógeno, que ensanchan y desplazan a menor frecuencia las bandas O–H y N–H.');

  const regions = REGIONS.map((r) => ({
    ...r,
    bands: bands.filter((bd) => bd.center >= r.lo && bd.center < r.hi).map((bd) => ({
      center: bd.center, range: bd.range, label: bd.label, labelText: bd.labelText, vib: bd.vib, groupName: bd.groupName, why: bd.why, height: bd.height, overtone: bd.overtone
    }))
  }));

  const ids = new Set(groups.map((g) => g.id));
  const absences = [];
  const hasOH = ['alcohol', 'phenol', 'carboxylic_acid'].some((i) => ids.has(i));
  const hasNH = bands.some((bd) => /N–H/.test(bd.vib) && bd.center > 3000);
  if (!hasOH && !hasNH) absences.push('No hay O–H ni N–H: la región 3600–3200 cm⁻¹ debe quedar libre de bandas anchas. Si aparece una banda ahí, sospeche humedad (agua) en la muestra.');
  if (!bands.some((bd) => /C=O/.test(bd.vib))) absences.push('No hay grupo carbonilo: no se espera ninguna banda intensa entre 1870 y 1630 cm⁻¹.');
  if (!bands.some((bd) => bd.center > 2000 && bd.center < 2300 && !bd.overtone)) absences.push('No hay triples enlaces ni sistemas acumulados: la región 2300–2000 cm⁻¹ debería estar vacía.');
  if (!ids.has('aromatic') && !ids.has('alkene') && !ids.has('heteroaromatic') && !ids.has('alkyne')) absences.push('No hay C–H sobre carbonos sp² o sp: no se esperan bandas C–H por encima de 3000 cm⁻¹.');

  const groupsOut = groups.map((g) => ({
    id: g.id,
    name: g.name,
    category: g.category,
    theory: g.theory,
    variants: g.variants,
    count: g.count,
    atoms: g.atoms,
    bands: g.bands.slice().sort((a, b) => b.center - a.center).map((bd) => ({ ...bd, text: fmtBand(bd) }))
  }));

  return { summary, regions, groups: groupsOut, absences };
}

/** Texto explicativo por pico experimental asignado. */
export function explainPeaks(assigned, { phase }) {
  return assigned.peaks.map((p) => {
    const head = `${p.x} cm⁻¹ — intensidad ${p.intensityText}, banda ${p.shape} (anchura ≈ ${p.width} cm⁻¹)`;
    let text;
    if (p.assignment) {
      const a = p.assignment;
      const dev = p.x - a.center;
      text = `${a.vib} — ${a.groupName}. Valor esperado ≈ ${a.center} cm⁻¹ (intervalo ${a.range[0]}–${a.range[1]}); desviación ${dev >= 0 ? '+' : ''}${dev} cm⁻¹. ${a.why}`;
      if (p.note) text += ' ' + p.note;
    } else if (p.correlation && p.correlation.length && p.x >= 1500) {
      text = `Posibles asignaciones según la tabla de correlación: ${p.correlation.map((c) => c.vib + (c.note ? ` (${c.note})` : '')).join('; ')}.` + (p.note ? ' ' + p.note : '');
    } else {
      text = p.note || 'Sin asignación clara.';
      if (p.correlation && p.correlation.length) text += ` Candidatos generales: ${p.correlation.map((c) => c.vib).join('; ')}.`;
    }
    if (p.satNote) text += ' ' + p.satNote;
    const region = REGIONS.find((r) => p.x >= r.lo && p.x < r.hi);
    return { x: p.x, head, text, region: region ? region.name : '', confidence: p.assignment ? p.assignment.score : 0 };
  });
}

export function explainMissing(missing, phase) {
  return missing.map((bd) => ({
    center: bd.center,
    text: `No se observa claramente la banda esperada a ${fmtBand(bd)}: ${bd.vib} (${bd.groupName}).` +
      (/O–H|N–H/.test(bd.vib) && phase === 'gas' ? ' En gas esta banda es aguda y débil; puede quedar por debajo del umbral de detección.' : ' Puede estar solapada con otra banda, ser más débil de lo previsto o indicar que la estructura propuesta no es correcta.')
  }));
}
