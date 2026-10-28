// Envoltorio de RDKit (WebAssembly) para analizar estructuras químicas.
import initRDKitModule from '@rdkit/rdkit';

let RDKit = null;

export async function initChem() {
  if (!RDKit) RDKit = await initRDKitModule();
  return RDKit;
}

export function rdkitVersion() {
  return RDKit ? RDKit.version() : null;
}

const SYMBOLS = {
  1: 'H', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 11: 'Na', 12: 'Mg', 14: 'Si', 15: 'P',
  16: 'S', 17: 'Cl', 19: 'K', 34: 'Se', 35: 'Br', 53: 'I', 3: 'Li', 13: 'Al', 50: 'Sn'
};

/** Cache de consultas SMARTS compiladas (no se liberan: son pocas y se reutilizan). */
const qCache = new Map();
function qmol(smarts) {
  let q = qCache.get(smarts);
  if (!q) {
    q = RDKit.get_qmol(smarts);
    if (!q || !q.is_valid()) throw new Error(`SMARTS inválido: ${smarts}`);
    qCache.set(smarts, q);
  }
  return q;
}

/**
 * Crea un contexto de análisis para una molécula.
 * input: SMILES o bloque MOL. Llamar a ctx.dispose() al terminar.
 */
export function createContext(input) {
  if (!RDKit) throw new Error('RDKit no inicializado');
  const text = String(input || '').trim();
  if (!text) throw new Error('No se recibió ninguna estructura.');
  const mol = RDKit.get_mol(text);
  if (!mol || !mol.is_valid()) {
    if (mol) mol.delete();
    throw new Error('La estructura no es válida (revise valencias y enlaces).');
  }

  const json = JSON.parse(mol.get_json());
  const defaults = json.defaults.atom;
  const m0 = json.molecules[0];
  const atoms = m0.atoms.map((a) => ({
    z: a.z ?? defaults.z,
    h: a.impHs ?? defaults.impHs,
    chg: a.chg ?? defaults.chg
  }));
  const bondDefaults = json.defaults.bond;
  const bonds = (m0.bonds || []).map((b) => ({
    atoms: b.atoms,
    bo: b.bo ?? bondDefaults.bo,
    stereo: b.stereo ?? bondDefaults.stereo
  }));
  if (atoms.length === 0) {
    mol.delete();
    throw new Error('La estructura está vacía. Dibuje una molécula.');
  }

  // Fórmula molecular (orden de Hill)
  const counts = {};
  for (const a of atoms) {
    const s = SYMBOLS[a.z] || `Z${a.z}`;
    counts[s] = (counts[s] || 0) + 1;
    if (a.h) counts.H = (counts.H || 0) + a.h;
  }
  const charge = atoms.reduce((s, a) => s + a.chg, 0);
  const order = Object.keys(counts).sort((a, b) => {
    const rank = (x) => (x === 'C' ? 0 : x === 'H' ? 1 : 2);
    return rank(a) - rank(b) || a.localeCompare(b);
  });
  const formula = order.map((s) => s + (counts[s] > 1 ? counts[s] : '')).join('') +
    (charge ? (charge > 0 ? `${charge > 1 ? charge : ''}+` : `${charge < -1 ? -charge : ''}-`) : '');

  // Índice de deficiencia de hidrógeno (grados de insaturación)
  const C = (counts.C || 0) + (counts.Si || 0);
  // los metales alcalinos se cuentan como H (sustituyen a un H ácido)
  const H = (counts.H || 0) + (counts.Li || 0) + (counts.Na || 0) + (counts.K || 0);
  const N = (counts.N || 0) + (counts.P || 0);
  const X = (counts.F || 0) + (counts.Cl || 0) + (counts.Br || 0) + (counts.I || 0);
  const ihd = (2 * C + 2 + N - H - X) / 2;

  const desc = JSON.parse(mol.get_descriptors());
  let inchi = null;
  let inchikey = null;
  try {
    inchi = mol.get_inchi();
    inchikey = inchi ? RDKit.get_inchikey_for_inchi(inchi) : null;
  } catch { /* InChI opcional */ }

  // Tamaño del anillo más pequeño de cada átomo (0 = acíclico)
  const ringSize = new Array(atoms.length).fill(0);
  for (let r = 8; r >= 3; r--) {
    for (const m of matchRaw(`[r${r}]`)) ringSize[m.atoms[0]] = r;
  }

  function matchRaw(smarts) {
    const res = mol.get_substruct_matches(qmol(smarts));
    const arr = JSON.parse(res || '[]');
    return Array.isArray(arr) ? arr : [];
  }

  const ctx = {
    mol,
    atoms,
    bonds,
    formula,
    counts,
    ihd,
    charge,
    smiles: mol.get_smiles(),
    inchi,
    inchikey,
    mw: desc.amw,
    exactMass: desc.exactmw,
    ringSize,
    /** Coincidencias completas: [{atoms:[...], bonds:[...]}] */
    matchRaw,
    /** Solo listas de átomos */
    match(smarts) {
      return matchRaw(smarts).map((m) => m.atoms);
    },
    has(smarts) {
      return mol.get_substruct_match(qmol(smarts)) !== '{}';
    },
    /** Átomos distintos que ocupan la posición idx del patrón */
    uniqueAtoms(smarts, idx = 0) {
      return [...new Set(matchRaw(smarts).map((m) => m.atoms[idx]))];
    },
    atomSet(smarts, idx = 0) {
      return new Set(matchRaw(smarts).map((m) => m.atoms[idx]));
    },
    svg(highlightAtoms = [], width = 360, height = 260) {
      const details = {
        width,
        height,
        bondLineWidth: 1.6,
        addStereoAnnotation: false,
        clearBackground: false
      };
      if (highlightAtoms.length) {
        const set = new Set(highlightAtoms);
        details.atoms = [...set];
        details.bonds = bonds
          .map((b, i) => (set.has(b.atoms[0]) && set.has(b.atoms[1]) ? i : -1))
          .filter((i) => i >= 0);
        details.highlightColour = [0.98, 0.62, 0.2];
        return mol.get_svg_with_highlights(JSON.stringify(details));
      }
      return mol.get_svg_with_highlights(JSON.stringify(details));
    },
    dispose() {
      try { mol.delete(); } catch { /* ignorar */ }
    }
  };
  return ctx;
}

/** Convierte una estructura a SMILES canónico (o null si no es válida). */
export function toSmiles(input) {
  const mol = RDKit.get_mol(String(input || ''));
  if (!mol) return null;
  try {
    return mol.is_valid() ? mol.get_smiles() : null;
  } finally {
    mol.delete();
  }
}
