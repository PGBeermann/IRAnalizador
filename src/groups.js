// Base de conocimiento de frecuencias de grupo para espectroscopía IR.
// Cada regla detecta un grupo funcional con SMARTS (RDKit) y devuelve
// las bandas esperadas con su explicación química.
//
// eps   = absortividad relativa aproximada por cada aparición del grupo
//         (se usa para simular intensidades; el espectro se normaliza al final)
// label = intensidad descrita en tablas de correlación (mf, f, m, d, md, v)
// fwhm  = anchura a media altura en cm⁻¹

export const INT_LABEL = {
  mf: 'muy fuerte', f: 'fuerte', m: 'media', d: 'débil', md: 'muy débil', v: 'variable'
};

function b(center, eps, fwhm, label, vib, why, opts = {}) {
  const half = opts.halfRange ?? Math.max(20, fwhm / 2);
  return {
    center,
    eps,
    fwhm,
    label,
    vib,
    why,
    range: opts.range || [Math.round(center - half), Math.round(center + half)],
    gas: opts.gas || null,
    overtone: !!opts.overtone
  };
}

// ---------- utilidades ----------
function conjNeighbors(c, carbonIdx) {
  // vecinos del carbono carbonílico que lo conjugan (C=C, arilo, C≡C)
  const s = new Set();
  for (const m of c.matchRaw('[CX3](=[O,S])-[$([#6]=[#6]),$(c),$([#6]#[#6])]')) {
    if (m.atoms[0] === carbonIdx) s.add(m.atoms[2]);
  }
  return s.size;
}
function isAromaticAtom(c, idx) {
  return c.__arom ? c.__arom.has(idx) : (c.__arom = c.atomSet('a')).has(idx);
}
const RING_CO_KETONE = { 3: 1815, 4: 1780, 5: 1745, 6: 1715, 7: 1702 };
const RING_LACTONE = { 4: 1840, 5: 1770, 6: 1735 };
const RING_LACTAM = { 4: 1745, 5: 1700, 6: 1670 };
const RING_NAME = { 3: 'tres', 4: 'cuatro', 5: 'cinco', 6: 'seis', 7: 'siete' };

const WHY = {
  conj: 'La conjugación con un enlace C=C o un anillo aromático deslocaliza los electrones π, reduce el orden de enlace C=O y desplaza la banda ~25–30 cm⁻¹ hacia números de onda menores.',
  strain: (n) => `El carbonilo está en un anillo de ${RING_NAME[n] || n} miembros: la tensión angular aumenta el carácter s del enlace σ C=O y obliga a acoplarse con los enlaces del anillo, lo que eleva la frecuencia.`,
  hbond: 'Los puentes de hidrógeno debilitan el enlace O–H y generan una distribución de entornos distintos; por eso la banda es ancha e intensa.',
  gasOH: 'En fase gaseosa (o en disolución muy diluida) no hay puentes de hidrógeno: se observa el O–H "libre" como una banda aguda a mayor número de onda.'
};

// ===================================================================
//  REGLAS
// ===================================================================
export const RULES = [
  // ------------------------------------------------ CARBONILOS
  {
    id: 'acid_chloride', name: 'Halogenuro de acilo', category: 'Carbonilos',
    theory: 'El halógeno, muy electronegativo, retira densidad electrónica por efecto inductivo y apenas cede electrones por resonancia; el enlace C=O se acorta y se fortalece, por lo que los halogenuros de acilo presentan la frecuencia C=O más alta entre los derivados de ácido (~1800 cm⁻¹).',
    detect(c) {
      return c.matchRaw('[CX3](=O)[F,Cl,Br,I]').map((m) => {
        const conj = conjNeighbors(c, m.atoms[0]);
        const nu = conj ? 1775 : 1800;
        const bands = [b(nu, 1.1, 18, 'mf', 'Estiramiento C=O (ν C=O)',
          `Carbonilo de halogenuro de acilo${conj ? ' conjugado' : ''}. ${conj ? WHY.conj : 'El efecto inductivo –I del halógeno eleva la frecuencia por encima de la de cetonas y ésteres.'}`,
          { range: [1815 - (conj ? 40 : 0), 1815], gas: { center: nu + 10 } })];
        if (conj) bands.push(b(1735, 0.35, 18, 'm', 'Banda de resonancia de Fermi',
          'Los cloruros de aroilo muestran a menudo una segunda banda junto al C=O por resonancia de Fermi con el sobretono de la vibración C–C(=O) cercana a 870 cm⁻¹.'));
        bands.push(b(nu > 1790 ? 950 : 880, 0.3, 30, 'm', 'Estiramiento C–C(=O)–X',
          'Vibración del esqueleto acilo acoplada con el enlace C–halógeno.', { range: [850, 1000] }));
        return { atoms: m.atoms, bands, variant: conj ? 'conjugado' : 'alifático' };
      });
    }
  },
  {
    id: 'anhydride', name: 'Anhídrido de ácido', category: 'Carbonilos',
    theory: 'Los dos carbonilos unidos por un oxígeno vibran acoplados: el modo simétrico (en fase) y el asimétrico (fuera de fase) dan DOS bandas C=O separadas ~60 cm⁻¹. En anhídridos acíclicos la banda de mayor frecuencia es la más intensa; en los cíclicos ocurre lo contrario.',
    detect(c) {
      const seen = new Set();
      const out = [];
      for (const m of c.matchRaw('[CX3](=O)[OX2][CX3](=O)')) {
        const o = m.atoms[2];
        if (seen.has(o)) continue;
        seen.add(o);
        const ring = c.ringSize[o];
        const conj = Math.min(1, conjNeighbors(c, m.atoms[0]) + conjNeighbors(c, m.atoms[3]));
        const d = conj ? 20 : 0;
        const bands = ring === 5
          ? [
            b(1860 - d, 0.45, 18, 'm', 'Estiramiento C=O simétrico (en fase)', 'Anhídrido cíclico de cinco miembros: la tensión del anillo eleva ambas frecuencias; el modo simétrico es el más débil.'),
            b(1780 - d, 1.1, 20, 'mf', 'Estiramiento C=O asimétrico (fuera de fase)', 'En anhídridos cíclicos la banda de menor frecuencia (asimétrica) es la más intensa.'),
            b(1230, 0.7, 40, 'f', 'Estiramiento C–O–C', 'Estiramiento del puente C–O–C del anillo.', { range: [1300, 1175] }),
            b(920, 0.5, 30, 'm', 'Estiramiento C–O–C del anillo', 'Modo del anillo característico de anhídridos cíclicos.', { range: [950, 880] })
          ]
          : [
            b(1820 - d, 1.0, 18, 'f', 'Estiramiento C=O simétrico (en fase)', 'En anhídridos acíclicos la banda de mayor frecuencia es la más intensa.' + (conj ? ' ' + WHY.conj : '')),
            b(1760 - d, 0.7, 18, 'f', 'Estiramiento C=O asimétrico (fuera de fase)', 'Segunda banda carbonílica debida al acoplamiento mecánico de los dos C=O a través del oxígeno puente.'),
            b(1045, 0.8, 40, 'f', 'Estiramiento C–O–C', 'Estiramiento del enlace C–O–C del anhídrido.', { range: [1000, 1300] })
          ];
        out.push({ atoms: m.atoms, bands, variant: ring === 5 ? 'cíclico' : 'acíclico' });
      }
      return out;
    }
  },
  {
    id: 'carboxylic_acid', name: 'Ácido carboxílico', category: 'Carbonilos',
    theory: 'En fase condensada los ácidos carboxílicos forman dímeros cíclicos unidos por dos puentes de hidrógeno muy fuertes. Por eso la banda O–H es extraordinariamente ancha (3300–2500 cm⁻¹), se superpone a los C–H y muestra "hombros" por combinaciones; el C=O del dímero aparece cerca de 1710 cm⁻¹. Completan la firma el estiramiento C–O (~1290 cm⁻¹) y la flexión O–H fuera del plano, ancha, cerca de 930 cm⁻¹.',
    detect(c) {
      return c.matchRaw('[CX3](=O)[OX2H1]').map((m) => {
        const conj = conjNeighbors(c, m.atoms[0]);
        const nu = conj ? 1688 : 1710;
        return {
          atoms: m.atoms,
          variant: conj ? 'conjugado / aromático' : 'alifático',
          bands: [
            b(3000, 0.32, 480, 'f', 'Estiramiento O–H (dímero con puente de H)',
              'Banda muy ancha característica de los dímeros de ácido. ' + WHY.hbond + ' Sobre ella aparecen pequeños picos de combinación entre 2700 y 2500 cm⁻¹.',
              { range: [2500, 3300], gas: { center: 3580, fwhm: 25, eps: 0.25, vib: 'Estiramiento O–H libre (monómero)', why: 'En fase gaseosa el ácido está parcialmente como monómero: O–H libre, agudo, cerca de 3580 cm⁻¹.' } }),
            b(nu, 1.0, 30, 'mf', 'Estiramiento C=O (ν C=O)',
              `Carbonilo de ácido dimerizado. ${conj ? WHY.conj : 'El puente de hidrógeno del dímero rebaja la frecuencia respecto al monómero (~1760 cm⁻¹).'}`,
              { range: [1680 - (conj ? 10 : 0), 1725], gas: { center: conj ? 1752 : 1780, why: 'En fase gaseosa domina el monómero sin puentes de hidrógeno: el C=O aparece ~1780 cm⁻¹.' } }),
            b(1420, 0.3, 30, 'm', 'Flexión O–H en el plano (acoplada con ν C–O)',
              'Modo mixto de flexión O–H y estiramiento C–O.', { range: [1395, 1440] }),
            b(1290, 0.75, 40, 'f', 'Estiramiento C–O (acoplado con δ O–H)',
              'El enlace C–O del grupo carboxilo tiene carácter parcial de doble enlace y vibra entre 1320 y 1210 cm⁻¹.', { range: [1210, 1320], gas: { center: 1180 } }),
            b(935, 0.3, 70, 'm', 'Flexión O–H fuera del plano',
              'Banda ancha de intensidad media, diagnóstica del dímero de ácido carboxílico.', { range: [900, 960] })
          ]
        };
      });
    }
  },
  {
    id: 'carboxylate', name: 'Carboxilato (sal)', category: 'Carbonilos',
    theory: 'En el anión carboxilato la carga se deslocaliza por igual sobre los dos oxígenos: no hay C=O ni C–O aislados, sino dos enlaces equivalentes de orden 1,5 que vibran acoplados en un modo asimétrico (~1580 cm⁻¹) y otro simétrico (~1410 cm⁻¹).',
    detect(c) {
      return c.matchRaw('[CX3](=O)[O-]').map((m) => ({
        atoms: m.atoms,
        bands: [
          b(1580, 1.0, 40, 'mf', 'Estiramiento asimétrico COO⁻', 'Los dos enlaces C–O equivalentes (orden 1,5) vibran fuera de fase.', { range: [1550, 1650] }),
          b(1410, 0.6, 35, 'f', 'Estiramiento simétrico COO⁻', 'Vibración en fase de los dos enlaces C–O del carboxilato.', { range: [1300, 1440] })
        ]
      }));
    }
  },
  {
    id: 'ester', name: 'Éster / lactona', category: 'Carbonilos',
    theory: 'El oxígeno alcoxi retira densidad por efecto inductivo más de lo que la cede por resonancia, por lo que el C=O de los ésteres (~1740 cm⁻¹) está por encima del de las cetonas. La "regla de los tres picos" ayuda a reconocerlos: C=O ~1740, C–C(=O)–O ~1240–1160 (muy fuerte) y O–C–C ~1100–1030 cm⁻¹.',
    detect(c) {
      const smarts = '[$([CX3](=O)[#6]),$([CX3H1](=O))](=O)[OX2][#6;!$([CX3]=O)]';
      return c.matchRaw(smarts).map((m) => {
        const C = m.atoms[0];
        const altC = m.atoms[3];
        const ring = c.ringSize[C];
        const conj = conjNeighbors(c, C);
        const alkoxyUnsat = isAromaticAtom(c, altC) || c.atomSet('[OX2][#6]=[#6]', 1).has(altC);
        const isFormate = c.atoms[C].h === 1;
        const isAcetate = c.match('[CH3][CX3](=O)[OX2]').some((x) => x[1] === C);
        let nu = RING_LACTONE[ring] || 1740;
        const notes = [];
        if (RING_LACTONE[ring] && ring !== 6) notes.push(WHY.strain(ring));
        if (ring) notes.push(`Lactona (éster cíclico) de ${RING_NAME[ring] || ring} miembros.`);
        if (isFormate) { nu = 1725; notes.push('Los formiatos (H–C(=O)–O–R) absorben algo por debajo de los ésteres alquílicos.'); }
        if (conj) { nu -= 20; notes.push('Éster α,β-insaturado o aromático: ' + WHY.conj); }
        if (alkoxyUnsat) { nu += 25; notes.push('El oxígeno alcoxi está unido a un vinilo o arilo: su par libre se deslocaliza hacia ese grupo y deja de donar al C=O, que se refuerza (≈1765 cm⁻¹, p. ej. acetato de fenilo).'); }
        if (!notes.length) notes.push('Éster alifático saturado: valor típico 1750–1735 cm⁻¹.');
        const coC = isAcetate ? 1240 : conj ? 1275 : isFormate ? 1185 : 1190;
        const bands = [
          b(nu, 1.0, 22, 'mf', 'Estiramiento C=O (ν C=O)', notes.join(' '),
            { range: [nu - 15, nu + 15], gas: { center: nu + 15 } }),
          b(coC, 0.95, 40, 'f', 'Estiramiento asimétrico C–C(=O)–O',
            isAcetate ? 'Los acetatos muestran esta banda muy intensa cerca de 1240 cm⁻¹, a menudo más intensa que el propio C=O.'
              : conj ? 'En benzoatos y ésteres conjugados aparece cerca de 1280 cm⁻¹ (con otra banda cerca de 1110 cm⁻¹).'
                : 'Banda muy intensa típica de ésteres de ácidos alifáticos de cadena larga (~1190 cm⁻¹).',
            { range: [1160, 1300] }),
          b(conj ? 1110 : 1050, 0.55, 35, 'f', 'Estiramiento O–C–C (alcoxi)',
            'Estiramiento del enlace O–C del fragmento alcohólico; completa la "regla de los tres picos" del éster.',
            { range: [1030, 1120] })
        ];
        return {
          atoms: m.atoms,
          bands,
          variant: ring ? `lactona (anillo de ${ring})` : isFormate ? 'formiato' : isAcetate ? 'acetato' : conj ? 'conjugado' : alkoxyUnsat ? 'éster de vinilo/arilo' : 'alifático'
        };
      });
    }
  },
  {
    id: 'carbonate', name: 'Carbonato orgánico', category: 'Carbonilos',
    theory: 'Dos oxígenos alcoxi unidos al mismo carbonilo refuerzan el efecto inductivo: el C=O aparece ~1750 cm⁻¹ (≈1800 cm⁻¹ en carbonatos cíclicos de cinco miembros) y el C–O–C da una banda muy intensa cerca de 1260 cm⁻¹.',
    detect(c) {
      return c.matchRaw('[OX2]([#6])[CX3](=O)[OX2][#6]').map((m) => {
        const ring5 = c.ringSize[m.atoms[2]] === 5;
        return {
          atoms: m.atoms,
          bands: [
            b(ring5 ? 1800 : 1750, 1.0, 22, 'mf', 'Estiramiento C=O', ring5 ? WHY.strain(5) : 'Carbonato acíclico.'),
            b(1260, 1.0, 40, 'mf', 'Estiramiento asimétrico O–C(=O)–O', 'Banda muy intensa del esqueleto carbonato.', { range: [1210, 1300] })
          ]
        };
      });
    }
  },
  {
    id: 'carbamate', name: 'Carbamato (uretano)', category: 'Carbonilos',
    theory: 'El carbamato combina rasgos de amida y de éster: C=O entre 1740 y 1690 cm⁻¹, N–H (si existe) cerca de 3330 cm⁻¹, banda tipo "amida II" ~1530 cm⁻¹ y C–O intenso ~1250 cm⁻¹.',
    detect(c) {
      return c.matchRaw('[NX3][CX3](=O)[OX2][#6]').map((m) => {
        const nh = c.atoms[m.atoms[0]].h > 0;
        const bands = [
          b(1705, 1.0, 25, 'mf', 'Estiramiento C=O', 'Carbonilo de carbamato: intermedio entre amida y éster.', { range: [1690, 1740] }),
          b(1250, 0.85, 40, 'f', 'Estiramiento C–O / C–N acoplado', 'Vibración del esqueleto O–C(=O)–N.', { range: [1200, 1300] })
        ];
        if (nh) {
          bands.push(b(3330, 0.3, 60, 'm', 'Estiramiento N–H', 'N–H de carbamato secundario, ensanchado por puentes de hidrógeno.', { range: [3250, 3450], gas: { center: 3450, fwhm: 25 } }));
          bands.push(b(1530, 0.6, 30, 'f', 'Flexión N–H + estiramiento C–N (tipo amida II)', 'Modo acoplado típico de amidas y carbamatos secundarios.', { range: [1510, 1550] }));
        }
        return { atoms: m.atoms, bands };
      });
    }
  },
  {
    id: 'urea', name: 'Urea', category: 'Carbonilos',
    theory: 'Dos nitrógenos donan densidad electrónica por resonancia al carbonilo, que adquiere mucho carácter de enlace simple: el C=O de las ureas aparece bajo, ~1660 cm⁻¹, acompañado de N–H y flexiones N–H ~1600 cm⁻¹.',
    detect(c) {
      return c.matchRaw('[NX3][CX3](=O)[NX3]').map((m) => {
        const nh = c.atoms[m.atoms[0]].h + c.atoms[m.atoms[3]].h;
        const bands = [b(1660, 1.0, 30, 'mf', 'Estiramiento C=O (amida I)', 'La doble donación por resonancia de los dos nitrógenos reduce el orden de enlace C=O.', { range: [1630, 1690] })];
        if (nh) {
          bands.push(b(3350, 0.35, 90, 'f', 'Estiramiento N–H', 'N–H asociado por puentes de hidrógeno.', { range: [3200, 3450] }));
          bands.push(b(1600, 0.5, 35, 'f', 'Flexión N–H (amida II)', 'Flexión N–H acoplada con el estiramiento C–N.', { range: [1550, 1640] }));
        }
        return { atoms: m.atoms, bands };
      });
    }
  },
  {
    id: 'imide', name: 'Imida', category: 'Carbonilos',
    theory: 'Como en los anhídridos, los dos carbonilos unidos al mismo nitrógeno se acoplan y producen dos bandas: una débil ~1770 cm⁻¹ (simétrica) y otra intensa ~1700 cm⁻¹ (asimétrica).',
    detect(c) {
      const seen = new Set();
      return c.matchRaw('[CX3](=O)[NX3][CX3](=O)').filter((m) => !seen.has(m.atoms[2]) && seen.add(m.atoms[2])).map((m) => {
        const nh = c.atoms[m.atoms[2]].h > 0;
        const bands = [
          b(1770, 0.4, 20, 'm', 'Estiramiento C=O simétrico', 'Modo en fase de los dos carbonilos de la imida.', { range: [1740, 1790] }),
          b(1705, 1.0, 25, 'mf', 'Estiramiento C=O asimétrico', 'Modo fuera de fase, el más intenso.', { range: [1680, 1725] })
        ];
        if (nh) bands.push(b(3200, 0.3, 80, 'm', 'Estiramiento N–H', 'N–H de imida, ácido y asociado por puentes de hidrógeno.', { range: [3100, 3300] }));
        return { atoms: m.atoms, bands };
      });
    }
  },
  {
    id: 'amide', name: 'Amida / lactama', category: 'Carbonilos',
    theory: 'El par libre del nitrógeno se deslocaliza hacia el carbonilo (resonancia N–C=O ↔ N⁺=C–O⁻), lo que debilita el C=O: la banda "amida I" aparece entre 1690 y 1630 cm⁻¹. Las amidas primarias muestran dos N–H (asimétrico y simétrico) y la flexión NH₂ (~1620 cm⁻¹); las secundarias, un único N–H (~3300 cm⁻¹) y la intensa banda "amida II" (~1550 cm⁻¹, flexión N–H + estiramiento C–N); las terciarias no presentan N–H.',
    detect(c) {
      const smarts = '[$([CX3](=O)[#6]),$([CX3H1](=O))](=O)[NX3]';
      const imideN = c.atomSet('[CX3](=O)[NX3][CX3](=O)', 2);
      return c.matchRaw(smarts).filter((m) => !imideN.has(m.atoms[2])).map((m) => {
        const C = m.atoms[0];
        const N = m.atoms[2];
        const nH = c.atoms[N].h;
        const ring = c.ringSize[C] && c.ringSize[N] ? c.ringSize[C] : 0;
        const conj = conjNeighbors(c, C);
        let nu = nH === 2 ? 1665 : nH === 1 ? 1650 : 1650;
        const notes = [];
        if (ring && RING_LACTAM[ring]) {
          nu = RING_LACTAM[ring];
          notes.push(`Lactama de ${RING_NAME[ring]} miembros.`);
          if (ring < 6) notes.push(WHY.strain(ring));
        }
        if (conj) { nu -= 12; notes.push('La conjugación con el grupo acilo (arilo o C=C) baja ligeramente la frecuencia.'); }
        notes.unshift(nH === 2 ? 'Amida primaria.' : nH === 1 ? 'Amida secundaria.' : 'Amida terciaria (sin N–H).');
        notes.push('La deslocalización del par libre del N reduce el orden de enlace C=O (banda amida I).');
        const bands = [b(nu, 1.0, 30, 'mf', 'Estiramiento C=O (banda amida I)', notes.join(' '),
          { range: [nu - 25, nu + 25], gas: { center: nu + 45, why: 'Sin puentes de hidrógeno (gas o disolución diluida) la banda amida I sube a ~1700 cm⁻¹.' } })];
        if (nH === 2) {
          bands.push(b(3350, 0.4, 60, 'm', 'Estiramiento N–H asimétrico (NH₂)', 'Los dos N–H del grupo NH₂ vibran fuera de fase.', { range: [3300, 3400], gas: { center: 3550, fwhm: 25 } }));
          bands.push(b(3180, 0.4, 70, 'm', 'Estiramiento N–H simétrico (NH₂)', 'Modo en fase; en estado sólido ambos N–H se ensanchan por puentes de hidrógeno.', { range: [3150, 3220], gas: { center: 3430, fwhm: 25 } }));
          bands.push(b(1625, 0.6, 30, 'f', 'Flexión NH₂ (banda amida II)', 'Flexión de tijera del NH₂; aparece como hombro o banda junto a la amida I.', { range: [1590, 1650] }));
          bands.push(b(1405, 0.4, 30, 'm', 'Estiramiento C–N', 'Estiramiento C–N con carácter parcial de doble enlace.', { range: [1380, 1420] }));
        } else if (nH === 1) {
          bands.push(b(3290, 0.4, 70, 'm', 'Estiramiento N–H', 'Un solo N–H, asociado por puentes de hidrógeno N–H···O=C.', { range: [3250, 3330], gas: { center: 3450, fwhm: 25 } }));
          bands.push(b(3080, 0.1, 40, 'd', 'Sobretono de amida II (2 × ~1550)', 'Banda débil por resonancia de Fermi con el sobretono de la amida II.', { range: [3060, 3100], overtone: true }));
          if (!(ring && ring < 7)) bands.push(b(1550, 0.75, 30, 'f', 'Banda amida II (δ N–H + ν C–N)', 'Exclusiva de amidas secundarias en conformación trans; en lactamas pequeñas (cis) no aparece.', { range: [1515, 1570] }));
          bands.push(b(1260, 0.25, 40, 'm', 'Banda amida III', 'Modo mixto ν C–N + δ N–H.', { range: [1230, 1300] }));
          bands.push(b(700, 0.25, 80, 'm', 'Flexión N–H fuera del plano', 'Banda ancha de "aleteo" N–H.', { range: [650, 750] }));
        } else {
          bands.push(b(1400, 0.35, 35, 'm', 'Estiramiento C–N', 'Estiramiento C–N de amida terciaria.', { range: [1370, 1420] }));
        }
        return { atoms: m.atoms, bands, variant: (ring ? 'lactama ' : '') + (nH === 2 ? 'primaria' : nH === 1 ? 'secundaria' : 'terciaria') };
      });
    }
  },
  {
    id: 'aldehyde', name: 'Aldehído', category: 'Carbonilos',
    theory: 'Además del C=O (~1730 cm⁻¹ alifático, ~1700 cm⁻¹ aromático), los aldehídos muestran el estiramiento C–H del grupo CHO como un DOBLETE (~2820 y ~2720 cm⁻¹) producido por resonancia de Fermi con el sobretono de la flexión C–H (~1390 cm⁻¹). La banda de 2720 cm⁻¹, aislada a la derecha de los C–H alquílicos, es muy diagnóstica.',
    detect(c) {
      return c.matchRaw('[CX3H1](=O)[#6]').map((m) => {
        const conj = conjNeighbors(c, m.atoms[0]);
        const nu = conj ? 1702 : 1730;
        return {
          atoms: m.atoms,
          variant: conj ? 'conjugado / aromático' : 'alifático',
          bands: [
            b(nu, 1.0, 22, 'mf', 'Estiramiento C=O (ν C=O)', conj ? 'Aldehído conjugado. ' + WHY.conj : 'Aldehído alifático saturado: algo por encima de las cetonas porque el H no aporta el efecto donador (+I) de un grupo alquilo.',
              { range: [nu - 20, nu + 15], gas: { center: nu + 15 } }),
            b(2820, 0.22, 25, 'm', 'Estiramiento C–H de aldehído (doblete de Fermi)', 'Componente de alta frecuencia del doblete de Fermi; puede quedar parcialmente oculta bajo los C–H alquílicos.', { range: [2800, 2860] }),
            b(2720, 0.2, 25, 'm', 'Estiramiento C–H de aldehído (doblete de Fermi)', 'Componente de baja frecuencia del doblete: la resonancia de Fermi con el sobretono 2 × δ C–H (~1390 cm⁻¹) desdobla el ν C–H en dos bandas. Muy diagnóstica de aldehídos.', { range: [2695, 2775] }),
            b(1390, 0.2, 25, 'd', 'Flexión C–H del grupo CHO', 'Flexión en el plano del C–H aldehídico; su sobretono participa en la resonancia de Fermi.', { range: [1380, 1410] }),
            b(2 * nu - 20, 0.04, 30, 'md', 'Sobretono 2 × ν C=O', 'Primer sobretono del C=O, débil, cerca del doble de su frecuencia.', { overtone: true })
          ]
        };
      });
    }
  },
  {
    id: 'ketone', name: 'Cetona', category: 'Carbonilos',
    theory: 'La cetona alifática saturada es la referencia de las frecuencias carbonílicas: ν C=O ≈ 1715 cm⁻¹, banda muy intensa porque el enlace C=O es muy polar y su estiramiento produce un gran cambio de momento dipolar. La conjugación la baja (~1690 cm⁻¹, arilcetonas) y la tensión de anillo la sube (ciclopentanona ~1745, ciclobutanona ~1780 cm⁻¹).',
    detect(c) {
      return c.matchRaw('[#6][CX3](=O)[#6]').filter((m, i, arr) => arr.findIndex((x) => x.atoms[1] === m.atoms[1]) === i).map((m) => {
        const C = m.atoms[1];
        const ring = c.ringSize[C];
        const conj = conjNeighbors(c, C);
        let nu = RING_CO_KETONE[ring] || 1715;
        const notes = [];
        if (ring && ring !== 6) notes.push(WHY.strain(ring));
        if (ring === 6) notes.push('Ciclohexanona: el anillo de seis miembros no está tensionado y absorbe como una cetona acíclica.');
        if (conj) {
          nu -= 25 * Math.min(conj, 2);
          notes.push((conj > 1 ? 'Cetona doblemente conjugada (p. ej., diarilcetona o quinona). ' : 'Cetona α,β-insaturada o arilcetona. ') + WHY.conj);
        }
        const halo = c.match('[CX3](=O)[CX4][F,Cl,Br]').some((x) => x[0] === C);
        if (halo) { nu += 18; notes.push('Un halógeno en α (efecto de campo/inductivo) eleva la frecuencia ~15–25 cm⁻¹.'); }
        if (!notes.length) notes.push('Cetona alifática saturada acíclica: valor de referencia 1715 cm⁻¹.');
        return {
          atoms: [m.atoms[0], C, m.atoms[2], m.atoms[3]],
          variant: (ring ? `cíclica (${ring})` : 'acíclica') + (conj ? ', conjugada' : ''),
          bands: [
            b(nu, 1.0, 22, 'mf', 'Estiramiento C=O (ν C=O)', notes.join(' '), { range: [nu - 15, nu + 15], gas: { center: nu + 15 } }),
            b(conj ? 1265 : 1220, 0.4, 40, 'm', 'Estiramiento asimétrico C–C(=O)–C', 'Vibración del esqueleto alrededor del carbonilo.', { range: [1100, 1300] }),
            b(2 * nu - 20, 0.04, 30, 'md', 'Sobretono 2 × ν C=O', 'Primer sobretono del carbonilo: banda débil cerca del doble de la frecuencia fundamental (3500–3350 cm⁻¹); no debe confundirse con un O–H.', { overtone: true })
          ]
        };
      });
    }
  },
  {
    id: 'ring_carbonyl', name: 'Carbonilo en heterociclo (tipo lactama/urea)', category: 'Carbonilos',
    theory: 'En heterociclos como piridonas, uracilos o xantinas (cafeína, teofilina) el C=O forma parte de un anillo con nitrógenos que donan densidad por resonancia: se comporta como una amida o urea cíclica, con bandas C=O intensas entre 1720 y 1640 cm⁻¹. Cuando hay dos carbonilos vecinos suelen verse dos bandas por acoplamiento.',
    detect(c) {
      const nAr = c.atomSet('n');
      return c.matchRaw('[#6;a]=O').map((m) => {
        const C = m.atoms[0];
        const nbrs = c.bonds.filter((bd) => bd.atoms.includes(C)).map((bd) => bd.atoms.find((x) => x !== C));
        const nN = nbrs.filter((x) => nAr.has(x)).length;
        const nu = nN >= 2 ? 1700 : nN === 1 ? 1660 : 1640;
        const kind = nN >= 2 ? 'tipo urea cíclica (N–C(=O)–N)' : nN === 1 ? 'tipo lactama (N–C(=O)–C)' : 'tipo pirona/tropona';
        const bands = [b(nu, 1.0, 25, 'mf', 'Estiramiento C=O en anillo heteroaromático', `Carbonilo ${kind}. La donación electrónica de los nitrógenos vecinos y la conjugación con el anillo bajan la frecuencia respecto a una cetona.`, { range: [nu - 25, nu + 25] })];
        if (nbrs.some((x) => nAr.has(x) && c.atoms[x].h > 0)) {
          bands.push(b(3160, 0.3, 150, 'm', 'Estiramiento N–H del anillo (asociado)', 'N–H de tipo lactama, unido por puentes de hidrógeno N–H···O=C: banda ancha 3250–3000 cm⁻¹.', { range: [3000, 3250] }));
        }
        return { atoms: m.atoms, variant: kind, bands };
      });
    }
  },
  {
    id: 'thiocarbonyl', name: 'Tiocarbonilo (C=S)', category: 'Azufre',
    theory: 'El azufre es mucho más pesado que el oxígeno y el enlace C=S es más débil, por lo que su estiramiento aparece a frecuencias mucho menores (1250–1020 cm⁻¹).',
    detect(c) {
      return c.matchRaw('[#6,#7][CX3]=[SX1]').map((m) => ({
        atoms: m.atoms, bands: [b(1120, 0.6, 40, 'f', 'Estiramiento C=S', 'Mayor masa reducida y menor constante de fuerza que el C=O.', { range: [1020, 1250] })]
      }));
    }
  },

  // ------------------------------------------------ O–H y C–O
  {
    id: 'alcohol', name: 'Alcohol', category: 'Oxigenados',
    theory: 'El O–H asociado por puentes de hidrógeno da una banda intensa y ancha centrada ~3350 cm⁻¹; en gas o disolución diluida aparece el O–H libre, agudo, ~3650 cm⁻¹. El estiramiento C–O (1260–1000 cm⁻¹) permite distinguir el tipo de alcohol: primario ~1050, secundario ~1100, terciario ~1150 cm⁻¹.',
    detect(c) {
      return c.matchRaw('[OX2H][CX4]').map((m) => {
        const C = m.atoms[1];
        const h = c.atoms[C].h;
        const type = h === 3 ? 'metanol' : h === 2 ? 'primario' : h === 1 ? 'secundario' : 'terciario';
        const co = { metanol: 1030, primario: 1050, secundario: 1100, terciario: 1150 }[type];
        return {
          atoms: m.atoms,
          variant: type,
          bands: [
            b(3340, 0.5, 260, 'f', 'Estiramiento O–H (asociado por puentes de H)', WHY.hbond,
              { range: [3200, 3550], gas: { center: 3655, fwhm: 22, eps: 0.2, why: WHY.gasOH, vib: 'Estiramiento O–H libre (sin puentes de H)' } }),
            b(co, 0.6, 40, 'f', `Estiramiento C–O (alcohol ${type})`,
              `En alcoholes ${type === 'metanol' ? 'como el metanol' : type + 's'} el estiramiento C–O aparece cerca de ${co} cm⁻¹; sube con la ramificación del carbono porque el C–O se acopla con más enlaces C–C.`,
              { range: [co - 30, co + 30] }),
            b(1350, 0.25, 70, 'm', 'Flexión O–H en el plano', 'Banda ancha de intensidad media, acoplada con flexiones C–H.', { range: [1260, 1420] }),
            b(660, 0.25, 110, 'm', 'Flexión O–H fuera del plano', 'Banda ancha en la zona de baja frecuencia, propia de alcoholes asociados.', { range: [600, 720], gas: { eps: 0.03, fwhm: 40, why: 'En gas, sin asociación, la torsión O–H aparece por debajo de 400 cm⁻¹ y esta banda prácticamente desaparece.' } })
          ]
        };
      });
    }
  },
  {
    id: 'phenol', name: 'Fenol', category: 'Oxigenados',
    theory: 'Los fenoles muestran el O–H asociado (~3350 cm⁻¹, ancho) y un estiramiento C–O más alto que en alcoholes (~1230 cm⁻¹) porque el oxígeno conjuga con el anillo y el enlace C(arilo)–O adquiere carácter parcial de doble enlace.',
    detect(c) {
      return c.matchRaw('[OX2H]c').map((m) => ({
        atoms: m.atoms,
        bands: [
          b(3350, 0.45, 250, 'f', 'Estiramiento O–H (asociado por puentes de H)', WHY.hbond, { range: [3200, 3550], gas: { center: 3655, fwhm: 22, eps: 0.2, why: WHY.gasOH, vib: 'Estiramiento O–H libre (sin puentes de H)' } }),
          b(1230, 0.8, 40, 'f', 'Estiramiento C(arilo)–O', 'La conjugación del par libre del oxígeno con el anillo refuerza el enlace C–O.', { range: [1180, 1260] }),
          b(1360, 0.3, 50, 'm', 'Flexión O–H en el plano', 'Flexión O–H acoplada con el estiramiento C–O.', { range: [1310, 1410] })
        ]
      }));
    }
  },
  {
    id: 'ether', name: 'Éter', category: 'Oxigenados',
    theory: 'Los éteres no tienen O–H ni C=O; su firma es el estiramiento asimétrico C–O–C, intenso, entre 1300 y 1000 cm⁻¹. En éteres dialquílicos aparece ~1120 cm⁻¹; en éteres aril-alquílicos se desdobla en ~1250 (asimétrico, C(arilo)–O con carácter de doble enlace) y ~1040 cm⁻¹ (simétrico).',
    detect(c) {
      const arom = c.atomSet('a');
      const vinyl = c.atomSet('[CX3]=[CX3]');
      const seen = new Set();
      const out = [];
      for (const m of c.matchRaw('[OX2;!r3]([#6;!$([#6]=[O,S,N])])[#6;!$([#6]=[O,S,N])]')) {
        const [o, a1, a2] = m.atoms;
        if (seen.has(o)) continue;
        seen.add(o);
        const ar = [a1, a2].filter((x) => arom.has(x)).length;
        const vi = [a1, a2].filter((x) => vinyl.has(x)).length;
        let bands;
        let variant;
        if (ar === 2) {
          variant = 'diarílico';
          bands = [b(1240, 0.9, 40, 'f', 'Estiramiento asimétrico C(arilo)–O–C(arilo)', 'Ambos enlaces C–O conjugados con anillos aromáticos.', { range: [1200, 1270] })];
        } else if (ar === 1) {
          variant = 'aril-alquílico';
          bands = [
            b(1250, 0.9, 35, 'f', 'Estiramiento asimétrico =C–O–C', 'El enlace C(arilo)–O tiene carácter parcial de doble enlace por resonancia y vibra a mayor frecuencia (p. ej., anisol 1245 cm⁻¹).', { range: [1230, 1275] }),
            b(1040, 0.55, 30, 'm', 'Estiramiento simétrico =C–O–C', 'Componente simétrica del grupo aril-alquil éter.', { range: [1010, 1075] })
          ];
        } else if (vi) {
          variant = 'vinílico';
          bands = [
            b(1220, 0.85, 35, 'f', 'Estiramiento asimétrico =C–O–C', 'Éter vinílico: conjugación del oxígeno con el C=C.', { range: [1200, 1240] }),
            b(1060, 0.4, 35, 'm', 'Estiramiento simétrico =C–O–C', 'Componente simétrica.', { range: [1020, 1080] })
          ];
        } else {
          const ring = c.ringSize[o];
          variant = ring ? `cíclico (${ring})` : 'dialquílico';
          const nu = ring === 5 ? 1070 : ring === 4 ? 980 : 1120;
          bands = [b(nu, 0.85, 40, 'f', 'Estiramiento asimétrico C–O–C', ring ? `Éter cíclico de ${RING_NAME[ring] || ring} miembros.` : 'Éter dialquílico: el modo asimétrico C–O–C es la banda más intensa del espectro en la región 1150–1085 cm⁻¹.', { range: [nu - 40, nu + 30] })];
        }
        out.push({ atoms: m.atoms, bands, variant });
      }
      return out;
    }
  },
  {
    id: 'epoxide', name: 'Epóxido (oxirano)', category: 'Oxigenados',
    theory: 'El anillo de tres miembros muy tensionado da bandas de "respiración" del anillo (~1250 cm⁻¹) y deformaciones asimétrica (~915 cm⁻¹) y simétrica (~840 cm⁻¹); los C–H del anillo absorben por encima de 3000 cm⁻¹.',
    detect(c) {
      return c.matchRaw('[OX2r3]1[#6r3][#6r3]1').map((m) => ({
        atoms: m.atoms,
        bands: [
          b(3045, 0.1, 25, 'd', 'Estiramiento C–H del anillo epóxido', 'Los C–H de anillos de tres miembros tienen más carácter s y vibran por encima de 3000 cm⁻¹.', { range: [3000, 3070] }),
          b(1255, 0.45, 30, 'm', 'Respiración del anillo (simétrico)', 'Todos los enlaces del anillo se estiran en fase.', { range: [1230, 1280] }),
          b(915, 0.6, 30, 'f', 'Deformación asimétrica del anillo', 'Banda intensa típica de epóxidos.', { range: [880, 950] }),
          b(840, 0.5, 30, 'f', 'Deformación simétrica del anillo', 'Segunda banda de deformación del oxirano.', { range: [810, 870] })
        ]
      }));
    }
  },

  // ------------------------------------------------ NITROGENADOS
  {
    id: 'amine', name: 'Amina', category: 'Nitrogenados',
    theory: 'El enlace N–H es menos polar que el O–H, así que sus bandas son más débiles y estrechas. Aminas primarias: dos bandas N–H (asimétrica y simétrica, ~3400 y ~3300 cm⁻¹) y flexión NH₂ ~1600 cm⁻¹. Secundarias: una sola banda N–H débil. Terciarias: ninguna. El estiramiento C–N aparece a 1250–1020 cm⁻¹ en aminas alifáticas y a 1340–1250 cm⁻¹ en aromáticas.',
    detect(c) {
      const arom = c.atomSet('a');
      const smarts = '[NX3;+0;!$(N[#6,#16,#15]=[O,S,N]);!$(N[#7,#8]);!$(N[SX4])]';
      return c.matchRaw(smarts).map((m) => {
        const N = m.atoms[0];
        const h = c.atoms[N].h;
        const nbrs = c.bonds.filter((bd) => bd.atoms.includes(N)).map((bd) => bd.atoms.find((x) => x !== N));
        if (!nbrs.length && h === 3) return null; // amoníaco
        const aryl = nbrs.some((x) => arom.has(x));
        const bands = [];
        let variant;
        if (h === 2) {
          variant = aryl ? 'primaria aromática' : 'primaria alifática';
          bands.push(b(aryl ? 3440 : 3370, 0.22, 45, 'm', 'Estiramiento N–H asimétrico (NH₂)', 'Primera de las dos bandas N–H típicas de las aminas primarias.' + (aryl ? ' En anilinas aparecen más altas porque el par libre conjuga con el anillo.' : ''), { range: aryl ? [3420, 3500] : [3330, 3400], gas: { center: (aryl ? 3500 : 3410), fwhm: 20 } }));
          bands.push(b(aryl ? 3360 : 3290, 0.2, 45, 'm', 'Estiramiento N–H simétrico (NH₂)', 'Segunda banda N–H; la pareja de bandas es diagnóstica de NH₂.', { range: aryl ? [3330, 3420] : [3250, 3330], gas: { center: (aryl ? 3415 : 3340), fwhm: 20 } }));
          bands.push(b(aryl ? 1620 : 1600, 0.45, 35, 'm', 'Flexión de tijera NH₂', 'Deformación angular H–N–H; en anilinas puede solaparse con las bandas del anillo.', { range: [1580, 1650] }));
          bands.push(b(aryl ? 1280 : 1080, 0.45, 40, aryl ? 'f' : 'm', 'Estiramiento C–N', aryl ? 'C(arilo)–N con carácter parcial de doble enlace por resonancia.' : 'C–N alifático.', { range: aryl ? [1250, 1340] : [1020, 1250] }));
          if (!aryl) bands.push(b(800, 0.35, 90, 'm', 'Aleteo NH₂ fuera del plano', 'Banda ancha característica de aminas primarias y secundarias líquidas.', { range: [665, 910] }));
        } else if (h === 1) {
          variant = aryl ? 'secundaria aromática' : 'secundaria alifática';
          bands.push(b(aryl ? 3410 : 3310, aryl ? 0.25 : 0.12, 40, aryl ? 'm' : 'd', 'Estiramiento N–H', 'Amina secundaria: una sola banda N–H, débil en aminas alifáticas.', { range: aryl ? [3380, 3450] : [3290, 3350] }));
          bands.push(b(aryl ? 1310 : 1130, 0.45, 40, aryl ? 'f' : 'm', 'Estiramiento C–N', aryl ? 'C(arilo)–N.' : 'C–N alifático.', { range: aryl ? [1250, 1340] : [1020, 1250] }));
          bands.push(b(720, 0.25, 80, 'm', 'Aleteo N–H fuera del plano', 'Banda ancha.', { range: [665, 750] }));
        } else {
          variant = aryl ? 'terciaria aromática' : 'terciaria alifática';
          bands.push(b(aryl ? 1350 : 1150, 0.45, 40, 'm', 'Estiramiento C–N', 'Amina terciaria: no presenta bandas N–H; sólo el estiramiento C–N.', { range: aryl ? [1310, 1380] : [1020, 1250] }));
          if (!aryl) bands.push(b(2780, 0.12, 30, 'd', 'Bandas de Bohlmann (C–H en α al N)', 'Los C–H antiperiplanares al par libre del N se debilitan por hiperconjugación y absorben por debajo de 2800 cm⁻¹.', { range: [2700, 2800] }));
        }
        return { atoms: [N, ...nbrs], bands, variant };
      }).filter(Boolean);
    }
  },
  {
    id: 'ammonium', name: 'Sal de amonio', category: 'Nitrogenados',
    theory: 'En las sales de amonio el N–H⁺ forma puentes de hidrógeno muy fuertes con el contraión: aparece una banda muy ancha y compleja entre 3200 y 2250 cm⁻¹, y flexiones N–H⁺ hacia 1600–1500 cm⁻¹.',
    detect(c) {
      return c.matchRaw('[NX4+;!H0]').map((m) => ({
        atoms: m.atoms,
        bands: [
          b(2950, 0.4, 350, 'f', 'Estiramiento N–H⁺ (asociado)', 'Banda muy ancha con estructura fina por combinaciones.', { range: [2250, 3200] }),
          b(1580, 0.35, 40, 'm', 'Flexión N–H⁺ asimétrica', 'Deformación del grupo NH₃⁺ / NH₂⁺.', { range: [1560, 1625] }),
          b(1500, 0.3, 40, 'm', 'Flexión N–H⁺ simétrica', 'Deformación simétrica.', { range: [1480, 1550] })
        ]
      }));
    }
  },
  {
    id: 'nitrile', name: 'Nitrilo', category: 'Nitrogenados',
    theory: 'El triple enlace C≡N es muy fuerte y los átomos son ligeros: absorbe en la región de los triples enlaces (2260–2210 cm⁻¹), una zona casi vacía del espectro, por lo que es muy diagnóstico. La conjugación lo baja ~25 cm⁻¹ y aumenta su intensidad.',
    detect(c) {
      const conjSet = c.atomSet('[CX2](#[NX1])[$(c),$([#6]=[#6])]');
      return c.matchRaw('[CX2]#[NX1]').map((m) => {
        const conj = conjSet.has(m.atoms[0]);
        return {
          atoms: m.atoms,
          variant: conj ? 'conjugado / aromático' : 'alifático',
          bands: [b(conj ? 2228 : 2250, conj ? 0.5 : 0.3, 14, conj ? 'f' : 'm', 'Estiramiento C≡N',
            conj ? 'Nitrilo conjugado: la deslocalización π baja la frecuencia y aumenta el cambio de dipolo (banda más intensa).' : 'Nitrilo alifático saturado: banda aguda de intensidad media cerca de 2250 cm⁻¹.',
            { range: conj ? [2215, 2240] : [2240, 2265], gas: { center: (conj ? 2235 : 2262) } })]
        };
      });
    }
  },
  {
    id: 'cumulated', name: 'Doble enlace acumulado (X=C=Y)', category: 'Nitrogenados',
    theory: 'Los sistemas acumulados (isocianatos, isotiocianatos, azidas, carbodiimidas, alenos) tienen un estiramiento asimétrico muy intenso en la región de 2300–1900 cm⁻¹, donde casi nada más absorbe.',
    detect(c) {
      const defs = [
        ['[NX2]=[CX2]=[OX1]', 2270, 'Estiramiento asimétrico N=C=O (isocianato)', 1.2, 40, [2250, 2285]],
        ['[NX2]=[CX2]=[SX1]', 2100, 'Estiramiento asimétrico N=C=S (isotiocianato)', 1.1, 90, [2040, 2150]],
        ['[NX2]=[NX2+]=[NX1-]', 2100, 'Estiramiento asimétrico N=N=N (azida)', 1.1, 30, [2080, 2170]],
        ['[NX2]=[CX2]=[NX2]', 2130, 'Estiramiento asimétrico N=C=N (carbodiimida)', 1.1, 30, [2100, 2150]],
        ['[CX3]=[CX2]=[CX3]', 1950, 'Estiramiento asimétrico C=C=C (aleno)', 0.4, 20, [1930, 1970]],
        ['[CX1-]#[NX2+]', 2140, 'Estiramiento C≡N⁺ (isonitrilo)', 0.9, 20, [2110, 2165]]
      ];
      const out = [];
      for (const [s, nu, vib, eps, w, range] of defs) {
        for (const m of c.matchRaw(s)) {
          out.push({ atoms: m.atoms, variant: vib.split('(')[1]?.replace(')', ''), bands: [b(nu, eps, w, 'mf', vib, 'Vibración asimétrica del sistema acumulado: los dos enlaces dobles se estiran fuera de fase y producen un gran cambio de momento dipolar.', { range })] });
        }
      }
      return out;
    }
  },
  {
    id: 'nitro', name: 'Nitro', category: 'Nitrogenados',
    theory: 'El grupo NO₂ tiene dos enlaces N–O equivalentes (orden 1,5) que vibran acoplados: estiramiento asimétrico muy intenso (~1550 alifático / ~1525 aromático) y simétrico (~1375 / ~1345 cm⁻¹). Ambas bandas son de las más intensas del espectro.',
    detect(c) {
      const arom = c.atomSet('a');
      return c.matchRaw('[#6][$([NX3+](=O)[O-]),$([NX3](=O)=O)]').map((m) => {
        const ar = arom.has(m.atoms[0]);
        return {
          atoms: m.atoms,
          variant: ar ? 'aromático' : 'alifático',
          bands: [
            b(ar ? 1525 : 1550, 1.0, 25, 'mf', 'Estiramiento asimétrico NO₂', ar ? 'La conjugación con el anillo baja ~25 cm⁻¹ respecto a los nitroalcanos.' : 'Nitroalcano.', { range: ar ? [1500, 1550] : [1530, 1570] }),
            b(ar ? 1348 : 1375, 0.9, 25, 'f', 'Estiramiento simétrico NO₂', 'Componente simétrica del par de bandas del nitro.', { range: ar ? [1330, 1370] : [1360, 1390] }),
            b(ar ? 850 : 870, 0.3, 25, 'm', 'Estiramiento C–N', 'Estiramiento C–NO₂ (región de huella dactilar).', { range: [830, 890] })
          ]
        };
      });
    }
  },
  {
    id: 'nitrate', name: 'Éster nitrato', category: 'Nitrogenados',
    theory: 'Los nitratos orgánicos (R–O–NO₂) muestran estiramientos NO₂ muy intensos a ~1630 y ~1280 cm⁻¹ y el O–N a ~860 cm⁻¹.',
    detect(c) {
      return c.matchRaw('[OX2][$([NX3+](=O)[O-]),$([NX3](=O)=O)]').map((m) => ({
        atoms: m.atoms,
        bands: [
          b(1630, 1.0, 25, 'mf', 'Estiramiento asimétrico NO₂', 'Nitrato orgánico.', { range: [1610, 1660] }),
          b(1280, 0.9, 25, 'f', 'Estiramiento simétrico NO₂', 'Componente simétrica.', { range: [1260, 1300] }),
          b(860, 0.6, 25, 'f', 'Estiramiento O–N', 'Enlace O–NO₂.', { range: [830, 870] })
        ]
      }));
    }
  },
  {
    id: 'imine', name: 'Imina / oxima', category: 'Nitrogenados',
    theory: 'El doble enlace C=N absorbe entre 1690 y 1640 cm⁻¹, con intensidad variable. En las oximas se añade un O–H ancho (~3250 cm⁻¹) y el estiramiento N–O (~945 cm⁻¹).',
    detect(c) {
      const out = [];
      const oximeC = c.atomSet('[CX3]=[NX2][OX2H]');
      for (const m of c.matchRaw('[CX3;!$(C([#7])=[#7])]=[NX2;!$(N=N);!$(N[#8;!H1])]')) {
        if (oximeC.has(m.atoms[0])) {
          out.push({
            atoms: m.atoms, variant: 'oxima', bands: [
              b(3250, 0.4, 200, 'f', 'Estiramiento O–H de oxima', WHY.hbond, { range: [3100, 3350], gas: { center: 3600, fwhm: 22, vib: 'Estiramiento O–H libre de oxima' } }),
              b(1665, 0.25, 25, 'm', 'Estiramiento C=N', 'C=N de oxima.', { range: [1640, 1690] }),
              b(945, 0.55, 30, 'f', 'Estiramiento N–O', 'Banda intensa característica de oximas.', { range: [930, 960] })
            ]
          });
        } else {
          out.push({ atoms: m.atoms, variant: 'imina', bands: [b(1655, 0.35, 25, 'v', 'Estiramiento C=N', 'Doble enlace C=N; menos polar que el C=O, por eso su intensidad es variable.', { range: [1630, 1690] })] });
        }
      }
      return out;
    }
  },
  {
    id: 'azo', name: 'Azo (N=N)', category: 'Nitrogenados',
    theory: 'El N=N apenas cambia el momento dipolar al vibrar (en compuestos trans simétricos es inactivo en IR); cuando se observa, es una banda débil ~1575 cm⁻¹.',
    detect(c) {
      return c.matchRaw('[#6][NX2]=[NX2][#6]').map((m) => ({ atoms: m.atoms, bands: [b(1575, 0.08, 25, 'd', 'Estiramiento N=N', 'Débil o ausente por la simetría del grupo.', { range: [1550, 1630] })] }));
    }
  },
  {
    id: 'nitroso', name: 'Nitroso', category: 'Nitrogenados',
    theory: 'El estiramiento N=O de los C-nitrosocompuestos aparece ~1550–1500 cm⁻¹.',
    detect(c) {
      return c.matchRaw('[#6][NX2]=[OX1]').map((m) => ({ atoms: m.atoms, bands: [b(1530, 0.7, 25, 'f', 'Estiramiento N=O', 'Grupo nitroso.', { range: [1500, 1600] })] }));
    }
  },

  // ------------------------------------------------ INSATURACIONES C–C
  {
    id: 'alkyne', name: 'Alquino', category: 'Hidrocarburos',
    theory: 'Los alquinos terminales se reconocen por dos bandas: el ≡C–H, agudo e intenso a ~3300 cm⁻¹ (carbono sp con 50 % de carácter s, enlace C–H corto y fuerte) y el C≡C débil a ~2120 cm⁻¹. En alquinos internos simétricos el C≡C no cambia el dipolo y puede no observarse.',
    detect(c) {
      const out = [];
      for (const m of c.matchRaw('[CX2H1]#[CX2]')) {
        out.push({
          atoms: m.atoms, variant: 'terminal', bands: [
            b(3305, 0.45, 16, 'f', 'Estiramiento ≡C–H', 'Carbono sp: enlace C–H muy fuerte, banda aguda e intensa; se distingue del O–H/N–H por su estrechez.', { range: [3270, 3330] }),
            b(2120, 0.12, 14, 'd', 'Estiramiento C≡C', 'Triple enlace terminal: débil porque el cambio de dipolo es pequeño.', { range: [2100, 2140] }),
            b(635, 0.5, 30, 'f', 'Flexión ≡C–H', 'Flexión del ≡C–H, intensa y algo ancha.', { range: [610, 700] }),
            b(1250, 0.06, 40, 'md', 'Sobretono 2 × δ ≡C–H', 'Sobretono de la flexión ≡C–H.', { overtone: true, range: [1220, 1280] })
          ]
        });
      }
      for (const m of c.matchRaw('[#6][CX2H0]#[CX2H0][#6]')) {
        const sym = c.atoms[m.atoms[0]].h === c.atoms[m.atoms[3]].h;
        out.push({ atoms: m.atoms, variant: 'interno', bands: [b(2230, sym ? 0.02 : 0.06, 14, 'md', 'Estiramiento C≡C (alquino interno)', sym ? 'Alquino interno casi simétrico: el estiramiento apenas cambia el dipolo y la banda es muy débil o no se observa.' : 'Alquino interno: banda débil.', { range: [2190, 2260] })] });
      }
      return out.filter((x, i, arr) => arr.findIndex((y) => y.atoms.slice().sort().join() === x.atoms.slice().sort().join()) === i);
    }
  },
  {
    id: 'alkene', name: 'Alqueno', category: 'Hidrocarburos',
    theory: 'Los =C–H de carbonos sp² absorben por ENCIMA de 3000 cm⁻¹ (3100–3010), a diferencia de los C–H alquílicos. El C=C aparece ~1680–1620 cm⁻¹ (débil si está muy sustituido o es simétrico). Las flexiones =C–H fuera del plano (1000–650 cm⁻¹) son intensas y revelan el patrón de sustitución: vinilo 990 + 910, vinilideno 890, trans 965, cis ~700, trisustituido ~820 cm⁻¹.',
    detect(c) {
      const conjSet = c.atomSet('[CX3]=[CX3]-[$([#6]=[#6,#7,#8]),$(c),$([#6]#*)]');
      const seen = new Set();
      const out = [];
      for (const m of c.matchRaw('[CX3;!a]=[CX3;!a]')) {
        const [a, d] = m.atoms;
        const key = [a, d].sort().join('-');
        if (seen.has(key)) continue;
        seen.add(key);
        const bond = c.bonds[m.bonds[0]];
        const ha = c.atoms[a].h;
        const hd = c.atoms[d].h;
        const hs = [ha, hd].sort((x, y) => y - x).join('');
        const conj = conjSet.has(a) || conjSet.has(d);
        const inRing = c.ringSize[a] && c.ringSize[d] && c.ringSize[a] <= 8;
        const cc = (nu, eps, label, why) => b(conj ? nu - 25 : nu, conj ? eps * 2.5 : eps, 18, conj ? 'm' : label, 'Estiramiento C=C', why + (conj ? ' La conjugación baja su frecuencia y aumenta su intensidad.' : ''), { range: conj ? [1590, 1640] : [nu - 20, nu + 20] });
        let bands = [];
        let variant;
        if (hs === '21') {
          variant = 'monosustituido (vinilo)';
          bands = [
            b(3080, 0.12, 20, 'm', 'Estiramiento =C–H (=CH₂ asimétrico)', 'Carbono sp²: más carácter s que en sp³, enlace C–H más fuerte, por eso por encima de 3000 cm⁻¹.', { range: [3070, 3100] }),
            b(3010, 0.08, 20, 'd', 'Estiramiento =C–H', 'Segundo =C–H del vinilo.', { range: [3000, 3030] }),
            cc(1642, 0.18, 'm', 'Alqueno monosustituido.'),
            b(990, 0.45, 18, 'f', 'Flexión =C–H fuera del plano (trans-CH=CH₂)', 'Una de las dos bandas que identifican al grupo vinilo.', { range: [985, 1000] }),
            b(910, 0.55, 18, 'f', 'Flexión =CH₂ fuera del plano (aleteo)', 'Segunda banda diagnóstica del grupo vinilo.', { range: [905, 920] }),
            b(1825, 0.04, 30, 'md', 'Sobretono 2 × 910', 'Sobretono del aleteo =CH₂.', { overtone: true })
          ];
        } else if (hs === '20') {
          variant = 'disustituido 1,1 (vinilideno)';
          bands = [
            b(3080, 0.12, 20, 'm', 'Estiramiento =CH₂', 'Carbono sp² terminal.', { range: [3070, 3095] }),
            cc(1655, 0.16, 'm', 'Alqueno 1,1-disustituido.'),
            b(890, 0.6, 18, 'f', 'Aleteo =CH₂ fuera del plano', 'Banda intensa característica del grupo vinilideno (R₂C=CH₂).', { range: [880, 900] }),
            b(1780, 0.04, 30, 'md', 'Sobretono 2 × 890', 'Sobretono del aleteo =CH₂.', { overtone: true })
          ];
        } else if (hs === '11') {
          const cis = inRing || bond.stereo === 'cis' || bond.stereo === 'Z';
          const known = inRing || ['cis', 'trans', 'Z', 'E'].includes(bond.stereo);
          variant = cis ? 'disustituido cis (Z)' : `disustituido trans (E)${known ? '' : ' — estereoquímica no definida, se supone trans'}`;
          bands = [b(3020, 0.1, 20, 'm', 'Estiramiento =C–H', 'Carbono sp².', { range: [3010, 3040] })];
          if (cis) {
            bands.push(cc(1655, 0.12, 'm', 'Alqueno cis-disustituido.'));
            bands.push(b(700, 0.45, 40, 'f', 'Flexión =C–H fuera del plano (cis)', 'Banda ancha característica de alquenos cis (730–665 cm⁻¹).', { range: [665, 730] }));
          } else {
            bands.push(cc(1670, 0.03, 'md', 'Alqueno trans-disustituido: casi simétrico, el C=C es muy débil o no se observa.'));
            bands.push(b(965, 0.6, 18, 'f', 'Flexión =C–H fuera del plano (trans)', 'Banda intensa y muy diagnóstica de alquenos trans.', { range: [960, 980] }));
          }
        } else if (hs === '10') {
          variant = 'trisustituido';
          bands = [
            b(3020, 0.06, 20, 'd', 'Estiramiento =C–H', 'Un solo =C–H.', { range: [3010, 3040] }),
            cc(1670, 0.08, 'd', 'Alqueno trisustituido.'),
            b(820, 0.35, 25, 'm', 'Flexión =C–H fuera del plano', 'Alqueno trisustituido (840–800 cm⁻¹).', { range: [800, 840] })
          ];
        } else if (hs === '00') {
          variant = 'tetrasustituido';
          bands = [cc(1670, 0.02, 'md', 'Alqueno tetrasustituido: sin =C–H y con un C=C casi inactivo.')];
        } else {
          variant = 'eteno';
          bands = [b(3100, 0.1, 20, 'm', 'Estiramiento =C–H', 'Carbono sp².'), b(950, 0.6, 20, 'f', 'Aleteo =CH₂', 'Flexión fuera del plano.')];
        }
        out.push({ atoms: [a, d], bands, variant: variant + (conj ? ', conjugado' : '') + (inRing ? ', en anillo' : '') });
      }
      return out;
    }
  },
  {
    id: 'aromatic', name: 'Anillo bencénico', category: 'Aromáticos',
    theory: 'Los anillos bencénicos muestran: ν C–H aromático débil (3100–3000 cm⁻¹); un grupo de bandas débiles de sobretonos y combinaciones (2000–1667 cm⁻¹) cuya forma depende de la sustitución; estiramientos C=C del anillo a ~1600, ~1580, ~1500 y ~1450 cm⁻¹; y flexiones C–H fuera del plano (900–675 cm⁻¹), intensas, que dependen del número de hidrógenos ADYACENTES y permiten deducir el patrón de sustitución.',
    detect(c) {
      const out = [];
      const subs = {
        '5': ['monosustituido', [[750, 0.6, 'Flexión C–H fuera del plano (5 H adyacentes)', 'Banda intensa de bencenos monosustituidos (770–730 cm⁻¹).', [730, 770]], [695, 0.55, 'Deformación del anillo fuera del plano', 'Segunda banda diagnóstica del benceno monosustituido (710–690 cm⁻¹).', [690, 710]]]],
        '4': ['orto-disustituido (1,2)', [[750, 0.6, 'Flexión C–H fuera del plano (4 H adyacentes)', 'Una sola banda intensa en 770–735 cm⁻¹: típico de sustitución orto.', [735, 770]]]],
        '1,3': ['meta-disustituido (1,3)', [[780, 0.5, 'Flexión C–H fuera del plano (3 H adyacentes)', 'Banda de los tres H contiguos del anillo meta-disustituido.', [750, 810]], [880, 0.3, 'Flexión C–H fuera del plano (H aislado)', 'El H aislado entre los dos sustituyentes absorbe en 900–860 cm⁻¹.', [860, 900]], [690, 0.45, 'Deformación del anillo fuera del plano', 'Tercera banda típica de sustitución meta.', [680, 710]]]],
        '2,2': ['para-disustituido (1,4)', [[830, 0.6, 'Flexión C–H fuera del plano (2 H adyacentes)', 'Una banda intensa en 860–800 cm⁻¹: típico de sustitución para.', [800, 860]]]],
        '3': ['1,2,3-trisustituido', [[780, 0.5, 'Flexión C–H fuera del plano (3 H adyacentes)', 'Tres H contiguos.', [760, 800]], [710, 0.35, 'Deformación del anillo', 'Banda adicional del patrón 1,2,3.', [685, 720]]]],
        '1,2': ['1,2,4-trisustituido', [[820, 0.5, 'Flexión C–H fuera del plano (2 H adyacentes)', 'Dos H contiguos.', [800, 860]], [880, 0.3, 'Flexión C–H fuera del plano (H aislado)', 'H aislado.', [860, 900]]]],
        '1,1,1': ['1,3,5-trisustituido', [[850, 0.45, 'Flexión C–H fuera del plano (H aislados)', 'Tres H aislados.', [810, 865]], [690, 0.4, 'Deformación del anillo', 'Banda adicional del patrón 1,3,5.', [675, 710]]]],
        '2': ['1,2,3,4-tetrasustituido', [[810, 0.45, 'Flexión C–H fuera del plano (2 H adyacentes)', 'Dos H contiguos.', [800, 860]]]],
        '1,1': ['tetrasustituido (1,2,3,5 / 1,2,4,5)', [[870, 0.4, 'Flexión C–H fuera del plano (H aislados)', 'H aislados.', [850, 900]]]],
        '1': ['pentasustituido', [[870, 0.35, 'Flexión C–H fuera del plano (H aislado)', 'Un solo H en el anillo.', [850, 900]]]],
        '': ['hexasustituido', []]
      };
      const rings = c.matchRaw('c1ccccc1');
      const membership = new Map();
      for (const r of rings) for (const a of r.atoms) membership.set(a, (membership.get(a) || 0) + 1);
      for (const m of rings) {
        const hasH = m.atoms.map((i) => c.atoms[i].h > 0);
        const fused = m.atoms.filter((a) => membership.get(a) > 1).length;
        const nH = hasH.filter(Boolean).length;
        // longitudes de corridas de H adyacentes (anillo circular)
        let runs = [];
        if (nH === 6) runs = [6];
        else if (nH > 0) {
          const start = hasH.findIndex((x) => !x);
          let run = 0;
          for (let k = 1; k <= 6; k++) {
            const v = hasH[(start + k) % 6];
            if (v) run++;
            else if (run) { runs.push(run); run = 0; }
          }
          if (run) runs.push(run);
        }
        runs.sort((x, y) => x - y);
        const key = nH === 6 ? 'benceno' : runs.join(',');
        let [variant, oop] = key === 'benceno'
          ? ['benceno sin sustituir', [[675, 0.7, 'Flexión C–H fuera del plano', 'Benceno: banda única muy intensa (~675 cm⁻¹).', [660, 690]]]]
          : subs[key] || ['polisustituido', []];
        if (fused >= 2 && key === '4') {
          variant = 'anillo fusionado (tipo naftaleno)';
          oop = [[780, 0.65, 'Flexión C–H fuera del plano (4 H adyacentes, anillo fusionado)', 'En naftalenos el anillo con cuatro H contiguos da una banda muy intensa cerca de 780 cm⁻¹.', [760, 800]]];
        } else if (fused >= 2) {
          variant += ', fusionado';
        }
        const bands = [
          b(3060, 0.035 * Math.max(1, nH), 25, 'd', 'Estiramiento C–H aromático', 'Carbono sp² aromático: los C–H absorben justo por encima de 3000 cm⁻¹ (3100–3000), con varias bandas débiles.', { range: [3000, 3100] }),
          b(3030, 0.03 * Math.max(1, nH), 22, 'd', 'Estiramiento C–H aromático', 'Componente adicional del ν C–H aromático.', { range: [3000, 3100] }),
          b(1600, 0.22, 14, 'm', 'Estiramiento C=C del anillo (8a)', 'Vibración de estiramiento del esqueleto aromático; su intensidad aumenta si hay un sustituyente conjugado o polar.', { range: [1585, 1615] }),
          b(1580, 0.12, 12, 'v', 'Estiramiento C=C del anillo (8b)', 'Componente degenerada de la anterior; se hace visible con sustituyentes conjugados.', { range: [1570, 1590] }),
          b(1495, 0.32, 12, 'm', 'Estiramiento C=C del anillo (19a)', 'Banda del esqueleto aromático, a menudo la más intensa del grupo.', { range: [1480, 1520] }),
          b(1455, 0.25, 14, 'm', 'Estiramiento C=C del anillo (19b)', 'Puede solaparse con la flexión de CH₂/CH₃.', { range: [1440, 1465] }),
          b(1030, 0.12, 15, 'd', 'Flexión C–H en el plano', 'Vibraciones C–H en el plano del anillo (1300–1000 cm⁻¹), débiles.', { range: [1000, 1300] }),
          b(1950, 0.03, 40, 'md', 'Sobretonos y combinaciones (2000–1667)', 'Bandas débiles de sobretonos de las flexiones fuera del plano; su patrón ("dedos de benceno") también refleja la sustitución.', { overtone: true, range: [1667, 2000] }),
          b(1810, 0.025, 40, 'md', 'Sobretonos y combinaciones (2000–1667)', 'Segundo componente del patrón de sobretonos aromáticos.', { overtone: true, range: [1667, 2000] })
        ];
        for (const [nu, eps, vib, why, range] of oop) bands.push(b(nu, eps, 22, 'f', vib, why, { range }));
        out.push({ atoms: m.atoms, bands, variant });
      }
      return out;
    }
  },
  {
    id: 'heteroaromatic', name: 'Heterociclo aromático', category: 'Aromáticos',
    theory: 'Los anillos heteroaromáticos (piridina, pirrol, furano, tiofeno…) muestran estiramientos C=C/C=N del anillo entre 1600 y 1400 cm⁻¹, C–H aromáticos ~3100–3000 cm⁻¹ y flexiones C–H fuera del plano intensas. Pirrol e indol presentan además un N–H ~3400 cm⁻¹.',
    detect(c) {
      const out = [];
      const defs = [
        ['n1ccccc1', 'piridina', [[1590, 0.35, 'Estiramiento C=C / C=N del anillo', [1570, 1610]], [1575, 0.3, 'Estiramiento C=C / C=N del anillo', [1560, 1590]], [1480, 0.3, 'Estiramiento del anillo', [1465, 1500]], [1435, 0.45, 'Estiramiento del anillo', [1420, 1450]], [750, 0.5, 'Flexión C–H fuera del plano', [700, 780]]]],
        ['[nH]1cccc1', 'pirrol / indol', [[3400, 0.4, 'Estiramiento N–H (heteroaromático)', [3350, 3480]], [1530, 0.25, 'Estiramiento del anillo', [1500, 1560]], [1470, 0.25, 'Estiramiento del anillo', [1450, 1490]], [740, 0.5, 'Flexión C–H fuera del plano', [720, 770]]]],
        ['o1cccc1', 'furano', [[1580, 0.2, 'Estiramiento C=C del anillo', [1560, 1600]], [1500, 0.25, 'Estiramiento C=C del anillo', [1480, 1510]], [1015, 0.5, 'Estiramiento C–O–C del anillo / respiración', [1000, 1030]], [740, 0.6, 'Flexión C–H fuera del plano', [720, 760]]]],
        ['s1cccc1', 'tiofeno', [[1520, 0.2, 'Estiramiento C=C del anillo', [1500, 1540]], [1410, 0.3, 'Estiramiento del anillo', [1400, 1440]], [700, 0.6, 'Flexión C–H fuera del plano', [690, 720]]]]
      ];
      for (const [s, variant, list] of defs) {
        for (const m of c.matchRaw(s)) {
          const bands = [b(3080, 0.06, 25, 'd', 'Estiramiento C–H heteroaromático', 'C–H sobre carbono sp² aromático.', { range: [3000, 3150] })];
          for (const [nu, eps, vib, range] of list) {
            bands.push(b(nu, eps, nu > 3000 ? 60 : 18, nu > 3000 ? 'm' : 'm', vib, nu > 3000 ? 'El N–H del pirrol forma puentes de hidrógeno más débiles que un O–H; banda relativamente aguda.' : `Vibración característica del anillo de ${variant}.`, { range, gas: nu > 3000 ? { center: 3530, fwhm: 20 } : null }));
          }
          out.push({ atoms: m.atoms, bands, variant });
        }
      }
      return out;
    }
  },

  // ------------------------------------------------ C–H ALQUÍLICOS
  {
    id: 'alkyl', name: 'C–H alquílicos (sp³)', category: 'Hidrocarburos',
    theory: 'Los C–H de carbonos sp³ absorben por DEBAJO de 3000 cm⁻¹ (2960–2850). El metilo da un estiramiento asimétrico ~2962 y uno simétrico ~2872 cm⁻¹; el metileno ~2926 y ~2853 cm⁻¹. En la región de flexión: ~1465 cm⁻¹ (tijera de CH₂ y flexión asimétrica de CH₃) y ~1375 cm⁻¹ ("paraguas" del CH₃). Un doblete 1385/1370 delata un isopropilo o gem-dimetilo; 1395/1365 un terc-butilo; una banda ~720 cm⁻¹ indica cadenas de cuatro o más CH₂.',
    detect(c) {
      const ch3 = c.uniqueAtoms('[CX4H3;!$(C[O,#7])]');
      const ch3O = c.uniqueAtoms('[CX4H3][OX2]');
      const ch3N = c.uniqueAtoms('[CX4H3][#7]');
      const ch2 = c.uniqueAtoms('[CX4H2]');
      const ch = c.uniqueAtoms('[CX4H1]');
      const acetyl = c.uniqueAtoms('[CX4H3][CX3]=O');
      const out = [];
      const nMe = ch3.length + ch3O.length + ch3N.length;
      if (nMe) {
        const bands = [
          b(2962, 0.11 * nMe, 18, 'f', 'Estiramiento C–H asimétrico del CH₃', 'Los tres C–H del metilo se estiran fuera de fase. Carbono sp³: por debajo de 3000 cm⁻¹.', { range: [2950, 2975] }),
          b(2872, 0.05 * Math.max(ch3.length, 0.3), 16, 'm', 'Estiramiento C–H simétrico del CH₃', 'Modo en fase de los tres C–H del metilo.', { range: [2860, 2885] }),
          b(1460, 0.05 * nMe, 18, 'm', 'Flexión asimétrica del CH₃', 'Deformación asimétrica; coincide con la tijera del CH₂.', { range: [1440, 1470] })
        ];
        const plainMe = ch3.length - acetyl.length;
        if (plainMe > 0) bands.push(b(1378, 0.045 * plainMe, 14, 'm', 'Flexión simétrica del CH₃ ("paraguas")', 'Banda característica del metilo, muy útil para distinguir CH₃ de CH₂.', { range: [1370, 1385] }));
        if (acetyl.length) bands.push(b(1360, 0.25 * acetyl.length, 16, 'f', 'Flexión simétrica del CH₃ unido a C=O', 'En metilcetonas y acetatos el "paraguas" del CH₃ se desplaza a ~1360 cm⁻¹ y se intensifica.', { range: [1350, 1370] }));
        out.push({ atoms: [...ch3, ...ch3O, ...ch3N], bands, variant: `${nMe} CH₃`, count: nMe });
      }
      if (ch3O.length) out.push({ atoms: ch3O, variant: 'O–CH₃', count: ch3O.length, bands: [b(2830, 0.09 * ch3O.length, 18, 'm', 'Estiramiento C–H simétrico del O–CH₃', 'El metoxilo tiene un estiramiento simétrico característico, más bajo (~2830 cm⁻¹) que el de un metilo unido a carbono.', { range: [2815, 2850] })] });
      if (ch3N.length) out.push({ atoms: ch3N, variant: 'N–CH₃', count: ch3N.length, bands: [b(2790, 0.08 * ch3N.length, 20, 'm', 'Estiramiento C–H simétrico del N–CH₃', 'El par libre del nitrógeno debilita los C–H del N-metilo (~2800–2780 cm⁻¹).', { range: [2770, 2820] })] });
      if (ch2.length) {
        const n = ch2.length;
        const alphaCO = c.uniqueAtoms('[CX4H2][CX3]=O').length;
        const bands = [
          b(2926, 0.1 * n, 18, 'f', 'Estiramiento C–H asimétrico del CH₂', 'Los dos C–H del metileno fuera de fase.', { range: [2915, 2940] }),
          b(2853, 0.07 * n, 16, 'm', 'Estiramiento C–H simétrico del CH₂', 'Modo en fase del metileno.', { range: [2845, 2865] }),
          b(1465, 0.05 * n, 16, 'm', 'Flexión de tijera del CH₂', 'Deformación angular H–C–H del metileno.', { range: [1455, 1475] })
        ];
        if (alphaCO) bands.push(b(1420, 0.08 * alphaCO, 16, 'm', 'Tijera del CH₂ en α a C=O', 'El metileno vecino a un carbonilo absorbe más bajo (~1420 cm⁻¹).', { range: [1410, 1430] }));
        if (c.has('[CH2X4]!@[CH2X4]!@[CH2X4]!@[CH2X4]')) bands.push(b(722, 0.1 + 0.015 * n, 14, 'm', 'Balanceo (rocking) del CH₂ en cadena', 'Cuatro o más CH₂ consecutivos en cadena abierta generan esta banda ~720 cm⁻¹ (desdoblada en sólidos cristalinos).', { range: [715, 730] }));
        if (c.has('[CH2r3]')) bands.push(b(3080, 0.1, 18, 'm', 'Estiramiento C–H del ciclopropano', 'La tensión del anillo de tres miembros da más carácter s a los C–H, que absorben por encima de 3000 cm⁻¹.', { range: [3040, 3100] }));
        out.push({ atoms: ch2, bands, variant: `${n} CH₂`, count: n });
      }
      if (ch.length) {
        out.push({ atoms: ch, variant: `${ch.length} C–H terciario`, count: ch.length, bands: [b(2890, 0.03 * ch.length, 20, 'd', 'Estiramiento C–H del metino', 'Banda débil, normalmente oculta por las de CH₃ y CH₂.', { range: [2880, 2900] }), b(1340, 0.03 * ch.length, 20, 'd', 'Flexión C–H del metino', 'Débil.', { range: [1330, 1350] })] });
      }
      const tbu = c.uniqueAtoms('[CX4]([CH3])([CH3])([CH3])', 0);
      const gem = c.uniqueAtoms('[CX4;!$(C([CH3])([CH3])[CH3])]([CH3])[CH3]', 0);
      if (tbu.length) out.push({ atoms: tbu, variant: 'terc-butilo', bands: [b(1395, 0.08, 10, 'm', 'Flexión simétrica CH₃ (doblete del terc-butilo)', 'El acoplamiento de los tres metilos desdobla el "paraguas" en 1395 (más débil) y 1365 cm⁻¹ (más intensa).', { range: [1385, 1400] }), b(1365, 0.14, 10, 'f', 'Flexión simétrica CH₃ (doblete del terc-butilo)', 'Componente intensa del doblete del terc-butilo.', { range: [1360, 1372] }), b(1250, 0.12, 20, 'm', 'Vibración del esqueleto C–(CH₃)₃', 'Banda esquelética del terc-butilo.', { range: [1240, 1260] })] });
      if (gem.length) out.push({ atoms: gem, variant: 'isopropilo / gem-dimetilo', bands: [b(1385, 0.07, 10, 'm', 'Flexión simétrica CH₃ (doblete gem-dimetilo)', 'Dos metilos sobre el mismo carbono acoplan sus "paraguas" y producen un doblete de intensidades parecidas (1385 y 1370 cm⁻¹).', { range: [1380, 1390] }), b(1368, 0.07, 10, 'm', 'Flexión simétrica CH₃ (doblete gem-dimetilo)', 'Segunda componente del doblete.', { range: [1365, 1375] }), b(1170, 0.06, 20, 'd', 'Vibración del esqueleto isopropilo', 'Banda esquelética.', { range: [1165, 1175] })] });
      return out;
    }
  },

  // ------------------------------------------------ HALÓGENOS
  {
    id: 'halide', name: 'Halogenuro de alquilo / arilo', category: 'Halogenados',
    theory: 'El estiramiento C–X baja al aumentar la masa del halógeno (ley de Hooke: ν ∝ √(k/μ)): C–F 1400–1000, C–Cl 850–550, C–Br 690–515 y C–I 600–500 cm⁻¹. Los C–F son muy intensos por la gran polaridad del enlace. Los grupos CH₂X dan además un aleteo CH₂ a 1300–1150 cm⁻¹.',
    detect(c) {
      const out = [];
      const cf3 = c.uniqueAtoms('[CX4](F)(F)F');
      for (const C of cf3) out.push({ atoms: [C], variant: 'CF₃', bands: [b(1320, 0.8, 30, 'f', 'Estiramiento C–F (CF₃)', 'El grupo CF₃ da varias bandas muy intensas por acoplamiento de los tres C–F.', { range: [1290, 1350] }), b(1170, 1.1, 35, 'mf', 'Estiramiento C–F asimétrico (CF₃)', 'Banda muy intensa.', { range: [1150, 1200] }), b(1130, 1.0, 30, 'mf', 'Estiramiento C–F simétrico (CF₃)', 'Banda muy intensa.', { range: [1100, 1150] })] });
      const defs = [
        ['[CX4;!$(C(F)(F)F)]F', 'C–F alifático', 1100, 0.9, 50, 'mf', [1000, 1400], 'El enlace C–F es muy polar: banda muy intensa en 1400–1000 cm⁻¹.'],
        ['cF', 'Ar–F', 1225, 0.8, 30, 'f', [1100, 1270], 'C(arilo)–F, acoplado con el anillo.'],
        ['[CX4]Cl', 'C–Cl alifático', 730, 0.6, 35, 'f', [550, 850], 'El Cl, más pesado que el F, baja la frecuencia del estiramiento C–X.'],
        ['cCl', 'Ar–Cl', 1090, 0.45, 20, 'm', [1035, 1100], 'En clorobencenos el C–Cl se acopla con el anillo y aparece ~1090 cm⁻¹.'],
        ['[CX4]Br', 'C–Br alifático', 600, 0.6, 35, 'f', [515, 690], 'El bromo, aún más pesado, desplaza la banda a 690–515 cm⁻¹.'],
        ['cBr', 'Ar–Br', 1070, 0.4, 20, 'm', [1030, 1075], 'C(arilo)–Br acoplado con el anillo.'],
        ['[#6]I', 'C–I', 520, 0.5, 35, 'f', [485, 600], 'El yodo es el más pesado: banda en el extremo de baja frecuencia.']
      ];
      for (const [s, variant, nu, eps, w, label, range, why] of defs) {
        const ms = c.matchRaw(s);
        if (!ms.length) continue;
        const bands = [b(nu, eps * ms.length, w, label, `Estiramiento ${variant}`, why, { range })];
        const ch2x = c.uniqueAtoms(s.replace('[CX4]', '[CH2X4]').replace('[CX4;!$(C(F)(F)F)]', '[CH2X4]'));
        if (ch2x.length && /Cl|Br|I/.test(variant) && s.startsWith('[CX4]')) {
          bands.push(b(variant.includes('Cl') ? 1265 : 1230, 0.25 * ch2x.length, 25, 'm', 'Aleteo del CH₂–X', 'Los metilenos unidos a halógeno (CH₂X) presentan un aleteo característico entre 1300 y 1150 cm⁻¹.', { range: [1150, 1300] }));
        }
        out.push({ atoms: ms.flatMap((m) => m.atoms), variant, count: ms.length, bands });
      }
      return out;
    }
  },

  // ------------------------------------------------ AZUFRE, FÓSFORO, SILICIO
  {
    id: 'sulfur', name: 'Grupos con azufre', category: 'Azufre',
    theory: 'El S–H (tiol) absorbe débil ~2570 cm⁻¹ en una región casi vacía. Los grupos S=O son muy polares: sulfóxidos ~1050, sulfonas 1350–1300 y 1160–1120, sulfonamidas ~1340 y ~1160, ésteres sulfonato ~1360 y ~1175 cm⁻¹ (bandas muy intensas).',
    detect(c) {
      const out = [];
      const defs = [
        ['[#6][SX2H]', 'tiol', [[2570, 0.1, 15, 'd', 'Estiramiento S–H', 'El S es pesado y el enlace S–H débil y poco polar: banda débil pero muy diagnóstica.', [2550, 2600]], [700, 0.1, 30, 'd', 'Estiramiento C–S', 'Débil.', [600, 750]]]],
        ['[#6][SX2][#6]', 'sulfuro', [[690, 0.12, 30, 'd', 'Estiramiento C–S', 'Banda débil, poco útil para el diagnóstico.', [600, 750]]]],
        ['[#6][SX2][SX2][#6]', 'disulfuro', [[510, 0.1, 20, 'd', 'Estiramiento S–S', 'Débil (S–S casi apolar).', [490, 540]]]],
        ['[#6][SX3](=O)[#6]', 'sulfóxido', [[1050, 0.9, 30, 'f', 'Estiramiento S=O', 'Enlace S=O muy polar: banda intensa ~1050 cm⁻¹.', [1030, 1070]]]],
        ['[#6][SX4](=O)(=O)[#6]', 'sulfona', [[1310, 0.9, 30, 'f', 'Estiramiento asimétrico SO₂', 'Los dos S=O acoplados dan dos bandas.', [1300, 1350]], [1150, 0.9, 30, 'f', 'Estiramiento simétrico SO₂', 'Componente simétrica.', [1120, 1160]]]],
        ['[SX4](=O)(=O)[NX3]', 'sulfonamida', [[1340, 0.9, 30, 'f', 'Estiramiento asimétrico SO₂', 'Sulfonamida.', [1310, 1370]], [1160, 0.9, 30, 'f', 'Estiramiento simétrico SO₂', 'Sulfonamida.', [1140, 1180]]]],
        ['[SX4](=O)(=O)[OX2][#6]', 'éster sulfonato', [[1360, 0.9, 30, 'f', 'Estiramiento asimétrico SO₂', 'Sulfonato (p. ej., tosilatos).', [1340, 1375]], [1175, 0.9, 30, 'f', 'Estiramiento simétrico SO₂', 'Sulfonato.', [1165, 1190]], [950, 0.5, 30, 'f', 'Estiramiento S–O–C', 'Banda del enlace S–O–C.', [900, 1000]]]],
        ['[SX4](=O)(=O)[OX2H]', 'ácido sulfónico', [[3000, 0.3, 450, 'f', 'Estiramiento O–H (asociado)', WHY.hbond, [2500, 3300]], [1350, 0.9, 30, 'f', 'Estiramiento asimétrico SO₂', 'Ácido sulfónico.', [1340, 1360]], [1160, 0.9, 30, 'f', 'Estiramiento simétrico SO₂', 'Ácido sulfónico.', [1150, 1180]]]],
        ['[SX4](=O)(=O)Cl', 'cloruro de sulfonilo', [[1375, 0.9, 25, 'f', 'Estiramiento asimétrico SO₂', 'Cloruro de sulfonilo.', [1365, 1390]], [1185, 0.9, 25, 'f', 'Estiramiento simétrico SO₂', 'Cloruro de sulfonilo.', [1170, 1195]]]]
      ];
      for (const [s, variant, list] of defs) {
        const ms = c.matchRaw(s);
        if (!ms.length) continue;
        const bands = list.map(([nu, eps, w, label, vib, why, range]) => b(nu, eps * ms.length, w, label, vib, why, { range }));
        if (variant === 'sulfonamida') {
          const nh = ms.reduce((s2, m) => s2 + c.atoms[m.atoms[3]].h, 0);
          if (nh) bands.push(b(3270, 0.3, 60, 'm', 'Estiramiento N–H (sulfonamida)', 'N–H ácido de la sulfonamida, asociado.', { range: [3200, 3390] }));
        }
        out.push({ atoms: ms.flatMap((m) => m.atoms), variant, count: ms.length, bands });
      }
      return out;
    }
  },
  {
    id: 'phosphorus', name: 'Grupos con fósforo', category: 'Otros heteroátomos',
    theory: 'El enlace P=O da una banda muy intensa entre 1300 y 1150 cm⁻¹; los enlaces P–O–C (alquilo) absorben muy fuerte entre 1050 y 1000 cm⁻¹; el P–H aparece ~2400 cm⁻¹.',
    detect(c) {
      const out = [];
      const po = c.matchRaw('[PX4]=O');
      if (po.length) out.push({ atoms: po.flatMap((m) => m.atoms), variant: 'P=O', bands: [b(1260, 0.9 * po.length, 35, 'f', 'Estiramiento P=O', 'Enlace P=O muy polar.', { range: [1150, 1300] })] });
      const poc = c.matchRaw('[PX4,PX3][OX2][CX4]');
      if (poc.length) out.push({ atoms: poc.flatMap((m) => m.atoms), variant: 'P–O–C', bands: [b(1030, 0.9 * Math.min(3, poc.length), 40, 'mf', 'Estiramiento P–O–C', 'Banda muy intensa de ésteres de fósforo.', { range: [970, 1050] })] });
      const ph = c.matchRaw('[P;!H0]');
      if (ph.length) out.push({ atoms: ph.flatMap((m) => m.atoms), variant: 'P–H', bands: [b(2400, 0.3, 25, 'm', 'Estiramiento P–H', 'Banda aguda en una región libre.', { range: [2350, 2440] })] });
      return out;
    }
  },
  {
    id: 'silicon', name: 'Grupos con silicio', category: 'Otros heteroátomos',
    theory: 'El Si–CH₃ da dos bandas muy características (~1260 aguda e intensa, y ~800 cm⁻¹); Si–O–C y Si–O–Si dan una banda ancha y muy intensa 1100–1000 cm⁻¹; el Si–H aparece ~2150 cm⁻¹.',
    detect(c) {
      const out = [];
      const sime = c.matchRaw('[Si][CH3]');
      if (sime.length) out.push({ atoms: sime.flatMap((m) => m.atoms), variant: 'Si–CH₃', bands: [b(1260, 0.6 * Math.min(3, sime.length), 12, 'f', 'Flexión simétrica CH₃ del Si–CH₃', 'Banda aguda e intensa típica de grupos sililo.', { range: [1245, 1265] }), b(800, 0.6 * Math.min(3, sime.length), 35, 'f', 'Balanceo CH₃ / estiramiento Si–C', 'Segunda banda característica del Si–CH₃.', { range: [750, 865] })] });
      const sio = c.matchRaw('[Si][OX2]');
      if (sio.length) out.push({ atoms: sio.flatMap((m) => m.atoms), variant: 'Si–O', bands: [b(1080, 1.0, 60, 'mf', 'Estiramiento Si–O–C / Si–O–Si', 'Banda ancha y muy intensa.', { range: [1000, 1110] })] });
      const sih = c.matchRaw('[Si;!H0]');
      if (sih.length) out.push({ atoms: sih.flatMap((m) => m.atoms), variant: 'Si–H', bands: [b(2150, 0.6, 20, 'f', 'Estiramiento Si–H', 'Banda intensa en la región de triples enlaces.', { range: [2100, 2260] })] });
      return out;
    }
  },
  {
    id: 'boron', name: 'Grupos con boro', category: 'Otros heteroátomos',
    theory: 'Los ácidos y ésteres borónicos muestran un estiramiento B–O muy intenso ~1350 cm⁻¹; los ácidos borónicos, además, un O–H ancho ~3300 cm⁻¹.',
    detect(c) {
      const bo = c.matchRaw('[BX3][OX2]');
      if (!bo.length) return [];
      const bands = [b(1350, 0.9, 35, 'f', 'Estiramiento B–O', 'Enlace B–O muy polar.', { range: [1310, 1380] })];
      if (c.has('[BX3][OX2H]')) bands.push(b(3280, 0.4, 220, 'f', 'Estiramiento O–H (ácido borónico)', WHY.hbond, { range: [3200, 3400] }));
      return [{ atoms: bo.flatMap((m) => m.atoms), variant: 'B–O', bands }];
    }
  }
];

// ------------------------------------------------ TABLA DE CORRELACIÓN GENERAL
// Se usa para interpretar espectros sin estructura conocida.
export const CORRELATION = [
  { lo: 3650, hi: 3580, shape: 'aguda', label: 'f', vib: 'O–H libre (sin puentes de H)', groups: ['alcohol', 'fenol'], note: 'Alcoholes/fenoles en gas o en disolución diluida.' },
  { lo: 3550, hi: 3200, shape: 'ancha', label: 'f', vib: 'O–H asociado por puentes de H', groups: ['alcohol', 'fenol'], note: 'Banda ancha e intensa típica de alcoholes y fenoles.' },
  { lo: 3500, hi: 3300, shape: 'media', label: 'm', vib: 'N–H (aminas, amidas)', groups: ['amina', 'amida'], note: 'Dos bandas = NH₂; una banda = N–H secundario.' },
  { lo: 3330, hi: 3270, shape: 'aguda', label: 'f', vib: '≡C–H de alquino terminal', groups: ['alquino'], note: 'Banda aguda e intensa; confirmar con C≡C ~2120 cm⁻¹.' },
  { lo: 3300, hi: 2500, shape: 'muy ancha', label: 'f', vib: 'O–H de ácido carboxílico (dímero)', groups: ['ácido carboxílico'], note: 'Banda muy ancha que se solapa con los C–H; confirmar con C=O ~1710 cm⁻¹.' },
  { lo: 3100, hi: 3000, shape: 'aguda', label: 'd', vib: 'C–H de carbono sp² (aromático o alqueno)', groups: ['aromático', 'alqueno'], note: 'Por encima de 3000 cm⁻¹ = C–H insaturado.' },
  { lo: 3000, hi: 2840, shape: 'aguda', label: 'f', vib: 'C–H de carbono sp³ (CH₃, CH₂, CH)', groups: ['alcano'], note: 'Por debajo de 3000 cm⁻¹ = C–H saturado.' },
  { lo: 2850, hi: 2800, shape: 'aguda', label: 'm', vib: 'C–H de aldehído (Fermi) / O–CH₃', groups: ['aldehído', 'éter metílico'], note: 'Junto a ~2720 cm⁻¹ indica aldehído.' },
  { lo: 2775, hi: 2695, shape: 'aguda', label: 'm', vib: 'C–H de aldehído (doblete de Fermi)', groups: ['aldehído'], note: 'Muy diagnóstica de aldehídos.' },
  { lo: 2600, hi: 2550, shape: 'aguda', label: 'd', vib: 'S–H de tiol', groups: ['tiol'], note: 'Débil, en región vacía.' },
  { lo: 2380, hi: 2320, shape: 'aguda', label: 'v', vib: 'CO₂ atmosférico (artefacto)', groups: ['artefacto'], note: 'Suele deberse a CO₂ del aire mal compensado en el fondo.' },
  { lo: 2280, hi: 2240, shape: 'aguda', label: 'mf', vib: 'N=C=O (isocianato) / C≡N alifático', groups: ['isocianato', 'nitrilo'], note: 'Isocianato: muy intensa; nitrilo: media.' },
  { lo: 2240, hi: 2210, shape: 'aguda', label: 'f', vib: 'C≡N conjugado / C≡C interno', groups: ['nitrilo', 'alquino'], note: '' },
  { lo: 2170, hi: 2080, shape: 'media', label: 'f', vib: 'Azida, isotiocianato, C≡C terminal, Si–H', groups: ['acumulado', 'alquino'], note: 'Banda intensa → azida/isotiocianato; débil → C≡C.' },
  { lo: 2000, hi: 1667, shape: 'débiles', label: 'md', vib: 'Sobretonos y combinaciones aromáticas', groups: ['aromático'], note: 'Patrón de "dedos de benceno".' },
  { lo: 1870, hi: 1760, shape: 'aguda', label: 'f', vib: 'C=O de anhídrido, halogenuro de acilo o lactona tensionada', groups: ['anhídrido', 'halogenuro de acilo', 'lactona'], note: 'Dos bandas separadas ~60 cm⁻¹ → anhídrido.' },
  { lo: 1760, hi: 1730, shape: 'aguda', label: 'mf', vib: 'C=O de éster o aldehído alifático', groups: ['éster', 'aldehído'], note: 'Éster: con C–O intenso 1300–1000; aldehído: con 2720 cm⁻¹.' },
  { lo: 1730, hi: 1700, shape: 'aguda', label: 'mf', vib: 'C=O de cetona, ácido o aldehído conjugado', groups: ['cetona', 'ácido carboxílico', 'aldehído'], note: '' },
  { lo: 1700, hi: 1630, shape: 'aguda', label: 'mf', vib: 'C=O conjugado o de amida (amida I)', groups: ['cetona conjugada', 'amida'], note: '' },
  { lo: 1680, hi: 1620, shape: 'aguda', label: 'v', vib: 'C=C de alqueno / C=N', groups: ['alqueno', 'imina'], note: 'Más débil que un C=O.' },
  { lo: 1650, hi: 1580, shape: 'media', label: 'm', vib: 'Flexión N–H (NH₂) / C=C aromático', groups: ['amina', 'aromático'], note: '' },
  { lo: 1600, hi: 1450, shape: 'aguda', label: 'm', vib: 'C=C del anillo aromático', groups: ['aromático'], note: 'Bandas ~1600, ~1500 y ~1450 cm⁻¹.' },
  { lo: 1570, hi: 1500, shape: 'aguda', label: 'mf', vib: 'NO₂ asimétrico / amida II', groups: ['nitro', 'amida'], note: 'Nitro: con otra banda 1390–1330.' },
  { lo: 1475, hi: 1440, shape: 'aguda', label: 'm', vib: 'Flexión CH₂ / CH₃ asimétrica', groups: ['alcano'], note: '' },
  { lo: 1390, hi: 1365, shape: 'aguda', label: 'm', vib: 'Flexión simétrica CH₃ ("paraguas")', groups: ['alcano'], note: 'Doblete → isopropilo o terc-butilo.' },
  { lo: 1390, hi: 1330, shape: 'aguda', label: 'f', vib: 'NO₂ simétrico', groups: ['nitro'], note: '' },
  { lo: 1360, hi: 1250, shape: 'media', label: 'f', vib: 'C–N aromático / C–O de ácido / SO₂ asim.', groups: ['amina aromática', 'ácido', 'sulfona'], note: '' },
  { lo: 1300, hi: 1150, shape: 'media', label: 'f', vib: 'C–O de éster (C–C(=O)–O) / C–O aril éter', groups: ['éster', 'éter'], note: '' },
  { lo: 1400, hi: 1000, shape: 'media', label: 'mf', vib: 'C–F', groups: ['fluoruro'], note: '' },
  { lo: 1150, hi: 1000, shape: 'media', label: 'f', vib: 'C–O de alcohol / éter / éster (O–C–C)', groups: ['alcohol', 'éter', 'éster'], note: '1°: ~1050, 2°: ~1100, 3°: ~1150.' },
  { lo: 1000, hi: 650, shape: 'aguda', label: 'f', vib: 'Flexiones =C–H y C–H aromáticas fuera del plano', groups: ['alqueno', 'aromático'], note: 'Revelan el patrón de sustitución.' },
  { lo: 850, hi: 550, shape: 'media', label: 'f', vib: 'C–Cl', groups: ['cloruro'], note: '' },
  { lo: 730, hi: 715, shape: 'aguda', label: 'm', vib: 'Balanceo (CH₂)n, n ≥ 4', groups: ['alcano'], note: '' },
  { lo: 690, hi: 515, shape: 'media', label: 'f', vib: 'C–Br', groups: ['bromuro'], note: '' },
  { lo: 680, hi: 660, shape: 'aguda', label: 'v', vib: 'CO₂ atmosférico (flexión) / C–H aromático', groups: ['artefacto', 'aromático'], note: '' }
];
