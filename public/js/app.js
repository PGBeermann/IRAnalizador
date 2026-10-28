'use strict';
/* IRAnalizador — lógica del cliente */

const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];
const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

const EXAMPLES = [
  ['Etanol', 'CCO'], ['Acetona', 'CC(C)=O'], ['Ácido acético', 'CC(=O)O'], ['Acetato de etilo', 'CCOC(C)=O'],
  ['Benzaldehído', 'O=Cc1ccccc1'], ['Ácido benzoico', 'OC(=O)c1ccccc1'], ['Acetonitrilo', 'CC#N'], ['Anilina', 'Nc1ccccc1'],
  ['Tolueno', 'Cc1ccccc1'], ['Ciclohexanona', 'O=C1CCCCC1'], ['Nitrobenceno', 'O=[N+]([O-])c1ccccc1'], ['1-Hexeno', 'C=CCCCC'],
  ['1-Hexino', 'C#CCCCC'], ['Fenol', 'Oc1ccccc1'], ['Acetanilida', 'CC(=O)Nc1ccccc1'], ['Paracetamol', 'CC(=O)Nc1ccc(O)cc1'],
  ['Aspirina', 'CC(=O)Oc1ccccc1C(=O)O'], ['Benzoato de metilo', 'COC(=O)c1ccccc1'], ['Butanamida', 'CCCC(N)=O'], ['Hexano', 'CCCCCC']
];

const REGIONS = [
  { lo: 2500, hi: 4000, name: 'X–H', color: '59,130,246' },
  { lo: 2000, hi: 2500, name: 'Triples', color: '168,85,247' },
  { lo: 1500, hi: 2000, name: 'Dobles', color: '234,88,12' },
  { lo: 400, hi: 1500, name: 'Huella dactilar', color: '100,116,139' }
];

const state = {
  jsme: null,
  sim: null,        // resultado de /api/simulate
  exp: null,        // resultado de /api/analyze o /api/upload
  structure: null,  // SMILES usado
  yAxis: 'T',
  group: null,      // id del grupo resaltado
  hoverX: null,
  search: null
};

// ------------------------------------------------------------------ utilidades
function setStatus(msg, type = 'info') {
  const el = $('#status');
  el.className = `status ${type}`;
  el.textContent = msg || '';
}
async function api(path, opts = {}) {
  const res = await fetch(path, opts);
  let data = null;
  try { data = await res.json(); } catch { /* sin cuerpo */ }
  if (!res.ok) throw new Error(data?.error || `Error ${res.status}`);
  return data;
}
const postJSON = (path, body) => api(path, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) });
function busy(btn, on) {
  if (!btn) return;
  btn.disabled = on;
}
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

// ------------------------------------------------------------------ JSME
window.jsmeOnLoad = function () {
  const box = $('#jsme_container');
  box.innerHTML = '';
  state.jsme = new JSApplet.JSME('jsme_container', `${box.clientWidth}px`, `${box.clientHeight}px`, { options: 'nocanonize,zoom,newlook' });
  let t;
  new ResizeObserver(() => {
    clearTimeout(t);
    t = setTimeout(() => state.jsme.setSize(`${box.clientWidth}px`, `${box.clientHeight}px`), 120);
  }).observe(box);
  if (state.structure) state.jsme.readGenericMolecularInput(state.structure);
  else {
    state.jsme.readGenericMolecularInput('CC(C)=O');
    simulate('CC(C)=O');
  }
};
function editorSmiles() {
  try { return state.jsme ? state.jsme.smiles().trim() : ''; } catch { return ''; }
}
function setEditor(smiles) {
  if (state.jsme && smiles) {
    try { state.jsme.readGenericMolecularInput(smiles); } catch { /* ignorar */ }
  }
}

// ------------------------------------------------------------------ acciones
async function simulate(smiles) {
  const structure = smiles || editorSmiles();
  if (!structure) return setStatus('Dibuje una estructura en el editor o elija un ejemplo.', 'err');
  const btn = $('#btnSimulate');
  busy(btn, true);
  setStatus('Simulando espectro…', 'busy');
  try {
    const r = await postJSON('/api/simulate', { structure, phase: $('#phase').value });
    state.sim = r;
    state.exp = null;
    state.structure = r.molecule.smiles;
    state.group = null;
    setStatus(`Espectro simulado para ${r.molecule.formula} (${r.bands.length} bandas previstas, fase ${r.phase}).`);
    renderAll();
  } catch (e) {
    setStatus(e.message, 'err');
  } finally {
    busy(btn, false);
  }
}

async function search(q, type = 'auto') {
  const out = $('#searchResult');
  setStatus(`Consultando PubChem y NIST WebBook para «${q}»…`, 'busy');
  out.innerHTML = '';
  try {
    const r = await api(`/api/search?q=${encodeURIComponent(q)}&type=${encodeURIComponent(type)}`);
    state.search = r;
    renderSearch(r);
    setStatus(r.spectra.length
      ? `${r.compound.name}: ${r.spectra.length} espectro(s) IR disponibles en NIST. Elija uno para analizarlo.`
      : `${r.compound.name}: ${r.nistError || 'NIST no tiene espectros IR de este compuesto.'} Puede simular el espectro a partir de la estructura.`,
    r.spectra.length ? 'info' : (r.nistError ? 'err' : 'info'));
  } catch (e) {
    setStatus(e.message, 'err');
  }
}

async function analyzeNist(id, index, btn) {
  const c = state.search?.compound;
  busy(btn, true);
  setStatus('Descargando el espectro de NIST y asignando bandas…', 'busy');
  try {
    const r = await postJSON('/api/analyze', { source: { type: 'nist', id, index }, structure: c?.smiles || null });
    state.exp = r;
    state.sim = null;
    state.group = null;
    if (c?.smiles) { state.structure = c.smiles; setEditor(c.smiles); }
    $('#phase').value = r.phase;
    setStatus(`Espectro de ${c?.name || id} (${r.meta.state || 'fase no indicada'}) analizado: ${r.peaks.length} picos detectados.`);
    renderAll();
  } catch (e) {
    setStatus(e.message, 'err');
  } finally {
    busy(btn, false);
  }
}

async function uploadSpectrum(ev) {
  ev.preventDefault();
  const file = $('#fileInput').files[0];
  if (!file) return;
  const fd = new FormData();
  fd.append('file', file);
  fd.append('yMode', $('#yMode').value);
  fd.append('phase', $('#phase').value);
  const smi = $('#useStructure').checked ? editorSmiles() : '';
  if (smi) fd.append('structure', smi);
  const btn = ev.submitter;
  busy(btn, true);
  setStatus(`Analizando ${file.name}…`, 'busy');
  try {
    const r = await api('/api/upload', { method: 'POST', body: fd });
    state.exp = r;
    state.sim = null;
    state.group = null;
    state.structure = smi || null;
    setStatus(`${file.name}: ${r.peaks.length} picos detectados${smi ? ' y asignados con la estructura dibujada' : ' (sin estructura: interpretación por tabla de correlación)'}.`);
    renderAll();
  } catch (e) {
    setStatus(e.message, 'err');
  } finally {
    busy(btn, false);
  }
}

async function highlightGroup(id) {
  state.group = state.group === id ? null : id;
  const groups = currentGroups();
  const g = groups.find((x) => x.id === state.group);
  $$('#groups .chip').forEach((el) => el.classList.toggle('on', el.dataset.id === state.group));
  renderChart();
  const mol = currentMolecule();
  if (!mol) return;
  try {
    const r = await postJSON('/api/depict', { structure: mol.smiles, atoms: g ? g.atoms : [] });
    $('#molSvg').innerHTML = r.svg;
  } catch { /* ignorar */ }
}

// ------------------------------------------------------------------ datos actuales
const currentGroups = () => (state.exp ? state.exp.groups : state.sim?.groups) || [];
const currentMolecule = () => (state.exp ? state.exp.molecule : state.sim?.molecule) || null;
const simSpectrum = () => (state.exp ? state.exp.simulation?.spectrum : state.sim?.spectrum) || null;
const simBands = () => (state.exp ? state.exp.simulation?.bands : state.sim?.bands) || [];

// ------------------------------------------------------------------ render
function renderAll() {
  renderChart();
  renderMolecule();
  renderGroups();
  renderDiagnosis();
  renderTable();
  renderExplanation();
}

function valueAt(spec, x, key) {
  const i = Math.round((x - spec.x[0]) / 2);
  return spec[key][Math.max(0, Math.min(spec.x.length - 1, i))];
}

function renderChart() {
  const el = $('#chart');
  if (typeof Plotly === 'undefined') return;
  const sim = simSpectrum();
  const exp = state.exp?.spectrum;
  if (!sim && !exp) return;
  const key = state.yAxis;
  const isT = key === 'T';
  const traces = [];
  if (exp && $('#showExp').checked) {
    traces.push({
      x: exp.x, y: exp[key], name: 'Experimental', mode: 'lines', line: { color: css('--exp'), width: 1.5 },
      hovertemplate: '%{x} cm⁻¹ · %{y:.2f}<extra>Experimental</extra>'
    });
  }
  if (sim && $('#showSim').checked) {
    traces.push({
      x: sim.x, y: sim[key], name: exp ? 'Simulado (estructura)' : 'Simulado', mode: 'lines',
      line: { color: css('--sim'), width: exp ? 1.2 : 1.6, dash: exp ? 'dot' : 'solid' },
      opacity: exp ? 0.85 : 1,
      hovertemplate: '%{x} cm⁻¹ · %{y:.2f}<extra>Simulado</extra>'
    });
  }
  const shapes = [];
  const annotations = [];
  if ($('#showRegions').checked) {
    for (const r of REGIONS) {
      shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: r.lo, x1: r.hi, y0: 0, y1: 1, fillcolor: `rgba(${r.color},0.06)`, line: { width: 0 }, layer: 'below' });
      annotations.push({ x: (r.lo + r.hi) / 2, y: 1, xref: 'x', yref: 'paper', yanchor: 'bottom', text: r.name, showarrow: false, font: { size: 10, color: `rgb(${r.color})` } });
    }
  }
  // grupo resaltado
  const g = currentGroups().find((x) => x.id === state.group);
  if (g) {
    for (const bd of g.bands) {
      const w = Math.max(10, Math.min(bd.fwhm, 160) / 2);
      shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: bd.center - w, x1: bd.center + w, y0: 0, y1: 1, fillcolor: 'rgba(245,158,11,0.22)', line: { width: 0 }, layer: 'below' });
    }
  }
  if (state.hoverX) {
    shapes.push({ type: 'line', xref: 'x', yref: 'paper', x0: state.hoverX, x1: state.hoverX, y0: 0, y1: 1, line: { color: css('--hl'), width: 2, dash: 'dash' } });
  }
  // etiquetas de picos
  if ($('#showLabels').checked) {
    let pts = [];
    if (exp && $('#showExp').checked) {
      pts = state.exp.peaks.filter((p) => !p.artifact).sort((a, b) => b.height - a.height).slice(0, 16)
        .map((p) => ({ x: p.x, y: valueAt(exp, p.x, key), c: css('--exp') }));
    } else if (sim && $('#showSim').checked) {
      pts = simBands().filter((b) => !b.overtone && b.height > 0.12).sort((a, b) => b.height - a.height).slice(0, 16)
        .map((b) => ({ x: b.center, y: valueAt(sim, b.center, key), c: css('--sim') }));
    }
    for (const p of pts) {
      annotations.push({ x: p.x, y: p.y, text: String(Math.round(p.x)), showarrow: false, yshift: isT ? -9 : 9, textangle: -90, xanchor: 'center', yanchor: isT ? 'top' : 'bottom', font: { size: 10, color: p.c } });
    }
  }
  const text = css('--text');
  const grid = css('--border');
  const layout = {
    margin: { l: 58, r: 14, t: 26, b: 48 },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { family: 'system-ui, Segoe UI, sans-serif', size: 12, color: text },
    xaxis: { title: { text: 'Número de onda (cm⁻¹)' }, range: [4000, 400], autorange: false, gridcolor: grid, zeroline: false, tickvals: el.clientWidth < 560 ? [4000, 3000, 2000, 1500, 1000, 500] : [4000, 3500, 3000, 2500, 2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400] },
    yaxis: isT
      ? { title: { text: 'Transmitancia (%)' }, range: [0, 112], gridcolor: grid, zeroline: false }
      : { title: { text: 'Absorbancia (normalizada)' }, range: [-0.03, 1.18], gridcolor: grid, zeroline: false },
    legend: { orientation: 'h', x: 0, y: -0.2 },
    shapes, annotations, hovermode: 'x unified'
  };
  if (!el.dataset.plot) { el.innerHTML = ''; el.dataset.plot = '1'; }
  Plotly.react(el, traces, layout, {
    responsive: true, displaylogo: false, locale: 'es',
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d'],
    toImageButtonOptions: { filename: 'espectro_IR', scale: 2 }
  });
  // título y similitud
  const mol = currentMolecule();
  const name = state.exp?.meta?.title || mol?.formula || '';
  $('#chartTitle').textContent = state.exp ? `Espectro IR experimental${name ? ' · ' + name : ''}` : `Espectro IR simulado${name ? ' · ' + name : ''}`;
  const sb = $('#similarity');
  if (state.exp && state.exp.similarity != null) {
    sb.textContent = `Similitud simulado/experimental: ${state.exp.similarity} %`;
    sb.title = 'Similitud coseno entre ambos espectros en absorbancia (suavizados). Orientativa: la simulación por frecuencias de grupo no reproduce la huella dactilar fina.';
    sb.classList.remove('hidden');
  } else sb.classList.add('hidden');
}

function renderMolecule() {
  const mol = currentMolecule();
  const card = $('#molCard');
  const meta = state.exp?.meta;
  if (!mol && !meta) { card.hidden = true; return; }
  card.hidden = false;
  $('#molSvg').innerHTML = mol ? mol.svg : '<p class="muted" style="padding:1rem">Sin estructura asociada</p>';
  const rows = [];
  const c = state.exp && state.search?.compound;
  if (c && meta?.database?.startsWith('NIST')) rows.push(['Nombre', esc(c.name)]);
  if (mol) {
    rows.push(['Fórmula', esc(mol.formula)], ['Masa molar', `${mol.mw} g/mol`], ['IHD', esc(mol.ihd)], ['SMILES', `<code>${esc(mol.smiles)}</code>`]);
    if (mol.inchikey) rows.push(['InChIKey', `<code>${esc(mol.inchikey)}</code>`]);
  }
  const phase = state.exp ? state.exp.phase : state.sim?.phase;
  if (phase) rows.push(['Fase considerada', esc(phase)]);
  if (meta) {
    if (meta.database) rows.push(['Fuente', meta.url ? `<a href="${esc(meta.url)}" target="_blank" rel="noopener">${esc(meta.database)}</a>` : esc(meta.database)]);
    if (meta.state) rows.push(['Estado de la muestra', esc(meta.state.length > 90 ? meta.state.slice(0, 87) + '…' : meta.state)]);
    if (meta.origin) rows.push(['Origen', esc(meta.origin)]);
    if (meta.instrument) rows.push(['Instrumento', esc(meta.instrument)]);
    if (meta.resolution) rows.push(['Resolución', esc(meta.resolution) + ' cm⁻¹']);
    if (meta.cas) rows.push(['CAS', esc(meta.cas)]);
  }
  $('#molInfo').innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
}

function renderGroups() {
  const groups = currentGroups();
  const card = $('#groupsCard');
  if (!groups.length) { card.hidden = true; return; }
  card.hidden = false;
  $('#groups').innerHTML = groups.map((g) => `<button class="chip${g.id === state.group ? ' on' : ''}" data-id="${esc(g.id)}">${esc(g.name)}${g.variants.length ? `<small>${esc(g.variants.join(' · '))}</small>` : ''}</button>`).join('');
}

function renderDiagnosis() {
  const card = $('#diagCard');
  const d = state.exp?.diagnosis;
  if (!d || !d.length) { card.hidden = true; return; }
  card.hidden = false;
  $('#diagnosis').innerHTML = d.map((f) => `<li class="${f.ok ? '' : 'no'}">${esc(f.text)}</li>`).join('');
}

const intTag = (code, text) => `<span class="tag ${code === 'mf' || code === 'f' ? 'f' : code === 'm' ? 'm' : ''}">${esc(text)}</span>`;

function renderTable() {
  const card = $('#tableCard');
  const thead = $('#peakTable thead');
  const tbody = $('#peakTable tbody');
  if (state.exp) {
    card.hidden = false;
    const hasStruct = !!state.exp.molecule;
    $('#tableTitle').textContent = 'Picos observados y su asignación';
    $('#tableHint').textContent = hasStruct
      ? 'Cada pico experimental se compara con las bandas esperadas para la estructura. Pulse una fila para ver la explicación.'
      : 'Sin estructura: se proponen asignaciones con la tabla de correlación general. Pulse una fila para ver la explicación.';
    thead.innerHTML = '<tr><th>ν̃ obs. (cm⁻¹)</th><th>Intensidad</th><th>Forma</th><th>Asignación</th><th>Grupo</th><th>ν̃ esperado</th><th>Δ</th><th>Confianza</th></tr>';
    tbody.innerHTML = state.exp.peaks.map((p, i) => {
      const a = p.assignment;
      const t = state.exp.peakText[i];
      const alt = p.alternatives?.length ? `<p class="muted">Alternativas: ${p.alternatives.map((x) => `${esc(x.vib)} (${esc(x.groupName)}, ~${x.center} cm⁻¹)`).join('; ')}.</p>` : '';
      const label = a ? esc(a.vib) : p.artifact ? '<span class="muted">Artefacto (CO₂)</span>' : p.correlation?.length && p.x >= 1500 ? `<span class="muted">${esc(p.correlation[0].vib)}?</span>` : `<span class="muted">${p.x < 1500 ? 'Huella dactilar' : 'Sin asignar'}</span>`;
      return `<tr class="main" data-x="${p.x}" data-i="${i}">
        <td class="num">${p.x}</td><td>${intTag(p.intensity, p.intensityText)}</td><td>${esc(p.shape)}</td>
        <td>${label}</td><td>${a ? esc(a.groupName) : ''}</td><td class="num">${a ? a.center : ''}</td>
        <td class="num">${a ? (p.x - a.center >= 0 ? '+' : '') + (p.x - a.center) : ''}</td>
        <td>${a ? `<span class="conf" title="${Math.round(a.score * 100)} %"><i style="width:${Math.round(a.score * 100)}%"></i></span>` : ''}</td></tr>
        <tr class="detail" hidden><td colspan="8"><b>${esc(t.head)}.</b> ${esc(t.text)}${alt}</td></tr>`;
    }).join('');
  } else if (state.sim) {
    card.hidden = false;
    $('#tableTitle').textContent = 'Bandas simuladas';
    $('#tableHint').textContent = 'Bandas previstas por frecuencias de grupo. Pulse una fila para ver por qué aparece en esa posición.';
    thead.innerHTML = '<tr><th>ν̃ (cm⁻¹)</th><th>Intervalo típico</th><th>Intensidad</th><th>Vibración</th><th>Grupo</th></tr>';
    tbody.innerHTML = state.sim.bands.map((b) => `<tr class="main" data-x="${b.center}">
        <td class="num">${b.center}</td><td class="num muted">${b.range[1]}–${b.range[0]}</td><td>${intTag(b.label, b.labelText)}${b.overtone ? ' <span class="tag">sobretono</span>' : ''}</td>
        <td>${esc(b.vib)}${b.n > 1 ? ` <span class="muted">(×${b.n})</span>` : ''}</td><td>${esc(b.groupName)}</td></tr>
        <tr class="detail" hidden><td colspan="5">${esc(b.why)}</td></tr>`).join('');
  } else {
    card.hidden = true;
  }
}

function bandHTML(b) {
  return `<div class="band"><b class="nu">${b.center} cm⁻¹</b> <span class="muted">(${b.range[1]}–${b.range[0]} cm⁻¹, ${esc(b.labelText)})</span> — <b>${esc(b.vib)}</b>${b.groupName ? ` · ${esc(b.groupName)}` : ''}<br>${esc(b.why)}</div>`;
}

function renderExplanation() {
  const card = $('#explainCard');
  const out = [];
  const ex = state.exp ? state.exp.explanation : state.sim?.explanation;
  if (ex) {
    out.push('<h3>Resumen del compuesto</h3>');
    out.push(ex.summary.map((p) => `<p>${esc(p)}</p>`).join(''));
  }
  if (state.exp) {
    const e = state.exp;
    out.push('<h3>Interpretación del espectro experimental</h3>');
    const assigned = e.peaks.filter((p) => p.assignment).length;
    out.push(`<p>Se detectaron <b>${e.peaks.length}</b> picos (fase: ${esc(e.phase)}${e.meta?.state ? `; muestra: ${esc(e.meta.state)}` : ''}).${e.molecule ? ` <b>${assigned}</b> se asignaron a bandas esperadas para la estructura.` : ''}${e.similarity != null ? ` La similitud global con el espectro simulado es del ${e.similarity} %.` : ''}</p>`);
    if (e.phase === 'gas') out.push('<p class="note">Espectro en fase gaseosa: las moléculas no forman puentes de hidrógeno (O–H y N–H aparecen agudos y a mayor número de onda) y cada banda puede mostrar estructura rotovibracional (ramas P, Q y R), por lo que algunas vibraciones aparecen como dos o tres picos cercanos.</p>');
    for (const region of ['Estiramientos X–H', 'Triples enlaces y acumulados', 'Dobles enlaces', 'Huella dactilar']) {
      const items = e.peakText.filter((t) => t.region === region);
      if (!items.length) continue;
      out.push(`<h4>${esc(region)}</h4><ul>${items.map((t) => `<li><b>${esc(t.head)}.</b> ${esc(t.text)}</li>`).join('')}</ul>`);
    }
    if (e.missing?.length) {
      out.push('<h4>Bandas esperadas que no se observan</h4><ul>' + e.missing.map((m) => `<li>${esc(m.text)}</li>`).join('') + '</ul>');
    }
  } else if (ex) {
    out.push('<h3>Bandas previstas por región espectral</h3>');
    for (const r of ex.regions) {
      if (!r.bands.length) continue;
      out.push(`<h4>${esc(r.name)} (${r.hi}–${r.lo} cm⁻¹)</h4><p class="note">${esc(r.text)}</p>`);
      out.push(r.bands.map(bandHTML).join(''));
    }
  }
  if (ex) {
    out.push('<h3>Grupos funcionales y sus bandas características</h3>');
    for (const g of ex.groups) {
      out.push(`<h4>${esc(g.name)}${g.variants.length ? ` <span class="muted">(${esc(g.variants.join('; '))})</span>` : ''}</h4>`);
      out.push(`<p class="theory">${esc(g.theory)}</p>`);
      if (state.exp) out.push(`<ul>${g.bands.filter((b) => !b.overtone).map((b) => `<li><b>${b.center} cm⁻¹</b> (${b.range[1]}–${b.range[0]}): ${esc(b.vib)}</li>`).join('')}</ul>`);
    }
    if (ex.absences?.length) {
      out.push('<h3>Lo que NO debe aparecer (diagnóstico por ausencia)</h3><ul>' + ex.absences.map((a) => `<li>${esc(a)}</li>`).join('') + '</ul>');
    }
  }
  if (!ex && state.exp) {
    out.push('<p class="note">Para una asignación completa, dibuje la estructura propuesta en el editor y vuelva a analizar el archivo con la opción «Asignar picos usando la estructura dibujada».</p>');
  }
  card.hidden = !out.length;
  $('#explanation').innerHTML = out.join('');
}

function renderSearch(r) {
  const c = r.compound;
  const spectra = r.spectra || [];
  $('#searchResult').innerHTML = `
    <div class="compound">
      <h3>${esc(c.name)}</h3>
      <div class="meta">${esc(c.formula)} · ${c.mw ? esc(c.mw.toFixed ? c.mw.toFixed(2) : c.mw) + ' g/mol' : ''}${c.cas.length ? ' · CAS ' + esc(c.cas[0]) : ''}</div>
      ${c.iupac && c.iupac !== c.name ? `<div class="meta">IUPAC: ${esc(c.iupac)}</div>` : ''}
      <div class="row gap wrap" style="margin-top:.5rem">
        <button class="btn small" data-act="load">Cargar en el editor</button>
        <button class="btn small" data-act="sim">Simular espectro</button>
      </div>
      <div class="links">
        <a href="${esc(c.pubchemUrl)}" target="_blank" rel="noopener">PubChem CID ${esc(c.cid)}</a>
        ${r.nistUrl ? `<a href="${esc(r.nistUrl)}" target="_blank" rel="noopener">NIST WebBook</a>` : ''}
      </div>
      ${spectra.length ? `<ul class="speclist">${spectra.map((s) => `
        <li><span title="${esc(s.title)}"><b>#${s.index}</b> ${esc(s.title.length > 110 ? s.title.slice(0, 107) + "…" : s.title)}${s.phase ? ` <span class="tag">${esc(s.phase)}</span>` : ''}</span>
        <button class="btn small primary" data-act="nist" data-id="${esc(s.id)}" data-index="${s.index}">Analizar</button></li>`).join('')}</ul>`
        : `<p class="hint" style="margin-top:.6rem">${esc(r.nistError || 'NIST WebBook no tiene espectros IR para este compuesto.')}</p>`}
    </div>`;
}

// ------------------------------------------------------------------ eventos
function bind() {
  $('#examples').innerHTML = EXAMPLES.map(([n, s]) => `<button class="chip" data-smiles="${esc(s)}">${esc(n)}</button>`).join('');
  $('#examples').addEventListener('click', (e) => {
    const b = e.target.closest('[data-smiles]');
    if (!b) return;
    setEditor(b.dataset.smiles);
    simulate(b.dataset.smiles);
  });
  $('#btnSimulate').addEventListener('click', () => simulate());
  $('#btnClear').addEventListener('click', () => state.jsme && state.jsme.reset());
  $('#btnFindExp').addEventListener('click', () => {
    const s = editorSmiles();
    if (!s) return setStatus('Dibuje primero una estructura.', 'err');
    $('#searchInput').value = s;
    search(s, 'smiles');
  });
  $('#phase').addEventListener('change', () => { if (state.sim && !state.exp) simulate(state.structure); });
  $('#searchForm').addEventListener('submit', (e) => { e.preventDefault(); const q = $('#searchInput').value.trim(); if (q) search(q); });
  $('#searchResult').addEventListener('click', (e) => {
    const b = e.target.closest('[data-act]');
    if (!b) return;
    const c = state.search?.compound;
    if (b.dataset.act === 'load' && c) setEditor(c.smiles);
    if (b.dataset.act === 'sim' && c) { setEditor(c.smiles); simulate(c.smiles); }
    if (b.dataset.act === 'nist') analyzeNist(b.dataset.id, Number(b.dataset.index), b);
  });
  $('#uploadForm').addEventListener('submit', uploadSpectrum);
  ['#showSim', '#showExp', '#showLabels', '#showRegions'].forEach((s) => $(s).addEventListener('change', renderChart));
  $$('.seg button').forEach((b) => b.addEventListener('click', () => {
    state.yAxis = b.dataset.y;
    $$('.seg button').forEach((x) => x.classList.toggle('on', x === b));
    renderChart();
  }));
  $('#groups').addEventListener('click', (e) => { const b = e.target.closest('[data-id]'); if (b) highlightGroup(b.dataset.id); });
  const tbody = $('#peakTable tbody');
  tbody.addEventListener('click', (e) => {
    const tr = e.target.closest('tr.main');
    if (!tr) return;
    const d = tr.nextElementSibling;
    d.hidden = !d.hidden;
    tr.classList.toggle('on', !d.hidden);
  });
  tbody.addEventListener('mouseover', (e) => {
    const tr = e.target.closest('tr.main');
    const x = tr ? Number(tr.dataset.x) : null;
    if (x !== state.hoverX) { state.hoverX = x; renderChart(); }
  });
  tbody.addEventListener('mouseleave', () => { state.hoverX = null; renderChart(); });
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', renderChart);
}

document.addEventListener('DOMContentLoaded', bind);
