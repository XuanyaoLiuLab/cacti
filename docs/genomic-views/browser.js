import igv from './vendor/igv.esm.min.js';

const el = id => document.getElementById(id);
const pageBase = new URL('.', import.meta.url);
const trackBase = new URL('tracks/', pageBase);
const trackURL = file => new URL(`${state.locus}/${file}`, trackBase).href;
const examples = {fig2b_h3k27ac: 'Figure 2B', fig2d_h3k36me3: 'Figure 2D'};
const state = {locus: 'fig2b_h3k27ac', summary: null, browser: null, busy: false};
// Available for browser QA and console inspection.
window.cactiPreview = state;

async function readJSON(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`Could not read ${path} (HTTP ${response.status}).`);
  return response.json();
}

function setBusy(value) {
  state.busy = value;
  document.querySelectorAll('.examples button, .controls select, #reset').forEach(node => { node.disabled = value; });
  if (!value) el('genotype').disabled = el('view').value !== 'individuals';
}

function status(message, error = false) {
  el('status').textContent = message;
  el('status').classList.toggle('error', error);
}

function focalLocus(summary) {
  return `${summary.chrom}:${summary.cwindow_start_0based + 1}-${summary.cwindow_end_0based_exclusive}`;
}

function selectedLocus() {
  const summary = state.summary;
  const region = el('region').value;
  if (region === 'cwindow') return focalLocus(summary);
  if (region === 'compare') return [focalLocus(summary), ...summary.regions.filter(r => r.region_type === 'unrelated_control').map(r => r.browser_coordinates)];
  return summary.regions.find(r => r.region_id === region).browser_coordinates;
}

function updateDetails() {
  const summary = state.summary;
  el('locus-title').textContent = `${examples[state.locus]} · ${summary.mark}`;
  el('locus-meta').textContent = `Focal SNP ${summary.chrom}:${summary.snp_position_1based.toLocaleString('en-US')} · ${summary.n_profiles} profiles from ${summary.n_donors} donors`;
  el('legend').replaceChildren();
  const groups = summary.tracks.filter(t => t.kind === 'mean');
  const oldGroup = el('genotype').value;
  el('genotype').replaceChildren(new Option('All three groups', 'all'));
  for (const group of groups) {
    const item = document.createElement('span');
    item.className = 'legend-item';
    const swatch = document.createElement('span');
    swatch.className = 'swatch';
    swatch.style.background = group.color;
    item.append(swatch, `${group.genotype} · ${group.n_profiles} profiles`);
    el('legend').append(item);
    el('genotype').add(new Option(group.genotype, group.genotype));
  }
  if (groups.some(g => g.genotype === oldGroup)) el('genotype').value = oldGroup;
  el('samples-link').href = trackURL('sample_metadata.tsv');
  el('regions-link').href = trackURL('regions.tsv');
  document.querySelectorAll('[data-example]').forEach(button => {
    const selected = button.dataset.example === state.locus;
    button.classList.toggle('selected', selected);
    button.setAttribute('aria-pressed', String(selected));
  });
}

function localTrack(track) {
  // Resolve exported sessions against this copy of the website (local or live).
  const url = new URL(track.url);
  const file = decodeURIComponent(url.pathname).split(`/${state.locus}/`)[1];
  if (!file || file.includes('..')) throw new Error('Invalid regional track URL');
  const result = {...track, url: trackURL(file)};
  if (track.type === 'wig') {
    const info = state.summary.tracks.find(t => t.file === file);
    if (!info) throw new Error(`Track metadata missing for ${file}`);
    result.id = file;
    result.genotype = info.genotype;
    result.name = info.kind === 'mean' ? `${info.genotype} mean · ${info.n_profiles} profiles` : `${info.genotype} | ${info.sample_id} | rep ${info.replicate_number}`;
    result.height = info.kind === 'mean' ? 120 : 40;
    result.autoscale = info.kind === 'mean';
  } else {
    if (url.pathname.endsWith('focal_snp.bed')) { result.name = 'Focal SNP'; result.color = '#263640'; result.height = 38; }
    if (url.pathname.endsWith('cwindow.bed')) { result.name = 'Focal cWindow'; result.color = '#438ac1'; result.height = 38; }
    if (url.pathname.endsWith('features.bed')) { result.name = 'Peak / segment positions'; result.height = 40; }
    if (url.pathname.endsWith('focal_segments.bed')) { result.name = '5-kb segments · QTL Z and P'; result.height = 78; }
    if (url.pathname.endsWith('nearby_windows.bed')) { result.name = 'Same SNP · cWindow PCO P'; result.height = 48; }
    if (url.pathname.endsWith('available_regions.bed')) { result.name = 'Exported coverage regions'; result.height = 32; result.color = '#88a39a'; }
  }
  return result;
}

async function renderBrowser() {
  const view = el('view').value;
  const session = await readJSON(trackURL(`session_${view}.json`));
  const group = el('genotype').value;
  let tracks = session.tracks.map(localTrack);
  if (view === 'individuals' && group !== 'all') tracks = tracks.filter(t => t.type !== 'wig' || t.genotype === group);
  const config = {
    // igv.js supports chromosome-size-only references: no invented bases and
    // no external genome requests. Dimensions are the validated hg38 sizes.
    reference: {id: 'hg38', name: 'GRCh38 / hg38', format: 'chromsizes', url: new URL('hg38.chrom.sizes', pageBase).href, wholeGenomeView: false},
    loadDefaultGenomes: false, showIdeogram: false, showChromosomeWidget: false,
    showSVGButton: true, showCursorTrackGuide: true, showCenterGuide: false,
    locus: selectedLocus(), tracks,
    // Numeric chromosome coordinates only; never send searches to an API.
    search: {url: new URL('coordinate-search-disabled.json', pageBase).href},
  };
  if (state.browser) {
    // Reuse the viewer when switching sessions so its DOM observers stay valid.
    await state.browser.loadSessionObject(config);
  } else {
    state.browser = await igv.createBrowser(el('igv-browser'), config);
    state.browser.on('locuschange', () => {
      if (!state.busy) status('Ready');
    });
  }
  // Session loading starts an asynchronous coverage update. Let it finish
  // before enabling another session switch or recomputing shared mean scales.
  const deadline = Date.now() + 30000;
  while (state.browser.trackViews.some(v => v.viewports.some(p => p.isLoading()))) {
    if (Date.now() > deadline) throw new Error('Track loading timed out. Please reset the view.');
    await new Promise(resolve => setTimeout(resolve, 25));
  }
  await state.browser.updateViews();
  el('scale-note').textContent = view === 'means' ? 'Shared y-axis across genotype means' : `Fixed y-axis: 0–${state.summary.individual_ymax.toLocaleString('en-US')}`;
  el('browser-scroll').scrollTop = 0;
  status(`Ready · ${tracks.filter(t => t.type === 'wig').length} coverage tracks loaded`);
}

async function update(action) {
  if (state.busy) return;
  setBusy(true);
  status('Loading tracks…');
  try { await action(); }
  catch (error) { console.error(error); status(error.message || String(error), true); }
  finally { setBusy(false); }
}

document.querySelectorAll('[data-example]').forEach(button => {
  button.addEventListener('click', () => update(async () => {
    state.locus = button.dataset.example;
    el('region').value = 'cwindow';
    el('genotype').value = 'all';
    state.summary = await readJSON(trackURL('export_summary.json'));
    updateDetails();
    await renderBrowser();
  }));
});
el('view').addEventListener('change', () => update(renderBrowser));
el('genotype').addEventListener('change', () => update(renderBrowser));
el('region').addEventListener('change', () => update(async () => {
  if (!state.browser) return renderBrowser();
  const locus = selectedLocus();
  await state.browser.search(Array.isArray(locus) ? locus.join(' ') : locus);
  status('Ready');
}));
el('reset').addEventListener('click', () => update(async () => {
  el('view').value = 'means'; el('region').value = 'cwindow'; el('genotype').value = 'all';
  await renderBrowser();
}));

update(async () => {
  state.summary = await readJSON(trackURL('export_summary.json'));
  updateDetails();
  await renderBrowser();
});
