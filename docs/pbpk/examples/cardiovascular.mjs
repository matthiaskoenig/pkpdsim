import {column, csv, defaultControls, lastBeat, sampleAt, solvePressure} from "./cardiovascular-model.mjs";
import {createSchematic} from "./cardiovascular-schematic.mjs";
import {createHeartSolver, normalizeHeartControls} from "./cardiovascular-solver.mjs";

const $ = id => document.getElementById(id);
const state = {mode: "pressure", window: "response", caseId: "baseline", cache: new Map(),
  serial: 0, pending: false, flowTime: 0, playing: false, raf: 0, cursors: {}};
const palette = name => getComputedStyle(document.documentElement).getPropertyValue(`--${name}`).trim();
const text = (id, value) => { $(id).textContent = value; };
const format = (value, digits = 2) => value.toLocaleString("en", {minimumFractionDigits: digits, maximumFractionDigits: digits});
const signed = (value, digits = 2) => {
  const rounded = Math.abs(value) < 0.5 * 10 ** -digits ? 0 : value;
  return `${rounded >= 0 ? "+" : "−"}${format(Math.abs(rounded), digits)}`;
};
const svgNS = "http://www.w3.org/2000/svg";

function element(name, attributes = {}, content) {
  const node = document.createElementNS(svgNS, name);
  for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, value);
  if (content !== undefined) node.textContent = content;
  return node;
}

function plot(id, curves, {xLabel, yLabel, zero = false, intervention = false, markers = [], cursor = false, point = false}) {
  const host = $(id);
  const width = Math.max(280, host.clientWidth);
  const height = 270;
  const margin = {left: 53, right: 16, top: 28, bottom: 46};
  const points = curves.flatMap(curve => curve.x.map((x, index) => [x, curve.y[index]]));
  let xmin = Math.min(...points.map(point => point[0]));
  let xmax = Math.max(...points.map(point => point[0]));
  let ymin = Math.min(...points.map(point => point[1]));
  let ymax = Math.max(...points.map(point => point[1]));
  if (zero) ymin = Math.min(0, ymin);
  if (xmax === xmin) { xmin -= 0.5; xmax += 0.5; }
  const padding = Math.max((ymax - ymin) * 0.1, Math.abs(ymax) * 0.015, 0.05);
  ymin = zero && ymin === 0 ? 0 : ymin - padding;
  ymax += padding;
  const sx = x => margin.left + (x - xmin) / (xmax - xmin) * (width - margin.left - margin.right);
  const sy = y => height - margin.bottom - (y - ymin) / (ymax - ymin) * (height - margin.top - margin.bottom);
  const svg = element("svg", {viewBox: `0 0 ${width} ${height}`, "aria-hidden": "true"});
  svg.append(element("text", {x: margin.left, y: 14, class: "axis-label"}, yLabel));
  const decimals = ymax - ymin < 3 ? 1 : 0;
  for (let i = 0; i <= 4; i++) {
    const y = ymin + (ymax - ymin) * i / 4;
    svg.append(element("line", {x1: margin.left, x2: width - margin.right, y1: sy(y), y2: sy(y), class: "grid"}));
    svg.append(element("text", {x: margin.left - 9, y: sy(y) + 4, "text-anchor": "end"}, format(y, decimals)));
    const x = xmin + (xmax - xmin) * i / 4;
    svg.append(element("text", {x: sx(x), y: height - margin.bottom + 22, "text-anchor": "middle"}, format(x, xmax - xmin < 5 ? 1 : 0)));
  }
  svg.append(element("text", {x: margin.left + (width - margin.left - margin.right) / 2, y: height - 5, "text-anchor": "middle", class: "axis-label"}, xLabel));
  if (intervention && xmin <= 24 && xmax >= 24) {
    svg.append(element("line", {x1: sx(24), x2: sx(24), y1: margin.top, y2: height - margin.bottom, stroke: palette("muted"), "stroke-dasharray": "3 5", opacity: "0.65"}));
    svg.append(element("text", {x: sx(24) + 5, y: margin.top + 12}, "Intervention"));
  }
  for (const curve of curves) {
    const d = curve.x.map((x, index) => `${index ? "L" : "M"}${sx(x).toFixed(2)},${sy(curve.y[index]).toFixed(2)}`).join(" ");
    svg.append(element("path", {d, stroke: palette(curve.color), class: "curve", ...(curve.dashed ? {"stroke-dasharray": "5 4"} : {})}));
  }
  for (const marker of markers) {
    svg.append(element("circle", {cx: sx(marker.x), cy: sy(marker.y), r: 5, fill: palette(marker.color), stroke: palette("paper"), "stroke-width": 2}));
  }
  const indicator = point
    ? element("circle", {r: 5, fill: palette("accent"), stroke: palette("paper"), "stroke-width": 2})
    : element("line", {y1: margin.top, y2: height - margin.bottom, stroke: palette("accent"), "stroke-width": 1.5, "stroke-dasharray": "3 3"});
  if (cursor || point) svg.append(indicator);
  host.replaceChildren(svg);
  return (x, y) => {
    if (point) { indicator.setAttribute("cx", sx(x)); indicator.setAttribute("cy", sy(y)); }
    else { indicator.setAttribute("x1", sx(x)); indicator.setAttribute("x2", sx(x)); }
  };
}

function metric(index, label, value, unit, baseline, digits = 2, context = "vs reference") {
  text(`metric-label-${index}`, label);
  text(`metric-${index}`, format(value, digits));
  text(`unit-${index}`, unit);
  text(`delta-${index}`, `${signed(value - baseline, digits)} ${unit} ${context}`);
}

const CONTROL_IDS = ["hepatic", "systemic", "portal", "renal", "arterial", "venous"];
const HEART_PARAMETER = {hepatic: "CV_fR_hepatic", systemic: "CV_fR_systemic", portal: "CV_fR_portal", renal: "CV_fR_ki", heartRate: "CV_HR", contractility: "CV_contractility"};

function readControls() {
  return {...Object.fromEntries(CONTROL_IDS.map(id => [id, Number($(id).value)])),
    heartRate: Number($("heart-rate").value), contractility: Number($("contractility").value)};
}

function heartKey(controls) {
  return JSON.stringify(normalizeHeartControls(controls));
}

function exampleControls(example) {
  const controls = {...defaultControls(state.data.pressure), heartRate: 75, contractility: 1};
  for (const [key, parameter] of Object.entries(HEART_PARAMETER)) {
    if (example.changes[parameter] !== undefined) controls[key] = example.changes[parameter];
  }
  return controls;
}

function applyControls(controls) {
  for (const id of CONTROL_IDS) if (controls[id] !== undefined) $(id).value = controls[id];
  if (controls.heartRate !== undefined) $("heart-rate").value = controls.heartRate;
  if (controls.contractility !== undefined) $("contractility").value = controls.contractility;
  controlsChanged(0);
}

function controlsChanged(delay = 160) {
  if (!state.data) return;
  const controls = readControls();
  const key = heartKey(controls);
  state.caseId = state.data.cases.find(example => heartKey(exampleControls(example)) === key)?.id ?? "custom";
  $("heart-case").value = state.caseId;
  if (state.mode === "heart") requestHeart(delay);
  else render();
}

function setSchematic() {
  if (state.schematicMode !== state.mode) {
    state.schematic = createSchematic($("circulation-schematic"), state.mode === "heart");
    state.schematicMode = state.mode;
  }
}

function organTable(flows, baseline, dynamic = false) {
  $("organ-flows").replaceChildren(...flows.map(organ => {
    const row = document.createElement("tr");
    const reference = baseline.find(item => item.id === organ.id).flow;
    const drainage = organ.id === "lu" ? dynamic ? "Pulmonary vein" : "Arterial pool"
      : organ.portal ? "Portal vein" : organ.id === "li" ? "Sinusoids" : "Systemic vein";
    for (const value of [organ.name, drainage, format(organ.flow, 3), `${signed(100 * (organ.flow / reference - 1), 1)}%`]) {
      const cell = document.createElement("td"); cell.textContent = value; row.append(cell);
    }
    return row;
  }));
}

function pressureView() {
  const reference = state.data.pressure;
  const controls = readControls();
  const result = solvePressure(reference, controls);
  const baseline = solvePressure(reference, defaultControls(reference));
  state.pressureResult = result;
  metric(0, "Portal pressure", result.ppo, "mmHg", baseline.ppo, 1);
  metric(1, "Portal flow", result.qpo, "l/min", baseline.qpo);
  metric(2, "Hepatic flow", result.qh, "l/min", baseline.qh);
  metric(3, "Cardiac output", result.co, "l/min", baseline.co);
  state.schematic.update({CV_P_ar: controls.arterial, CV_P_ve: controls.venous,
    CV_P_po: result.ppo, CV_P_li: result.pli, Qpo: result.qpo, Qh: result.qh,
    Qha: result.qha, ...Object.fromEntries(result.flows.map(organ => [`Q${organ.id}`, organ.flow]))});
  text("pressure-insight", `Portal pressure is ${format(result.ppo, 1)} mmHg (${signed(result.ppo - baseline.ppo, 1)} from reference). ` +
    `The portal-to-systemic venous gradient is ${format(result.ppo - controls.venous, 1)} mmHg. ` +
    `The hepatic artery contributes ${format(100 * result.qha / result.qh, 0)}% of hepatic flow.`);
  const factors = Array.from({length: 71}, (_, i) => 1 + i / 10);
  const scan = factors.map(hepatic => solvePressure(reference, {...controls, hepatic}));
  plot("resistance-chart", [
    {x: factors, y: scan.map(row => row.ppo), color: "portal"},
    {x: factors, y: scan.map(row => row.pli), color: "hepatic"},
  ], {xLabel: "Hepatic outflow resistance (× reference)", yLabel: "Pressure (mmHg)", markers: [
    {x: controls.hepatic, y: result.ppo, color: "portal"},
    {x: controls.hepatic, y: result.pli, color: "hepatic"},
  ]});
  organTable([...result.flows, {id: "lu", name: "Lung", flow: result.co}],
    [...baseline.flows, {id: "lu", flow: baseline.co}]);
  text("table-method", "The liver row is hepatic arterial flow. Total hepatic flow also includes portal supply. Lung flow equals total systemic flow in this steady model.");
  text("calculation-method", "Live algebraic solution of the SBML resistance network. The same four resistance settings carry into the dynamic heart model.");
}

function remember(key, example) {
  state.cache.delete(key);
  state.cache.set(key, example);
  if (state.cache.size > 20) state.cache.delete(state.cache.keys().next().value);
}

function adoptHeart(example) {
  state.heartResult = example;
  const time = column(example.beats, "time");
  const duration = time.at(-1) - time[0];
  state.flowTime = Math.min(state.flowTime, duration);
  $("flow-time").max = String(duration);
  state.pending = false;
  state.solveError = null;
}

function requestHeart(delay = 160) {
  clearTimeout(state.solveTimer);
  state.abort?.abort();
  const serial = ++state.serial;
  const controls = normalizeHeartControls(readControls());
  const key = heartKey(controls);
  state.solveError = null;
  if (state.cache.has(key)) {
    adoptHeart(state.cache.get(key));
    render();
    return;
  }
  state.pending = true;
  render();
  state.solveTimer = setTimeout(async () => {
    const controller = new AbortController();
    state.abort = controller;
    try {
      state.solverPromise ??= createHeartSolver();
      const solver = await state.solverPromise;
      state.solver = solver;
      if (serial !== state.serial) return;
      const started = performance.now();
      const example = await solver.simulate(controls, {signal: controller.signal});
      if (serial !== state.serial) return;
      example.elapsed_ms = performance.now() - started;
      remember(key, example);
      adoptHeart(example);
    } catch (error) {
      if (serial !== state.serial || error.name === "AbortError") return;
      state.pending = false;
      state.solveError = error.message;
      if (!state.solver) state.solverPromise = null;
      console.error("Cardiovascular IDA simulation:", error);
    }
    if (serial === state.serial) render();
  }, delay);
}

function heartView() {
  const example = state.heartResult;
  const baseline = state.data.cases.find(item => item.id === "baseline");
  const preset = state.data.cases.find(item => item.id === state.caseId);
  text("case-description", preset ? preset.description : "Combined slider settings, applied after the shared 24-second reference settling interval.");
  $("solver-retry").hidden = !state.solveError;
  text("solver-status", state.solveError ? `Could not solve these settings: ${state.solveError}. Previous result remains visible.`
    : state.pending ? "Solving the current sliders with SUNDIALS IDA in a Web Worker. Previous result remains visible."
    : example?.elapsed_ms ? `Current settings solved in the browser with IDA in ${format(example.elapsed_ms / 1000, 2)} s.`
    : "Ready. Reference IDA solution for the current settings; changing any slider runs the compiled SBML in your browser.");
  if (!example) return;
  const summary = example.summary.after;
  const before = example.summary.before;
  metric(0, "Mean portal pressure", summary.mean_portal_pressure, "mmHg", before.mean_portal_pressure, 1, "vs before");
  metric(1, "Mean cardiac output", summary.mean_cardiac_output, "l/min", before.mean_cardiac_output, 2, "vs before");
  metric(2, "Stroke volume", summary.stroke_volume, "ml", before.stroke_volume, 1, "vs before");
  metric(3, "Heart rate", summary.heart_rate, "/min", before.heart_rate, 0, "vs before");
  text("heart-insight", `After settling: ${format(summary.arterial_systolic, 0)}/${format(summary.arterial_diastolic, 0)} mmHg systemic pressure. ` +
    `Mean hepatic flow ${format(summary.mean_hepatic_flow)} l/min. Intervention at 24 s.`);
  const beats = state.window === "beats";
  const selectedSeries = beats ? example.beats : example.series;
  const baselineSeries = beats ? baseline.beats : baseline.series;
  function time(series) {
    const values = column(series, "time");
    return beats ? values.map(value => value - values[0]) : values;
  }
  state.cursors = {};
  for (const [id, variable, label] of [
    ["arterial-chart", "CV_P_ar", "Pressure (mmHg)"],
    ["portal-chart", "CV_P_po", "Pressure (mmHg)"],
    ["flow-chart", "CO", "Flow (l/min)"],
  ]) {
    state.cursors[id] = plot(id, [
      {x: time(baselineSeries), y: column(baselineSeries, variable), color: "baseline", dashed: true},
      {x: time(selectedSeries), y: column(selectedSeries, variable), color: "accent"},
    ], {xLabel: beats ? "Time within final 3 beats (s)" : "Time (s)", yLabel: label, intervention: !beats, zero: variable === "CO", cursor: true});
  }
  const referenceBeat = lastBeat(baseline.beats, baseline.heart_rate);
  const selectedBeat = lastBeat(example.beats, example.heart_rate);
  state.cursors.pv = plot("pv-chart", [
    {x: column(referenceBeat, "CV_V_lv"), y: column(referenceBeat, "CV_P_lv"), color: "baseline", dashed: true},
    {x: column(selectedBeat, "CV_V_lv"), y: column(selectedBeat, "CV_P_lv"), color: "accent"},
  ], {xLabel: "Left ventricular volume (ml)", yLabel: "LV pressure (mmHg)", zero: true, point: true});
  const flows = summary.mean_flows;
  const referenceFlows = before.mean_flows;
  const allOrgans = [...state.data.pressure.organs, {id: "lu", name: "Lung"}];
  organTable(allOrgans.map(organ => ({...organ, flow: flows[organ.id === "li" ? "Qha" : `Q${organ.id}`]})),
    allOrgans.map(organ => ({...organ, flow: referenceFlows[organ.id === "li" ? "Qha" : `Q${organ.id}`]})), true);
  text("table-method", "Mean flows over the final three beats, compared with the same circulation before the intervention. The liver row is hepatic arterial flow; portal flow is additional.");
  text("conservation", `Maximum total blood-volume error: ${example.summary.maximum_blood_volume_error.toExponential(1)} l. ` +
    "Readouts use complete beats on the high-resolution IDA solution. Replay follows the final three beats; instantaneous flow sets every edge width.");
  text("calculation-method", "The compiled SBML and sbmlsim-core WASM runtime use SUNDIALS IDA in a Web Worker. Sliders may be combined freely. All runs start from the same reference state and change settings at 24 s.");
  displayFrame();
}

function displayFrame() {
  if (state.mode !== "heart" || !state.heartResult) return;
  const series = state.heartResult.beats;
  const start = series.values[0][series.columns.indexOf("time")];
  const row = sampleAt(series, start + state.flowTime);
  state.schematic.update(row);
  $("flow-time").value = String(state.flowTime);
  text("flow-clock", `${format(state.flowTime, 2)} s`);
  const x = state.window === "beats" ? state.flowTime : start + state.flowTime;
  for (const name of ["arterial-chart", "portal-chart", "flow-chart"]) state.cursors[name]?.(x);
  state.cursors.pv?.(row.CV_V_lv, row.CV_P_lv);
  $("play-flow").setAttribute("aria-pressed", String(state.playing));
  text("play-flow", state.playing ? "Pause" : "Play");
}

function stopAnimation() {
  if (state.raf) cancelAnimationFrame(state.raf);
  state.raf = 0;
  state.lastFrame = null;
}

function startAnimation() {
  if (state.mode !== "heart" || !state.playing || !state.heartResult || state.raf) return;
  const animate = now => {
    state.raf = 0;
    if (state.mode !== "heart" || !state.playing) return;
    const series = state.heartResult.beats;
    const times = column(series, "time");
    const duration = times.at(-1) - times[0];
    if (!state.pending && state.lastFrame !== null) {
      state.flowTime = (state.flowTime + Math.min((now - state.lastFrame) / 1000, 0.2) * Number($("flow-speed").value)) % duration;
      displayFrame();
    }
    state.lastFrame = now;
    state.raf = requestAnimationFrame(animate);
  };
  state.lastFrame = null;
  state.raf = requestAnimationFrame(animate);
}

function render() {
  if (!state.data) return;
  const pressure = state.mode === "pressure";
  $("pressure-results").hidden = !pressure;
  $("pressure-boundaries").hidden = !pressure;
  $("heart-controls").hidden = pressure;
  $("heart-results").hidden = pressure;
  $("flow-playback").hidden = pressure;
  $("download").disabled = !pressure && (state.pending || !!state.solveError);
  for (const button of document.querySelectorAll("[data-mode]")) button.setAttribute("aria-pressed", button.dataset.mode === state.mode);
  for (const button of document.querySelectorAll("[data-window]")) button.setAttribute("aria-pressed", button.dataset.window === state.window);
  const controls = readControls();
  for (const id of CONTROL_IDS) {
    const dimensional = id === "arterial" || id === "venous";
    text(`${id}-value`, `${format(controls[id], dimensional ? 0 : id === "systemic" ? 2 : 1)}${dimensional ? " mmHg" : "×"}`);
    $(id).setAttribute("aria-valuetext", `${controls[id]} ${dimensional ? "millimeters of mercury" : "times reference resistance"}`);
  }
  text("heart-rate-value", `${format(controls.heartRate, 0)}/min`);
  text("contractility-value", `${format(controls.contractility)}×`);
  text("control-method", pressure ? "Live algebraic solution" : "Live SUNDIALS IDA simulation");
  text("control-explanation", pressure ? "Prescribed pressure boundaries drive the steady network. Change any resistance and inspect the regional flows."
    : "The same resistance sliders now control the beating heart. Settings are applied at 24 s, then the closed-loop model settles; pressures are outputs.");
  text("schematic-description", pressure ? "The arterial and venous boundaries supply parallel tissue branches and the liver's dual blood supply. All edges share one flow-width scale."
    : "The systemic and portal beds join a closed pulmonary loop with two ventricles and four valves. Play or scrub the final three beats to inspect instantaneous flows and pressures.");
  setSchematic();
  if (pressure) { stopAnimation(); pressureView(); }
  else { heartView(); startAnimation(); }
  resizeFrame();
}

let previousHeight = 0;
function resizeFrame() {
  if (window.parent === window) return;
  const height = Math.ceil($("app").getBoundingClientRect().height);
  if (Math.abs(height - previousHeight) < 2) return;
  previousHeight = height;
  window.parent.postMessage({pkpdFigure: true, height}, window.location.origin);
}

function theme() {
  let scheme;
  try { scheme = window.parent.document.body.getAttribute("data-md-color-scheme"); } catch { /* Standalone or cross-origin embedding. */ }
  document.documentElement.dataset.theme = scheme === "slate" || (!scheme && matchMedia("(prefers-color-scheme: dark)").matches) ? "dark" : "light";
  render();
}

async function load() {
  $("load-error").hidden = true;
  $("loading").hidden = false;
  try {
    const response = await fetch("cardiovascular-data.json");
    if (!response.ok) throw new Error(`Reference data request failed: ${response.status}`);
    const data = await response.json();
    if (data.schema_version !== 1 || !data.cases?.[0]?.beats.columns.includes("CV_Q_tv")) throw new Error("The circulation data needs regeneration with the current model.");
    state.data = data;
    state.cache.clear();
    for (const example of data.cases) remember(heartKey(exampleControls(example)), example);
    adoptHeart(data.cases.find(example => example.id === "baseline"));
    const custom = document.createElement("option"); custom.value = "custom"; custom.textContent = "Custom slider settings"; custom.disabled = true;
    $("heart-case").replaceChildren(custom, ...data.cases.map(example => {
      const option = document.createElement("option"); option.value = example.id; option.textContent = example.label; return option;
    }));
    $("workspace").hidden = false;
    $("loading").hidden = true;
    $("app").setAttribute("aria-busy", "false");
    applyControls({...defaultControls(data.pressure), heartRate: 75, contractility: 1});
  } catch (error) {
    $("loading").hidden = true;
    $("load-error").hidden = false;
    $("app").setAttribute("aria-busy", "false");
    console.error("Cardiovascular example:", error);
    resizeFrame();
  }
}

for (const button of document.querySelectorAll("[data-mode]")) button.addEventListener("click", () => {
  state.mode = button.dataset.mode;
  if (state.mode === "heart") {
    if (!state.hasPlayed) {
      state.playing = !matchMedia("(prefers-reduced-motion: reduce)").matches;
      state.hasPlayed = true;
    }
    requestHeart(0);
  } else {
    clearTimeout(state.solveTimer); state.abort?.abort(); ++state.serial; state.pending = false;
    render();
  }
});
for (const button of document.querySelectorAll("[data-window]")) button.addEventListener("click", () => { state.window = button.dataset.window; render(); });
for (const id of [...CONTROL_IDS, "heart-rate", "contractility"]) $(id).addEventListener("input", () => controlsChanged());
for (const button of document.querySelectorAll("[data-preset]")) button.addEventListener("click", () => {
  const controls = {...defaultControls(state.data.pressure), heartRate: 75, contractility: 1};
  if (button.dataset.preset === "hepatic") controls.hepatic = 4;
  if (button.dataset.preset === "vasodilator") controls.systemic = 0.7;
  applyControls(controls);
});
$("reset").addEventListener("click", () => applyControls({...defaultControls(state.data.pressure), heartRate: 75, contractility: 1}));
$("heart-case").addEventListener("change", event => {
  const example = state.data.cases.find(item => item.id === event.target.value);
  if (example) applyControls(exampleControls(example));
});
$("retry").addEventListener("click", load);
$("solver-retry").addEventListener("click", () => { state.solver?.close(); state.solver = null; state.solverPromise = null; requestHeart(0); });
$("play-flow").addEventListener("click", () => {
  state.playing = !state.playing;
  if (state.playing) startAnimation(); else stopAnimation();
  displayFrame();
});
$("flow-time").addEventListener("input", () => {
  state.playing = false; stopAnimation(); state.flowTime = Number($("flow-time").value); displayFrame();
});
$("flow-speed").addEventListener("change", () => { state.lastFrame = null; });
$("download").addEventListener("click", () => {
  let series, units, filename;
  if (state.mode === "pressure") {
    const result = state.pressureResult;
    const organs = [...result.flows, {id: "lu", flow: result.co}];
    series = {columns: ["organ", "flow", "portal_pressure", "sinusoidal_pressure"], values: organs.map(organ => [organ.id, organ.flow, result.ppo, result.pli])};
    units = {flow: "l/min", portal_pressure: "mmHg", sinusoidal_pressure: "mmHg"};
    filename = "pressure-network.csv";
  } else {
    const example = state.heartResult;
    series = state.window === "beats" ? example.beats : example.series;
    units = state.data.provenance.units;
    filename = `cardiovascular-${state.caseId}-${state.window}.csv`;
  }
  const url = URL.createObjectURL(new Blob([csv(series, units)], {type: "text/csv;charset=utf-8"}));
  const link = document.createElement("a"); link.href = url; link.download = filename; link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});

let resizeTimer;
new ResizeObserver(() => { clearTimeout(resizeTimer); resizeTimer = setTimeout(render, 80); }).observe($("app"));
try {
  if (window.parent !== window) new MutationObserver(theme).observe(window.parent.document.body, {attributes: true, attributeFilter: ["data-md-color-scheme"]});
} catch { /* Cross-origin embedding uses the system theme. */ }
window.addEventListener("pagehide", () => { clearTimeout(state.solveTimer); state.abort?.abort(); state.solver?.close(); stopAnimation(); });
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", theme);
theme();
load();
