// The liver lobe explorer: the clickable liver, the controls, the plots and
// the download. The numbers are either the precomputed courses of
// `lobes-data.json`, which scripts/lobes_examples.py generated with SUNDIALS
// IDA, or a live run of the same compiled SBML in a Web Worker; which one is
// shown is stated on the page. No equation of a lobed liver lives in this
// file.

import {STATES, caseChanges, csv, exposure, sampleAt, stateOfCase, trim} from "./lobes-model.mjs";
import {createLiver} from "./lobes-schematic.mjs";
import {CONTROLS, createLobeSolver, normalizeControls} from "./lobes-solver.mjs";

const DAY = 24 * 60;
const SPECIES = ["rat", "mouse", "human"];
const SVG = "http://www.w3.org/2000/svg";
const CURVES = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#9c755f", "#76b7b2", "#af7aa1", "#edc948"];

function element(tag, attributes = {}, text) {
  const node = document.createElement(tag);
  for (const [name, value] of Object.entries(attributes)) {
    if (name === "class") node.className = value;
    else node.setAttribute(name, String(value));
  }
  if (text !== undefined) node.textContent = text;
  return node;
}

function svgNode(tag, attributes = {}, text) {
  const node = document.createElementNS(SVG, tag);
  for (const [name, value] of Object.entries(attributes)) node.setAttribute(name, String(value));
  if (text !== undefined) node.textContent = text;
  return node;
}

/** A line plot of several courses over the days, axes and legend included. */
export function linePlot(host, {title, unit, times, curves, marker}) {
  const width = 460, height = 232, left = 54, bottom = 30, top = 36, right = 10;
  const all = curves.flatMap(curve => curve.values).filter(Number.isFinite);
  const max = Math.max(...all, 0) || 1;
  const min = Math.min(...all, 0);
  const end = times.at(-1) ?? 1;
  const x = time => left + ((width - left - right) * time) / (end || 1);
  const y = value => height - bottom - ((height - bottom - top) * (value - min)) / (max - min || 1);
  const svg = svgNode("svg", {class: "core-diagram plot", viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": title});
  svg.append(svgNode("title", {}, title));
  svg.append(svgNode("text", {x: left, y: 12, class: "cd-title"}, `${title} [${unit}]`));
  for (let day = 0; day <= Math.round(end / DAY); day++) {
    svg.append(svgNode("line", {class: "plot-grid", x1: x(day * DAY), x2: x(day * DAY), y1: top, y2: height - bottom}));
    svg.append(svgNode("text", {class: "cd-sub", x: x(day * DAY), y: height - bottom + 14, "text-anchor": "middle"}, String(day)));
  }
  svg.append(svgNode("text", {class: "cd-sub", x: width - right, y: height - 4, "text-anchor": "end"}, "days"));
  for (const bound of [min, min + (max - min) / 2, max]) {
    svg.append(svgNode("line", {class: "plot-grid", x1: left, x2: width - right, y1: y(bound), y2: y(bound)}));
    svg.append(svgNode("text", {class: "cd-sub", x: left - 6, y: y(bound) + 4, "text-anchor": "end"}, trim(bound)));
  }
  curves.forEach((curve, index) => {
    const points = curve.values.map((value, i) => `${x(times[i])},${y(value)}`).join(" ");
    svg.append(svgNode("polyline", {class: "plot-curve", points, style: `stroke:${CURVES[index % CURVES.length]}`}));
    const step = (width - left - right) / Math.max(curves.length, 1);
    svg.append(svgNode("text", {class: "cd-sub plot-legend", x: left + index * step, y: 26, style: `fill:${CURVES[index % CURVES.length]}`}, curve.label));
  });
  if (Number.isFinite(marker)) svg.append(svgNode("line", {class: "plot-marker", x1: x(marker), x2: x(marker), y1: top, y2: height - bottom}));
  host.replaceChildren(svg);
  return svg;
}

/** A bar plot of a value per lobe at the shown time. */
export function barPlot(host, {title, unit, bars}) {
  const width = 460, height = 206, left = 54, bottom = 24, top = 36, right = 10;
  const max = Math.max(...bars.flatMap(bar => bar.values), 0) || 1;
  const step = (width - left - right) / Math.max(bars.length, 1);
  const y = value => height - bottom - ((height - bottom - top) * value) / max;
  const svg = svgNode("svg", {class: "core-diagram plot", viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": title});
  svg.append(svgNode("title", {}, title));
  svg.append(svgNode("text", {x: left, y: 12, class: "cd-title"}, `${title} [${unit}]`));
  svg.append(svgNode("line", {class: "plot-axis", x1: left, x2: width - right, y1: height - bottom, y2: height - bottom}));
  for (const bound of [max / 2, max]) {
    svg.append(svgNode("line", {class: "plot-grid", x1: left, x2: width - right, y1: y(bound), y2: y(bound)}));
    svg.append(svgNode("text", {class: "cd-sub", x: left - 6, y: y(bound) + 4, "text-anchor": "end"}, trim(bound)));
  }
  bars.forEach((bar, index) => {
    const slot = left + index * step;
    // the bars of a lobe are a centred group, so one bar per lobe is as wide
    // as two are together and the slots stay the same whatever is plotted
    const parts = Math.max(bar.values.length, 1);
    const barWidth = Math.max(Math.min((step - 8) / parts - 2, 34), 1);
    const start = slot + (step - (barWidth * parts + 2 * (parts - 1))) / 2;
    bar.values.forEach((value, part) => {
      svg.append(svgNode("rect", {
        class: "plot-bar", x: start + part * (barWidth + 2), y: y(value), width: barWidth,
        height: Math.max(height - bottom - y(value), 0), style: `fill:${CURVES[part % CURVES.length]}`,
      }));
    });
    svg.append(svgNode("text", {class: "cd-sub", x: slot + step / 2, y: height - bottom + 14, "text-anchor": "middle"}, bar.label));
  });
  bar_legend(svg, bars[0]?.parts ?? [], left, 26);
  host.replaceChildren(svg);
  return svg;
}

function bar_legend(svg, parts, left, y) {
  parts.forEach((part, index) => {
    svg.append(svgNode("text", {class: "cd-sub", x: left + index * 74, y, style: `fill:${CURVES[index % CURVES.length]}`}, part));
  });
}

/** The explorer of the lobes of the liver on the page it is hosted in. */
export function createApp(root, {data, solverFactory = createLobeSolver} = {}) {
  const parts = {
    species: root.querySelector("#species"),
    liver: root.querySelector("#liver"),
    removed: root.querySelector("#removed"),
    legend: root.querySelector("#legend"),
    presets: root.querySelector("#presets"),
    sliders: root.querySelector("#sliders"),
    status: root.querySelector("#status"),
    reset: root.querySelector("#reset"),
    masses: root.querySelector("#plot-masses"),
    flows: root.querySelector("#plot-flows"),
    remnant: root.querySelector("#plot-remnant"),
    exposure: root.querySelector("#plot-exposure"),
    scrub: root.querySelector("#scrub"),
    clock: root.querySelector("#clock"),
    download: root.querySelector("#download"),
    provenance: root.querySelector("#provenance"),
  };
  const state = {
    species: "rat",
    states: {},
    controls: normalizeControls({}),
    time: 0,
    series: null,
    source: "precomputed",
    solver: null,
    running: false,
    again: false,
  };
  let liver = null;

  const of_species = () => data.species[state.species];
  const tokens = () => of_species().lobes.map(lobe => lobe.token);
  const reset = () => { state.states = Object.fromEntries(tokens().map(token => [token, "intact"])); };

  const matching = () => {
    const changes = JSON.stringify(caseChanges(state.states));
    const [name] = Object.entries(of_species().cases).find(([, entry]) => JSON.stringify(entry.changes) === changes) ?? [];
    return name;
  };

  const defaults = () => Object.entries(CONTROLS).every(([key, control]) => state.controls[key] === control.value);

  const show = () => {
    const sample = sampleAt(state.series, state.time);
    // a resected lobe has no mass in the liver any more, so it holds no area
    // of the diagram; it is listed beside it and stays one click from coming
    // back. What the page shows is the state the reader chose, which a run
    // which is still on its way or which failed cannot contradict.
    const present = token => state.states[token] !== "resected";
    const shares = Object.fromEntries(
      tokens().map(token => [token, present(token) ? (sample[`f_vol_${token}`] ?? 0) : 0]),
    );
    const gone = tokens().filter(token => !present(token));
    parts.removed.replaceChildren(...gone.map(token => {
      const chip = element("button", {class: "chip", type: "button", "data-lobe": token, title: "put the lobe back"}, `${token} resected`);
      chip.addEventListener("click", () => void pick(token));
      return chip;
    }));
    parts.removed.hidden = gone.length === 0;
    liver.update({
      shares,
      states: state.states,
      values: Object.fromEntries(tokens().map(token => [token, `${trim(1000 * (sample[`Vli_${token}`] ?? 0))} ml`])),
      flows: Object.fromEntries(tokens().flatMap(token => [[`ha_${token}`, sample[`Qha_${token}`]], [`po_${token}`, sample[`Qpo_${token}`]]])),
      flowLabel: `${trim(1000 * (sample.Qh ?? 0))} ml/min`,
    });
    const times = state.series.time ?? [];
    const shown = tokens().filter(present);
    linePlot(parts.masses, {
      title: "volume of the lobes", unit: "ml", times, marker: state.time,
      curves: shown.map(token => ({label: token, values: (state.series[`Vli_${token}`] ?? []).map(v => 1000 * v)})),
    });
    barPlot(parts.flows, {
      title: `blood flow per lobe at day ${trim(state.time / DAY)}`, unit: "ml/min",
      bars: shown.map(token => ({
        label: token, parts: ["hepatic artery", "portal vein"],
        values: [1000 * (sample[`Qha_${token}`] ?? 0), 1000 * (sample[`Qpo_${token}`] ?? 0)],
      })),
    });
    linePlot(parts.remnant, {
      title: "the liver which is there", unit: "fraction of the intact liver", times, marker: state.time,
      curves: [{label: "remnant", values: state.series.f_remnant ?? []}],
    });
    // the exposure of a lobe is the area under its own outflow over the days
    // which are shown, the readout a study of a regional intervention takes
    barPlot(parts.exposure, {
      title: "exposure of the lobes over the five days", unit: "mM*min",
      bars: shown.map(token => ({
        label: token, parts: ["area under the concentration leaving the lobe"],
        values: [exposure(times, state.series[`Cout_li_${token}_drg`] ?? [])],
      })),
    });
    parts.status.textContent = `${state.source === "live" ? "Live run of the compiled model" : "Precomputed reference run"}: ${describe()}.`;
    parts.clock.textContent = `day ${trim(state.time / DAY)}`;
    parts.provenance.textContent = `${data.provenance.solver} of sbmlsim-core ${data.provenance.sbmlsim_core_version}, model ${of_species().sbml_sha256.slice(0, 12)}. ${data.provenance.notice}`;
  };

  const describe = () => {
    const changed = Object.entries(state.states).filter(([, value]) => value !== "intact");
    if (changed.length === 0) return "the intact liver";
    return changed.map(([token, value]) => `${token} ${STATES.find(s => s.id === value).short}`).join(", ");
  };

  const render = async () => {
    const name = matching();
    if (name !== undefined && defaults()) {
      state.series = of_species().cases[name].series;
      state.source = "precomputed";
      show();
      return;
    }
    // a click while a run is on its way is not lost: it is run once the
    // running one is back, so the page ends on the state the reader chose
    if (state.running) {
      state.again = true;
      return;
    }
    state.running = true;
    parts.status.textContent = "Solving the lobed liver with SUNDIALS IDA...";
    try {
      state.solver ??= await solverFactory({species: SPECIES});
      const result = await state.solver.simulate(state.species, state.controls, caseChanges(state.states));
      state.series = result.series;
      state.source = "live";
      show();
    } catch (error) {
      parts.status.textContent = `The live run failed: ${error.message} The page keeps the last result.`;
    } finally {
      state.running = false;
    }
    if (state.again) {
      state.again = false;
      await render();
    }
  };

  const pick = token => {
    const order = STATES.map(entry => entry.id);
    const next = order[(order.indexOf(state.states[token]) + 1) % order.length];
    state.states = {...state.states, [token]: next};
    return render();
  };

  // The listeners of the elements of the page itself are attached once; the
  // children which `build` replaces carry their own and go with them.
  const wire = () => {
    parts.scrub.min = 0;
    parts.scrub.max = String(data.provenance.end_time);
    parts.scrub.step = String(data.provenance.end_time / (data.provenance.points - 1));
    parts.scrub.addEventListener("input", () => { state.time = Number(parts.scrub.value); show(); });
    parts.reset.addEventListener("click", () => {
      reset();
      state.controls = normalizeControls({});
      state.time = 0;
      parts.scrub.value = "0";
      build();
      void render();
    });
    parts.download.addEventListener("click", () => {
      const columns = ["time", "f_remnant", ...tokens().flatMap(token => [`Vli_${token}`, `Qha_${token}`, `Qpo_${token}`, `Cout_li_${token}_drg`])];
      const text = csv(state.series, columns.filter(column => column in state.series));
      const url = URL.createObjectURL(new Blob([text], {type: "text/csv"}));
      const link = element("a", {href: url, download: `lobes-${state.species}-${describe().replaceAll(/\W+/g, "-")}.csv`});
      link.click();
      URL.revokeObjectURL(url);
    });
  };

  const build = () => {
    parts.species.replaceChildren(...SPECIES.map(name => {
      const button = element("button", {class: name === state.species ? "tab active" : "tab", type: "button"}, name);
      button.addEventListener("click", () => void setSpecies(name));
      return button;
    }));
    parts.legend.replaceChildren(...STATES.map(entry => {
      const item = element("span", {class: "legend-item", "data-state": entry.id, title: entry.description});
      item.append(element("span", {class: "legend-swatch"}), element("span", {class: "legend-label"}, entry.label));
      return item;
    }));
    parts.sliders.replaceChildren(...Object.entries(CONTROLS).map(([key, control]) => {
      const row = element("label", {class: "slider"});
      const value = element("output", {class: "slider-value"}, String(state.controls[key]));
      const input = element("input", {type: "range", min: control.min, max: control.max, step: control.step, value: state.controls[key]});
      input.addEventListener("input", () => {
        state.controls = normalizeControls({...state.controls, [key]: Number(input.value)});
        value.textContent = String(state.controls[key]);
        void render();
      });
      row.append(element("span", {class: "slider-label"}, control.label), input, value);
      return row;
    }));
  };

  const setSpecies = name => {
    state.species = name;
    reset();
    state.time = 0;
    parts.scrub.value = "0";
    liver = createLiver(parts.liver, of_species().lobes, {onPick: pick});
    build();
    buildPresets();
    return render();
  };

  const buildPresets = () => {
    parts.presets.replaceChildren(...Object.entries(of_species().cases).map(([name, entry]) => {
      const button = element("button", {class: "preset", type: "button"}, entry.label);
      button.addEventListener("click", () => {
        state.states = stateOfCase(entry, tokens());
        void render();
      });
      return button;
    }));
  };

  reset();
  wire();
  build();
  buildPresets();
  liver = createLiver(parts.liver, of_species().lobes, {onPick: pick});
  return {
    state,
    pick,
    render,
    setSpecies,
    setTime(time) { state.time = time; show(); },
    close() { state.solver?.close(); },
  };
}

/** Start the explorer in the page which hosts it. */
export async function start(root = document) {
  const response = await fetch(new URL("lobes-data.json", import.meta.url));
  if (!response.ok) throw new Error(`Could not load the reference data: HTTP ${response.status}.`);
  const data = await response.json();
  const app = createApp(root, {data});
  await app.render();
  const theme = () => {
    const parent = globalThis.parent?.document?.body?.dataset?.mdColorScheme;
    const dark = parent === "slate" || (parent === undefined && globalThis.matchMedia?.("(prefers-color-scheme: dark)").matches);
    document.documentElement.dataset.theme = dark ? "dark" : "light";
  };
  theme();
  globalThis.matchMedia?.("(prefers-color-scheme: dark)").addEventListener?.("change", theme);
  // the page is embedded in the documentation as a frame, which follows the
  // height the page posts; the liver and the plots change it, so it is posted
  // whenever the layout settles and not only once
  let posted = 0;
  const height = () => {
    const now = document.documentElement.scrollHeight;
    if (Math.abs(now - posted) < 2) return;
    posted = now;
    globalThis.parent?.postMessage?.({pkpdFigure: true, height: now}, "*");
  };
  height();
  globalThis.addEventListener?.("resize", height);
  if (globalThis.ResizeObserver) new globalThis.ResizeObserver(() => height()).observe(document.body);
  return app;
}
