// The liver lobe explorer: the clickable liver, the controls, the plots and
// the download. The numbers are either the precomputed courses of
// `lobes-data.json`, which scripts/lobes_examples.py generated with SUNDIALS
// IDA, or a live run of the same compiled SBML in a Web Worker; which one is
// shown is stated on the page. No equation of a lobed liver lives in this
// file.

import {STATES, VIEWS, caseChanges, csv, duration, exposure, lobeColors, planSummary, sampleAt, stateOfCase, trim} from "./lobes-model.mjs";
import {createLiver} from "./lobes-schematic.mjs";
import {CONTROLS, SIMULATION, createLobeSolver, inflowChanges, normalizeControls, normalizeSimulation} from "./lobes-solver.mjs";

const DAY = 24 * 60;
const SPECIES = ["rat", "mouse", "human"];
const SEXES = ["male", "female"];
const SVG = "http://www.w3.org/2000/svg";
// the days of the course one second of the animation covers
const PLAY_DAYS_PER_SECOND = 0.8;
// the operations a reader plans: the intact liver and the three states a lobe
// is put in, in the order of the states, so the row of the planner, the shades
// of the legend and the groups of the named plans are read in one order
const PROCEDURES = STATES.map(state => ({
  id: state.id,
  label: state.label,
  short: state.short,
  // what the button of the operation says, and the name of the group of named
  // plans it heads
  button: {
    intact: "intact",
    portal: "PVL - portal vein ligation",
    arterial: "HAL - hepatic artery ligation",
    resected: "resection",
  }[state.id] ?? state.label,
  group: {intact: "intact", portal: "PVL", arterial: "HAL", resected: "resection"}[state.id] ?? state.short,
}));

/** What a named plan is called under the operation it is grouped with.
 *
 * The case names of the data say the operation and the lobes it acts on
 * (`pvl70`, `right_hepatectomy`); the group above them already says the
 * operation, so the button says what is left of the name: the share of the
 * mass, the side or the lobe. A name which is nothing but its operation keeps
 * it.
 */
export function presetLabel(label, kind) {
  const rest = label.replaceAll(/pvl|hal|hepatectomy|resection|intact/gi, " ").replaceAll(/[\s_]+/g, " ").trim();
  if (rest === "") return kind === "intact" ? "whole liver" : label.replaceAll("_", " ");
  return /^\d+$/.test(rest) ? `${rest} %` : rest;
}

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
  const width = 460, height = 196, left = 50, bottom = 28, top = 34, right = 10;
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
    const colour = curve.colour ?? "var(--lobe-other)";
    const points = curve.values.map((value, i) => `${x(times[i])},${y(value)}`).join(" ");
    svg.append(svgNode("polyline", {class: "plot-curve", points, style: `stroke:${colour}`}));
    const step = (width - left - right) / Math.max(curves.length, 1);
    svg.append(
      svgNode("rect", {class: "plot-swatch", x: left + index * step, y: 19, width: 8, height: 8, rx: 2, style: `fill:${colour}`}),
      svgNode("text", {class: "cd-sub plot-legend", x: left + index * step + 11, y: 26}, curve.label),
    );
  });
  if (Number.isFinite(marker)) svg.append(svgNode("line", {class: "plot-marker", x1: x(marker), x2: x(marker), y1: top, y2: height - bottom}));
  host.replaceChildren(svg);
  return svg;
}

/** A bar plot of a value per lobe at the shown time. */
export function barPlot(host, {title, unit, bars}) {
  const width = 460, height = 176, left = 50, bottom = 24, top = 34, right = 10;
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
      // the hue of a bar is its lobe, so one bar per lobe and a group of two
      // read the same way; the parts of a group are told apart by the shade
      svg.append(svgNode("rect", {
        class: `plot-bar plot-bar-${part}`, x: start + part * (barWidth + 2), y: y(value), width: barWidth,
        height: Math.max(height - bottom - y(value), 0), rx: 2,
        style: `fill:${bar.colour ?? "var(--lobe-other)"}`,
      }));
    });
    svg.append(svgNode("text", {class: "cd-sub", x: slot + step / 2, y: height - bottom + 14, "text-anchor": "middle"}, bar.label));
  });
  bar_legend(svg, bars[0]?.parts ?? [], left, 26);
  host.replaceChildren(svg);
  return svg;
}

function bar_legend(svg, parts, left, y) {
  // the parts of a group share the hue of their lobe, so the legend shows the
  // two shades of one neutral swatch and never a second set of hues
  let x = left;
  parts.forEach((part, index) => {
    svg.append(svgNode("rect", {class: `plot-swatch plot-bar-${index}`, x, y: y - 7, width: 8, height: 8, rx: 2}));
    const label = svgNode("text", {class: "cd-sub", x: x + 11, y}, part);
    svg.append(label);
    x += 11 + 6.2 * part.length + 14;
  });
}

/** The explorer of the lobes of the liver on the page it is hosted in. */
export function createApp(root, {data, anatomy = null, loadData = defaultLoadData, solverFactory = createLobeSolver} = {}) {
  const parts = {
    species: root.querySelector("#species"),
    sex: root.querySelector("#sex"),
    liver: root.querySelector("#liver"),
    removed: root.querySelector("#removed"),
    legend: root.querySelector("#legend"),
    presets: root.querySelector("#presets"),
    sliders: root.querySelector("#sliders"),
    status: root.querySelector("#status"),
    reset: root.querySelector("#reset"),
    kind: root.querySelector("#plan-kind"),
    lobes: root.querySelector("#plan-lobes"),
    summary: root.querySelector("#plan-summary"),
    sim: root.querySelector("#sim"),
    spinner: root.querySelector("#spinner"),
    runtime: root.querySelector("#runtime"),
    play: root.querySelector("#play"),
    views: root.querySelector("#views"),
    body: root.querySelector("#body"),
    masses: root.querySelector("#plot-masses"),
    flows: root.querySelector("#plot-flows"),
    remnant: root.querySelector("#plot-remnant"),
    exposure: root.querySelector("#plot-exposure"),
    scrub: root.querySelector("#scrub"),
    clock: root.querySelector("#clock"),
    download: root.querySelector("#download"),
    provenance: root.querySelector("#provenance"),
  };
  // the dataset of a sex is held once it was loaded: the page starts with the
  // one it was given and fetches the other when a reader asks for it
  const datasets = {[data.provenance.sex ?? "male"]: data};
  const solvers = {};
  const state = {
    species: "rat",
    sex: data.provenance.sex ?? "male",
    states: {},
    controls: normalizeControls({}),
    simulation: normalizeSimulation({}),
    time: 0,
    series: null,
    source: "precomputed",
    running: false,
    again: false,
    plan: {kind: "portal"},
    view: "volume",
    playing: false,
    ran: null,
  };
  let liver = null;
  let colours = {};
  let frame = null;

  const dataset = () => datasets[state.sex];
  const of_species = () => dataset().species[state.species];
  // the sexes of the species: the parameter database has the female of the
  // human, and the male alone of the rat and of the mouse
  const sexesOf = species => dataset().species[species]?.sexes ?? ["male"];
  const tokens = () => of_species().lobes.map(lobe => lobe.token);
  const reset = () => { state.states = Object.fromEntries(tokens().map(token => [token, "intact"])); };
  // the end of the course is the run which is drawn, not the one which was
  // precomputed: a reader who asks for ten days scrubs over ten days
  const endTime = () => state.series?.time?.at(-1) ?? dataset().provenance.end_time;

  // a set of changes as one string, so two sets are the same when they set the
  // same parameters to the same values whatever order they were written in
  const signature = changes => Object.entries(changes ?? {}).map(([key, value]) => `${key}=${value}`).sort().join(";");

  // what a run changes: the gates the inflow controls ask for, with the
  // ligations of the plan on top, so a ligated branch stays closed
  const changesOfRun = () => ({...inflowChanges(tokens(), state.states, state.controls), ...caseChanges(state.states)});

  const matching = () => {
    const changes = signature(caseChanges(state.states));
    const [name] = Object.entries(of_species().cases).find(([, entry]) => signature(entry.changes) === changes) ?? [];
    return name;
  };

  // a precomputed case is the run of the default physiology over the days and
  // the dose it was precomputed with: a reader who moves any of them is asking
  // for another run, which the compiled model answers here
  const defaults = () =>
    Object.entries(CONTROLS).every(([key, control]) => state.controls[key] === control.value)
    && Object.entries(SIMULATION).every(([key, setting]) => state.simulation[key] === setting.value);

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
    // a resected lobe holds no area of the liver, so its box sits beside the
    // drawing: the colour of the lobe, shaded as what was done to it
    parts.removed.replaceChildren(...gone.map(token => {
      const colour = colours[token] ?? "var(--lobe-other)";
      const chip = element("button", {class: "chip", type: "button", "data-lobe": token, title: "put the lobe back"});
      chip.append(
        element("span", {
          class: "lobe-dot",
          style: `background:repeating-linear-gradient(135deg, ${colour} 0 2px, transparent 2px 5px), color-mix(in srgb, ${colour} 22%, transparent)`,
        }),
        element("span", {}, `${token} resected`),
      );
      chip.addEventListener("click", () => void pick(token));
      return chip;
    }));
    parts.removed.hidden = gone.length === 0;
    // the area of the whole liver is the mass which is there, measured against
    // the intact liver of the same individual, so a hepatectomy shrinks the
    // drawing and the regrowth fills it again
    const present_volume = tokens().filter(present).reduce((sum, token) => sum + (sample[`Vli_${token}`] ?? 0), 0);
    const reference = sample.Vli_ref ?? intactVolume();
    liver.update({
      shares,
      scale: reference > 0 ? present_volume / reference : 1,
      states: state.states,
      colours,
      values: Object.fromEntries(tokens().map(token => [
        token,
        state.view === "pressure"
          ? `${trim(sample[`CV_P_li_${token}`] ?? 0)} mmHg`
          : `${trim(1000 * (sample[`Vli_${token}`] ?? 0))} ml`,
      ])),
      flows: Object.fromEntries(tokens().flatMap(token => [
        [`ha_${token}`, sample[`Qha_${token}`]],
        [`po_${token}`, sample[`Qpo_${token}`]],
        [`hv_${token}`, sample[`Qh_${token}`]],
      ])),
      flowLabel: `${trim(1000 * (sample.Qh ?? 0))} ml/min`,
    });
    const times = state.series.time ?? [];
    const shown = tokens().filter(present);
    const day = trim(state.time / DAY);
    if (state.view === "pressure") {
      // the same run read as the pressure network: what drives the flows
      linePlot(parts.masses, {
        title: "pressure of the hepatic inflow", unit: "mmHg", times, marker: state.time,
        curves: [
          {label: "portal vein", colour: "var(--edge-po)", values: state.series.CV_P_po ?? []},
          {label: "sinusoids", colour: "var(--edge-hv)", values: state.series.CV_P_li ?? []},
        ],
      });
      // the pressure of every lobe, which is where a ligation is read off:
      // a deportalised lobe sits at the pressure of its artery alone
      barPlot(parts.flows, {
        title: `pressure per lobe at day ${day}`, unit: "mmHg",
        bars: shown.map(token => ({
          label: token, parts: ["pressure of the lobe"], colour: colours[token],
          values: [sample[`CV_P_li_${token}`] ?? 0],
        })),
      });
      linePlot(parts.remnant, {
        title: "the hepatic inflow of the liver", unit: "ml/min", times, marker: state.time,
        curves: [{label: "hepatic inflow", colour: "var(--accent)", values: (state.series.Qh ?? []).map(v => 1000 * v)}],
      });
      barPlot(parts.exposure, {
        title: `perfusion per volume at day ${day}`, unit: "1/min",
        bars: shown.map(token => ({
          label: token, parts: ["flow per volume of the lobe"], colour: colours[token],
          values: [sample[`q_${token}`] ?? 0],
        })),
      });
    } else {
      linePlot(parts.masses, {
        title: "volume of the lobes", unit: "ml", times, marker: state.time,
        curves: shown.map(token => ({label: token, colour: colours[token], values: (state.series[`Vli_${token}`] ?? []).map(v => 1000 * v)})),
      });
      barPlot(parts.flows, {
        title: `blood flow per lobe at day ${day}`, unit: "ml/min",
        bars: shown.map(token => ({
          label: token, parts: ["hepatic artery", "portal vein"], colour: colours[token],
          values: [1000 * (sample[`Qha_${token}`] ?? 0), 1000 * (sample[`Qpo_${token}`] ?? 0)],
        })),
      });
      linePlot(parts.remnant, {
        title: "the liver which is there", unit: "fraction of the intact liver", times, marker: state.time,
        curves: [{label: "remnant", colour: "var(--accent)", values: state.series.f_remnant ?? []}],
      });
      // the exposure of a lobe is the area under its own outflow over the days
      // which are shown, the readout a study of a regional intervention takes
      barPlot(parts.exposure, {
        title: `exposure of the lobes over ${trim(endTime() / DAY)} days`, unit: "mM*min",
        bars: shown.map(token => ({
          label: token, parts: ["area under the concentration leaving the lobe"], colour: colours[token],
          values: [exposure(times, state.series[`Cout_li_${token}_drg`] ?? [])],
        })),
      });
    }
    parts.status.textContent = `${state.source === "live" ? "Live run of the compiled model" : "Precomputed reference run"}: ${describe()}.`;
    parts.clock.textContent = `day ${trim(state.time / DAY)}`;
    const span = endTime();
    parts.scrub.max = String(span);
    parts.scrub.step = String(span / Math.max((times.length || 2) - 1, 1));
    parts.scrub.value = String(state.time);
    showPlan();
    parts.provenance.textContent = [
      `${dataset().provenance.solver} of sbmlsim-core ${dataset().provenance.sbmlsim_core_version}, ${state.sex} ${state.species}, model ${of_species().sbml_sha256.slice(0, 12)}.`,
      dataset().provenance.notice,
      anatomy?.attribution,
    ].filter(Boolean).join(" ");
  };

  /** The animal of the species, from the physiology figures of the library.
   *
   * The drawing is data of the page (`lobes-anatomy.json`), so the explorer
   * draws the same animal the documentation draws and carries no second one;
   * a page served without the drawings simply has no figure.
   */
  const showBody = () => {
    const drawn = anatomy?.bodies?.[state.species]?.[state.sex];
    parts.body.hidden = !drawn;
    if (!drawn) { parts.body.replaceChildren(); return; }
    const svg = svgNode("svg", {
      class: "body-figure", viewBox: drawn.box.join(" "), role: "img",
      "aria-label": `the ${state.sex} ${state.species} of this run, its liver marked`,
    });
    const group = svgNode("g", {});
    group.innerHTML = drawn.svg;
    svg.append(group);
    // the liver of the animal, the organ this page splits into lobes
    const reach = Math.max(drawn.box[2], drawn.box[3]) / 14;
    svg.append(
      svgNode("circle", {class: "body-liver", cx: drawn.liver[0], cy: drawn.liver[1], r: reach}),
      svgNode("circle", {class: "body-liver-dot", cx: drawn.liver[0], cy: drawn.liver[1], r: reach / 4}),
    );
    parts.body.replaceChildren(svg);
  };

  /** The animal of a species as the icon of its button.
   *
   * The drawing is the one of the dashboard (`lobes-anatomy.json`), at the
   * sex which is picked where the database has it: the page carries no second
   * silhouette, and without the drawings a button is its name alone.
   */
  const speciesIcon = species => {
    const bodies = anatomy?.bodies?.[species];
    const drawn = bodies?.[state.sex] ?? bodies?.male ?? Object.values(bodies ?? {})[0];
    if (!drawn) return null;
    const svg = svgNode("svg", {class: "species-icon", viewBox: drawn.box.join(" "), "aria-hidden": "true", focusable: "false"});
    const group = svgNode("g", {});
    group.innerHTML = drawn.svg;
    svg.append(group);
    return svg;
  };

  /** The volume of the intact liver of the species, the area the frame holds. */
  const intactVolume = () => {
    const intact = of_species().cases.intact?.series;
    if (!intact) return 0;
    return tokens().reduce((sum, token) => sum + (intact[`Vli_${token}`]?.[0] ?? 0), 0);
  };

  /** The planner: the procedure, the lobes it acts on and what it costs. */
  const showPlan = () => {
    const summary = planSummary(of_species().lobes, state.states);
    for (const button of parts.kind.children ?? []) {
      const kind = button.getAttribute?.("data-kind");
      if (kind) button.setAttribute("aria-pressed", String(kind === state.plan.kind));
    }
    for (const chip of parts.lobes.children ?? []) {
      const token = chip.getAttribute?.("data-lobe");
      if (!token) continue;
      const chosen = state.states[token] ?? "intact";
      chip.setAttribute("aria-pressed", String(chosen !== "intact"));
      chip.setAttribute("data-state", chosen);
    }
    const procedure = PROCEDURES.find(entry => entry.id === summary.kind);
    const what = summary.kind === undefined
      ? "No procedure: the intact liver."
      : summary.kind === "mixed"
        ? `Mixed procedures on ${summary.lobes.join(", ")}.`
        : `${procedure.button[0].toUpperCase()}${procedure.button.slice(1)} of ${summary.lobes.join(", ")}.`;
    parts.summary.textContent = `${what} It acts on ${percent(summary.affected)} of the liver mass and leaves a remnant of ${percent(summary.remnant)}.`;
  };

  const percent = share => `${trim(Math.round(1000 * share) / 10)} %`;

  const describe = () => {
    const changed = Object.entries(state.states).filter(([, value]) => value !== "intact");
    if (changed.length === 0) return "the intact liver";
    return changed.map(([token, value]) => `${token} ${STATES.find(s => s.id === value).short}`).join(", ");
  };

  const render = async ({force = false, replay = true} = {}) => {
    const name = matching();
    if (!force && name !== undefined && defaults()) {
      state.series = of_species().cases[name].series;
      state.source = "precomputed";
      state.ran = null;
      parts.runtime.textContent = "";
      if (replay) restart(); else show();
      return;
    }
    // a click while a run is on its way is not lost: it is run once the
    // running one is back, so the page ends on the state the reader chose
    if (state.running) {
      state.again = true;
      return;
    }
    state.running = true;
    busy(true);
    parts.status.textContent = "Solving the lobed liver with SUNDIALS IDA...";
    const started = clock();
    try {
      solvers[state.sex] ??= await solverFactory({
        species: SPECIES.filter(name => sexesOf(name).includes(state.sex)),
        sex: state.sex,
      });
      const result = await solvers[state.sex].simulate(state.species, state.controls, changesOfRun(), {simulation: state.simulation});
      state.series = result.series;
      state.source = "live";
      state.ran = clock() - started;
      parts.runtime.textContent = `ran in ${duration(state.ran)}`;
      if (replay) restart(); else show();
    } catch (error) {
      parts.runtime.textContent = "";
      parts.status.textContent = `The live run failed: ${error.message} The page keeps the last result.`;
    } finally {
      state.running = false;
      busy(false);
    }
    if (state.again) {
      state.again = false;
      await render({replay});
    }
  };

  /** A run which is back is played from its first day.
   *
   * The page has no run button: a change of the plan, of the physiology or of
   * the settings runs the model and the course it produced is played at once,
   * which is what an intervention looks like. A reader who takes the scrub
   * bar pauses it and stays where they put it.
   */
  const restart = () => {
    state.time = 0;
    show();
    if (!state.playing) play();
  };

  /** The spinner says whether the model is running. */
  const busy = running => { parts.spinner.hidden = !running; };

  const clock = () => (globalThis.performance?.now ? globalThis.performance.now() : Date.now());

  /** Walk the days of the course, which is what an intervention looks like. */
  const advance = elapsed => {
    if (!state.playing) return;
    const end = endTime();
    const next = state.time + (elapsed / 1000) * PLAY_DAYS_PER_SECOND * DAY;
    if (next >= end) {
      state.time = end;
      pause();
      show();
      return;
    }
    state.time = next;
    show();
  };

  const play = () => {
    if (state.playing) return pause();
    if (state.time >= endTime()) state.time = 0;
    state.playing = true;
    parts.play.setAttribute("aria-pressed", "true");
    parts.play.textContent = "Pause";
    let last = clock();
    const step = () => {
      if (!state.playing) return;
      const now = clock();
      advance(now - last);
      last = now;
      if (state.playing) frame = schedule(step);
    };
    frame = schedule(step);
    show();
  };

  const pause = () => {
    state.playing = false;
    if (frame !== null) cancel(frame);
    frame = null;
    parts.play.setAttribute("aria-pressed", "false");
    parts.play.textContent = "Play";
  };

  const schedule = step => (globalThis.requestAnimationFrame ? globalThis.requestAnimationFrame(step) : null);
  const cancel = handle => { if (handle !== null && globalThis.cancelAnimationFrame) globalThis.cancelAnimationFrame(handle); };

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
    parts.scrub.addEventListener("input", () => { pause(); state.time = Number(parts.scrub.value); show(); });
    parts.play.addEventListener("click", () => play());
    // the views are two readings of one run, so a switch never runs the model
    parts.views.replaceChildren(...VIEWS.map(view => {
      const button = element("button", {
        class: "tab", type: "button", "data-view": view.id,
        "aria-pressed": String(view.id === state.view), title: view.description,
      }, view.label);
      button.addEventListener("click", () => {
        state.view = view.id;
        for (const other of parts.views.children) other.setAttribute("aria-pressed", String(other.getAttribute("data-view") === view.id));
        show();
      });
      return button;
    }));
    parts.reset.addEventListener("click", () => {
      reset();
      state.plan = {kind: "portal"};
      state.controls = normalizeControls({});
      state.simulation = normalizeSimulation({});
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
    // the sex is a subselection of the species: the parameter database has the
    // female of the human and the male alone of the rat and of the mouse, so
    // the row below the animals offers what the species which is picked has
    parts.sex.replaceChildren(...SEXES.map(name => {
      const has = sexesOf(state.species).includes(name);
      const button = element("button", {
        class: name === state.sex ? "tab active" : "tab", type: "button",
        "aria-pressed": String(name === state.sex),
        title: has
          ? `the reference ${name} ${state.species} of the parameter database`
          : `the parameter database has no reference physiology of the ${name} ${state.species}`,
      }, name);
      if (!has) button.setAttribute("disabled", "");
      else button.addEventListener("click", () => void setSex(name));
      return button;
    }));
    // a species is picked by its animal, the drawing of the physiology figures
    // the dashboard shows beside the liver; a page served without the drawings
    // has the name alone
    parts.species.replaceChildren(...SPECIES.map(name => {
      const button = element("button", {
        class: name === state.species ? "species active" : "species", type: "button",
        "data-species": name, "aria-pressed": String(name === state.species),
        title: `the reference ${name} of the parameter database`,
      });
      const icon = speciesIcon(name);
      if (icon) button.append(icon);
      button.append(element("span", {class: "species-name"}, name));
      button.addEventListener("click", () => void setSpecies(name));
      return button;
    }));
    // the procedures of the planner and the lobes each one can act on
    parts.kind.replaceChildren(...PROCEDURES.map(procedure => {
      const button = element("button", {
        class: "tab operation", type: "button", "data-kind": procedure.id, "data-state": procedure.id,
        "aria-pressed": String(procedure.id === state.plan.kind), title: procedure.label,
      });
      button.append(element("span", {class: "operation-shade"}), element("span", {}, procedure.button));
      button.addEventListener("click", () => {
        // the lobes already chosen move to the operation which was picked, and
        // the intact liver is the operation which puts every lobe back
        const chosen = Object.entries(state.states).filter(([, value]) => value !== "intact").map(([token]) => token);
        state.plan = {kind: procedure.id};
        state.states = procedure.id === "intact"
          ? Object.fromEntries(tokens().map(token => [token, "intact"]))
          : {...state.states, ...Object.fromEntries(chosen.map(token => [token, procedure.id]))};
        void render();
      });
      return button;
    }));
    parts.lobes.replaceChildren(...of_species().lobes.map(lobe => {
      const chip = element("button", {
        class: "lobe-chip", type: "button", "data-lobe": lobe.token,
        "aria-pressed": "false", title: `${lobe.name} (${trim(100 * lobe.mass_fraction)} % of the liver mass)`,
      });
      chip.append(
        element("span", {class: "lobe-dot", style: `background:${colours[lobe.token]}`}),
        element("span", {class: "lobe-chip-token"}, lobe.token),
      );
      chip.addEventListener("click", () => {
        const now = state.states[lobe.token] ?? "intact";
        state.states = {...state.states, [lobe.token]: now === "intact" ? state.plan.kind : "intact"};
        void render();
      });
      return chip;
    }));
    // the legend names the lobes by their colour and the states by their shade
    const lobeItems = of_species().lobes.map(lobe => {
      const item = element("span", {class: "legend-item", title: lobe.name});
      item.append(
        element("span", {class: "legend-swatch", style: `background:${colours[lobe.token]}`}),
        element("span", {class: "legend-label"}, lobe.token),
      );
      return item;
    });
    const stateItems = STATES.filter(entry => entry.id !== "resected").map(entry => {
      const item = element("span", {class: "legend-item", "data-state": entry.id, title: entry.description});
      item.append(element("span", {class: "legend-shade"}), element("span", {class: "legend-label"}, entry.label));
      return item;
    });
    parts.legend.replaceChildren(...lobeItems, element("span", {class: "legend-break"}), ...stateItems);
    // the settings of one simulation: how long it runs, how finely it is read
    // out and the dose of the marker drug the exposure is read from
    parts.sim.replaceChildren(...Object.entries(SIMULATION).map(([key, setting]) => {
      const row = element("label", {class: "sim-setting", title: setting.label});
      const select = element("select", {class: "sim-select", "data-setting": key, "aria-label": setting.label});
      for (const choice of setting.choices) {
        const option = element("option", {value: String(choice)}, String(choice));
        if (choice === state.simulation[key]) option.setAttribute("selected", "");
        select.append(option);
      }
      select.value = String(state.simulation[key]);
      select.addEventListener("change", () => {
        state.simulation = normalizeSimulation({...state.simulation, [key]: Number(select.value)});
        void render();
      });
      row.append(element("span", {class: "sim-label"}, setting.short ?? setting.label), select);
      return row;
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

  const setSex = async name => {
    // the parameter database has no female rat and no female mouse, so the
    // sex of the page can only be one the species of the page has
    if (name === state.sex || !sexesOf(state.species).includes(name)) return render();
    pause();
    datasets[name] ??= await loadData(name);
    state.sex = name;
    return setSpecies(state.species);
  };

  const setSpecies = name => {
    pause();
    state.species = name;
    // a species the sex of the page has no physiology of falls back to the
    // male, which every species of the database has
    if (!sexesOf(name).includes(state.sex)) state.sex = "male";
    colours = lobeColors(of_species().lobes);
    reset();
    state.time = 0;
    parts.scrub.value = "0";
    liver = createLiver(parts.liver, of_species().lobes, {onPick: pick});
    showBody();
    build();
    buildPresets();
    return render();
  };

  /** The named plans of the species, grouped by the operation they are.
   *
   * A reader who looks for a hepatectomy should not read a row of names to
   * find one: the plans sit under the operation they belong to, in the order
   * of the operations and with the shade of that operation, so the row of the
   * planner, the legend of the drawing and these buttons agree.
   */
  const buildPresets = () => {
    const groups = new Map([...PROCEDURES.map(procedure => [procedure.id, []]), ["mixed", []]]);
    for (const [name, entry] of Object.entries(of_species().cases)) {
      const states = stateOfCase(entry, tokens());
      const kind = planSummary(of_species().lobes, states).kind ?? "intact";
      (groups.get(kind) ?? groups.get("mixed")).push({name, entry, states});
    }
    const nodes = [];
    for (const [kind, cases] of groups) {
      if (cases.length === 0) continue;
      const procedure = PROCEDURES.find(entry => entry.id === kind);
      const group = element("div", {class: "preset-group", "data-state": kind});
      group.append(
        element("span", {class: "preset-group-label", title: procedure?.label ?? "several operations at once"}, procedure?.group ?? "mixed"),
      );
      for (const {name, entry, states} of cases) {
        const button = element("button", {
          class: "preset", type: "button", "data-state": kind, "data-case": name,
          title: `${entry.label}: the reference run of ${name}`,
        }, presetLabel(entry.label, kind));
        button.addEventListener("click", () => {
          state.states = states;
          const picked = planSummary(of_species().lobes, state.states).kind;
          state.plan = {kind: picked === undefined || picked === "mixed" ? state.plan.kind : picked};
          void render();
        });
        group.append(button);
      }
      nodes.push(group);
    }
    parts.presets.replaceChildren(...nodes);
  };

  colours = lobeColors(of_species().lobes);
  reset();
  wire();
  build();
  buildPresets();
  liver = createLiver(parts.liver, of_species().lobes, {onPick: pick});
  showBody();
  return {
    state,
    pick,
    render,
    run: () => render({force: true}),
    play,
    pause,
    advance,
    setSpecies,
    setSex,
    setTime(time) { pause(); state.time = time; show(); },
    close() { pause(); Object.values(solvers).forEach(solver => solver?.close()); },
  };
}

/** Start the explorer in the page which hosts it. */
/** The dataset of a sex: the male beside the page, the female under `female/`. */
export async function defaultLoadData(sex) {
  const path = sex === "male" ? "lobes-data.json" : `${sex}/lobes-data.json`;
  const response = await fetch(new URL(path, import.meta.url));
  if (!response.ok) throw new Error(`Could not load the ${sex} reference data: HTTP ${response.status}.`);
  return response.json();
}

export async function start(root = document) {
  const response = await fetch(new URL("lobes-data.json", import.meta.url));
  if (!response.ok) throw new Error(`Could not load the reference data: HTTP ${response.status}.`);
  const data = await response.json();
  // the drawings are a figure beside the numbers: a page without them works
  const drawings = await fetch(new URL("lobes-anatomy.json", import.meta.url)).catch(() => null);
  const anatomy = drawings?.ok ? await drawings.json() : null;
  const app = createApp(root, {data, anatomy});
  await app.render();
  // the explorer is both a page of its own and a figure of the documentation:
  // on its own it takes the window, in a frame it is as tall as it needs to be
  document.documentElement.dataset.framed = String(globalThis.parent !== undefined && globalThis.parent !== globalThis);
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
