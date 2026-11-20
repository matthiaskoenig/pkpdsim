// The pure helpers of the liver lobe explorer: the states of a lobe and the
// parameters they change, the columns of a precomputed case and the areas,
// samples and text a reader downloads. No model of a lobed liver lives here;
// every number comes from the compiled SBML, run by SUNDIALS IDA in a Web
// Worker (lobes-solver.mjs) or precomputed by scripts/lobes_examples.py.

export const STATES = [
  {id: "intact", label: "intact", short: "intact", description: "Both branches open: the lobe is perfused as the lobe table says."},
  {id: "portal", label: "portal vein ligated", short: "PVL", description: "The portal branch is closed: the lobe keeps its artery, stays perfused below the reference and atrophies."},
  {id: "arterial", label: "hepatic artery ligated", short: "HAL", description: "The arterial branch is closed: the lobe keeps its portal inflow."},
  {id: "resected", label: "resected", short: "resected", description: "Both branches closed and the lobe gone: no capacity, no growth, and the remnant is what is left."},
];

const PARAMETERS = {
  intact: () => ({}),
  portal: token => ({[`g_po_${token}`]: 0}),
  arterial: token => ({[`g_ha_${token}`]: 0}),
  resected: token => ({[`g_ha_${token}`]: 0, [`g_po_${token}`]: 0, [`f_present_${token}`]: 0}),
};

/** The parameter changes of a state per lobe, as `pkpdsim.pbpk.interventions` writes them. */
export function caseChanges(states) {
  const changes = {};
  for (const [token, state] of Object.entries(states)) {
    const of_state = PARAMETERS[state];
    if (!of_state) throw new RangeError(`'${state}' is no state of a lobe.`);
    Object.assign(changes, of_state(token));
  }
  return changes;
}

/** The state per lobe of a precomputed case, so a preset fills the clickable liver. */
export function stateOfCase(entry, tokens) {
  const changes = entry?.changes ?? {};
  return Object.fromEntries(tokens.map(token => {
    if (`f_present_${token}` in changes) return [token, "resected"];
    if (`g_po_${token}` in changes && `g_ha_${token}` in changes) return [token, "resected"];
    if (`g_po_${token}` in changes) return [token, "portal"];
    if (`g_ha_${token}` in changes) return [token, "arterial"];
    return [token, "intact"];
  }));
}

/** The trapezoid area under a course, the exposure of a lobe over the shown days. */
export function exposure(times, values) {
  let area = 0;
  for (let i = 1; i < times.length; i++) area += 0.5 * (times[i] - times[i - 1]) * (values[i] + values[i - 1]);
  return area;
}

/** The row of a course at or before a time, every column of it by its name. */
export function sampleAt(series, time) {
  const times = series?.time ?? [];
  if (times.length === 0) return {};
  let row = 0;
  while (row + 1 < times.length && times[row + 1] <= time + 1e-9) row++;
  return Object.fromEntries(Object.entries(series).map(([column, values]) => [column, values[row]]));
}

/** The columns of a course as comma separated text, the download of the page. */
export function csv(series, columns) {
  for (const column of columns) if (!(column in series)) throw new RangeError(`The course has no column '${column}'.`);
  const rows = series[columns[0]].length;
  const lines = [columns.join(",")];
  for (let row = 0; row < rows; row++) lines.push(columns.map(column => series[column][row]).join(","));
  return `${lines.join("\n")}\n`;
}

/** A number for a label: four significant digits, an exponent when it is small. */
export function trim(value) {
  if (!Number.isFinite(value)) return "-";
  if (value === 0) return "0";
  const magnitude = Math.abs(value);
  if (magnitude >= 1e-3 && magnitude < 1e6) return String(Number(value.toPrecision(4)));
  return value.toExponential(2).replace("e", "e").replace("e+", "e");
}

export const LOBE_SLOTS = 8;

/** The colour of every lobe of a set, by its place in the lobe table.
 *
 * One slot of the palette per lobe (`--lobe-1` to `--lobe-8` of lobes.css,
 * stepped for the light and for the dark theme), keyed by the order of the
 * table and never by the order of what is plotted: a lobe keeps its colour
 * when another one is resected and leaves the diagram. A set longer than the
 * palette gets the neutral slot for the lobes past it rather than a repeated
 * hue, which would say two lobes are one.
 */
export function lobeColors(lobes) {
  return Object.fromEntries(lobes.map((lobe, index) => [
    lobe.token,
    index < LOBE_SLOTS ? `var(--lobe-${index + 1})` : "var(--lobe-other)",
  ]));
}

/** What a plan does to the liver: the mass it covers and the mass left behind.
 *
 * `affected` is the share of the liver mass the procedure acts on, which is
 * what the names of the procedures of the literature count (a `pvl70`
 * deportalises 70 % of the mass); `remnant` is the share which is still in the
 * liver afterwards, so a ligation leaves it at 1 and a resection takes the
 * mass out. `kind` is the single state of the plan, `"mixed"` when the lobes
 * carry more than one and `undefined` for the intact liver.
 */
export function planSummary(lobes, states) {
  const mass = Object.fromEntries(lobes.map(lobe => [lobe.token, lobe.mass_fraction]));
  let affected = 0, remnant = 0;
  const kinds = new Set();
  const chosen = [];
  for (const lobe of lobes) {
    const state = states?.[lobe.token] ?? "intact";
    const share = mass[lobe.token] ?? 0;
    if (state === "intact") { remnant += share; continue; }
    kinds.add(state);
    chosen.push(lobe.token);
    affected += share;
    if (state !== "resected") remnant += share;
  }
  const kind = kinds.size === 0 ? undefined : kinds.size === 1 ? [...kinds][0] : "mixed";
  return {affected, remnant, kind, lobes: chosen};
}

/** A run time for a reader: milliseconds below a second, seconds above. */
export function duration(ms) {
  if (!Number.isFinite(ms)) return "-";
  if (ms < 1000) return `${Math.round(ms)} ms`;
  return `${Number((ms / 1000).toPrecision(2))} s`;
}

/** The two readings of one run: what grows, and what drives it.
 *
 * A view is no second simulation: the series of a run carries the volumes and
 * the pressures alike, so the switch only changes what the dashboard reads
 * out of the result which is already there.
 */
export const VIEWS = [
  {id: "volume", label: "volumes", description: "What the remodelling does: the volume of every lobe, the remnant and the exposure."},
  {id: "pressure", label: "pressures and flows", description: "What drives it: the pressure of the inflow, the pressure of every lobe and its perfusion per volume."},
];
