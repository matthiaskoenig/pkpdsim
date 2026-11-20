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
