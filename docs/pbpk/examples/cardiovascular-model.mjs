// Exact reduction of the steady pressure network. Dynamic traces are loaded
// separately from the reference solutions generated with sbmlsim/IDA.
export function defaultControls(reference) {
  return {hepatic: 1, systemic: 1, portal: 1, renal: 1,
    arterial: reference.arterial_pressure, venous: reference.venous_pressure};
}

export function solvePressure(reference, controls) {
  const {hepatic, systemic, portal, renal, arterial, venous} = controls;
  if (![hepatic, systemic, portal, renal, arterial, venous].every(Number.isFinite)
      || Math.min(hepatic, systemic, portal, renal) <= 0 || arterial <= venous) {
    throw new RangeError("Positive resistances and arterial pressure above venous pressure are required.");
  }
  const conductances = reference.organs.map(organ => ({...organ,
    conductance: 1 / (organ.resistance * systemic * (organ.id === "ki" ? renal : 1))}));
  const splanchnic = conductances.filter(organ => organ.portal)
    .reduce((sum, organ) => sum + organ.conductance, 0);
  const hepaticArtery = conductances.find(organ => organ.id === "li").conductance;
  const gp = 1 / (reference.portal_resistance * portal);
  const gh = 1 / (reference.hepatic_resistance * hepatic);
  // Nodal conservation, with potentials shifted by venous pressure to avoid
  // cancellation when calculating small pressure differences.
  const a = splanchnic + gp;
  const d = gp + hepaticArtery + gh;
  const determinant = splanchnic * d + gp * (hepaticArtery + gh);
  const head = arterial - venous;
  const ppo = venous + head * (splanchnic * d + gp * hepaticArtery) / determinant;
  const pli = venous + head * (hepaticArtery * a + gp * splanchnic) / determinant;
  const flows = conductances.map(organ => ({...organ,
    flow: organ.conductance * (arterial - (organ.portal ? ppo : organ.id === "li" ? pli : venous))}));
  const qpo = (ppo - pli) * gp;
  const qh = (pli - venous) * gh;
  const qha = flows.find(organ => organ.id === "li").flow;
  const co = flows.reduce((sum, organ) => sum + organ.flow, 0);
  return {ppo, pli, qpo, qh, qha, co, flows,
    systemicFlow: flows.filter(organ => !organ.portal && organ.id !== "li")
      .reduce((sum, organ) => sum + organ.flow, 0),
    portalResidual: flows.filter(organ => organ.portal)
      .reduce((sum, organ) => sum + organ.flow, 0) - qpo,
    hepaticResidual: qha + qpo - qh};
}

export function column(series, name) {
  const index = series.columns.indexOf(name);
  if (index < 0) throw new Error(`Missing reference column: ${name}`);
  return series.values.map(row => row[index]);
}

export function lastBeat(series, heartRate) {
  const timeIndex = series.columns.indexOf("time");
  const end = series.values.at(-1)[timeIndex];
  return {...series, values: series.values.filter(row => row[timeIndex] >= end - 60 / heartRate - 1e-8)};
}

export function sampleAt(series, time) {
  const timeIndex = series.columns.indexOf("time");
  const rows = series.values;
  let low = 0, high = rows.length;
  while (low < high) {
    const middle = Math.floor((low + high) / 2);
    if (rows[middle][timeIndex] <= time) low = middle + 1;
    else high = middle;
  }
  const left = rows[Math.max(0, low - 1)];
  const right = rows[Math.min(rows.length - 1, low)];
  const span = right[timeIndex] - left[timeIndex];
  const weight = span > 0 ? Math.max(0, Math.min(1, (time - left[timeIndex]) / span)) : 0;
  return Object.fromEntries(series.columns.map((name, index) => [name, left[index] + weight * (right[index] - left[index])]));
}

export function csv(series, units = {}) {
  const quote = value => `"${String(value).replaceAll('"', '""')}"`;
  return [series.columns.map(name => quote(units[name] ? `${name} [${units[name]}]` : name)).join(","),
    ...series.values.map(row => row.join(","))].join("\n") + "\n";
}
