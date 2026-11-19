// Anatomy comes from the physiology figure's Body, with the core-diagram SVG
// vocabulary used by sbmlsim-core. Overlays represent the hydraulic topology.
const SVG = "http://www.w3.org/2000/svg";
const SYSTEMIC = ["Qhe", "Qbr", "Qmu", "Qsk", "Qad", "Qbo", "Qgo"];
const SPLANCHNIC = ["Qst", "Qsi", "Qla", "Qsp", "Qpa"];
export function flowWidth(flow, scale = 1) {
  return 0.7 + (Number.isFinite(flow) ? 1.7 * scale * Math.sqrt(Math.abs(flow)) : 0);
}
export function schematicValues(row) {
  const sum = keys => keys.reduce((total, key) => total + (row[key] ?? 0), 0);
  return {...row, systemic: sum(SYSTEMIC), splanchnic: sum(SPLANCHNIC),
    arterialTotal: sum([...SYSTEMIC, ...SPLANCHNIC, "Qki", "Qha"]),
    venousTotal: sum([...SYSTEMIC, "Qki", "Qh"])};
}
function node(tag, attributes = {}, text) {
  const element = document.createElementNS(SVG, tag);
  for (const [key, value] of Object.entries(attributes)) element.setAttribute(key, value);
  if (text !== undefined) element.textContent = text;
  return element;
}
const number = value => Number.isFinite(value) ? value.toFixed(2) : "-";
export function createSchematic(host, dynamic = false, anatomy) {
  const id = dynamic ? "heart" : "pressure";
  const svg = node("svg", {viewBox: "0 0 720 630", class: "core-diagram cv-schematic", role: "img", "aria-labelledby": `${id}-title ${id}-description`});
  svg.append(node("title", {id: `${id}-title`}, "Whole-body cardiovascular circulation"));
  svg.append(node("desc", {id: `${id}-description`}, "Male or female anatomy from the physiology figures. Directed arterial, venous and portal flows show the model topology, not vessel anatomy. Width follows square root of absolute instantaneous flow, on one adjustable scale. Pressures: mmHg; flows: liters per minute."));
  const defs = node("defs");
  for (const color of ["arterial", "venous", "portal", "hepatic", "pulmonary"]) {
    const marker = node("marker", {id: `${id}-${color}-arrow`, viewBox: "0 0 10 10", refX: 9, refY: 5, markerWidth: 6, markerHeight: 6, orient: "auto-start-reverse", markerUnits: "userSpaceOnUse"});
    marker.append(node("path", {d: "M0 0 L10 5 L0 10 z", class: `cv-arrow cv-${color}`})); defs.append(marker);
  }
  svg.append(defs);
  const scale = 560 / anatomy.height;
  const body = node("g", {class: "cv-body", transform: `translate(360 20) scale(${scale})`});
  body.innerHTML = anatomy.svg; svg.append(body);
  const anchors = Object.fromEntries(Object.entries(anatomy.anchors).map(([key, [x, y]]) => [key, [360 + scale * x, 20 + scale * y]]));
  const edges = node("g", {class: "cv-edges"}), labels = node("g", {class: "cv-labels"});
  const edgeElements = [], readouts = [], legendLines = [];
  function edge(name, flow, d, color) {
    const group = node("g", {"data-edge": name});
    group.append(node("path", {d, class: "cv-edge-track"}));
    const line = node("path", {d, class: `cv-edge cv-${color}`});
    const title = node("title"); line.append(title); group.append(line); edges.append(group);
    edgeElements.push({name, flow, line, title, color});
  }
  function readout(name, title, x, y, key) {
    const group = node("g", {class: "cv-compartment", "data-node": name});
    group.append(node("rect", {x: x - 75, y: y - 23, width: 150, height: 46, rx: 6}));
    group.append(node("text", {x, y: y - 5, "text-anchor": "middle", class: "cv-node-title"}, title));
    const value = node("text", {x, y: y + 14, "text-anchor": "middle", class: "cv-node-value"});
    group.append(value); labels.append(group); readouts.push({value, key, unit: "mmHg"});
  }
  function organ(key, title, y, right = false) {
    const [x, ay] = anchors[key], lx = right ? 540 : 15;
    labels.append(node("path", {d: `M${x} ${ay} L${right ? 515 : 205} ${y - 3} H${right ? 540 : 170}`, class: "cv-leader"}));
    labels.append(node("text", {x: lx, y, class: "cv-label"}, title));
    const value = node("text", {x: lx, y: y + 16, class: "cd-sub cv-flow-value", "data-organ": key});
    labels.append(value); readouts.push({value, key: key === "li" ? "Qh" : `Q${key}`, unit: "l/min"});
    if (key === "li") return;
    edge(key === "ki" ? "Renal arterial supply" : `${title} arterial supply`, `Q${key}`, `M265 ${ay} Q${(265 + x) / 2} ${ay - 8} ${x} ${ay}`, "arterial");
    const portal = ["st", "si", "la", "sp", "pa"].includes(key);
    edge(key === "ki" ? "Renal venous return" : `${title} ${portal ? "portal drainage" : "venous return"}`, `Q${key}`, portal ? `M${x} ${ay} Q480 ${ay} 480 300` : `M${x} ${ay} Q${(460 + x) / 2} ${ay + 8} 460 ${ay}`, portal ? "portal" : "venous");
  }
  edge("Systemic arterial branches", "arterialTotal", "M265 50 V510", "arterial");
  edge("Systemic venous return", "venousTotal", "M460 510 V50", "venous");
  for (const [key, title, y] of [["br", "Brain", 43], ["he", "Heart tissue", 230], ["ki", "Kidney", 278], ["sk", "Skin", 330], ["ad", "Adipose", 376], ["mu", "Muscle", 422], ["bo", "Bone", 468], ["go", "Gonads", 514]]) organ(key, title, y);
  for (const [key, title, y] of [["st", "Stomach", 230], ["pa", "Pancreas", 268], ["sp", "Spleen", 306], ["si", "Small intestine", 344], ["la", "Large intestine", 382]]) organ(key, title, y, true);
  const [x, y] = anchors.li;
  edge("Portal organ drainage", "splanchnic", "M480 300 V280", "portal");
  edge("Portal vein to sinusoids", "Qpo", `M480 280 V${y + 18} L${x} ${y}`, "portal");
  edge("Hepatic arterial supply", "Qha", `M265 ${y} H${x}`, "hepatic");
  edge("Hepatic venous return", "Qh", `M${x} ${y} Q430 ${y - 10} 460 ${y}`, "hepatic");
  readout("po", "Portal vein", 620, 440, "CV_P_po");
  readout("li", "Sinusoids", 620, 496, "CV_P_li");
  labels.append(node("path", {d: `M545 440 H505 V280 H480 M545 496 H500 V${y} H${x}`, class: "cv-leader"}));
  organ("li", "Total hepatic flow", 550, true);
  if (dynamic) {
    edge("Tricuspid valve", "CV_Q_tv", "M460 160 H520 V85 H545", "venous");
    edge("Pulmonary valve", "CV_Q_pv", "M620 62 V22 H100 V87", "pulmonary");
    const [lx, ly] = anchors.lu;
    edge("Pulmonary artery to lung", "Qlu", `M175 110 H205 L${lx} ${ly}`, "pulmonary");
    edge("Lung to pulmonary vein", "Qlu", `M${lx} ${ly} L205 168 H175`, "arterial");
    edge("Mitral valve", "CV_Q_mv", "M175 168 H215 V204 H520 V145 H545", "arterial");
    edge("Aortic valve", "CO", "M620 168 V210 H265 V180", "arterial");
    readout("rv", "Right ventricle", 620, 85, "CV_P_rv");
    readout("lv", "Left ventricle", 620, 145, "CV_P_lv");
    readout("pa", "Pulmonary artery", 100, 110, "CV_P_pa");
    readout("pv", "Pulmonary vein", 100, 168, "CV_P_pv");
    for (const [label, lx, ly] of [["TV", 500, 149], ["PV", 350, 15], ["MV", 435, 198], ["AV", 540, 225]]) labels.append(node("text", {x: lx, y: ly, class: "cd-sub"}, label));
  } else {
    readout("ar", "Arterial pressure", 100, 120, "CV_P_ar");
    readout("ve", "Venous pressure", 620, 120, "CV_P_ve");
    edge("Arterial boundary", "arterialTotal", "M175 120 H265", "arterial");
    edge("Venous boundary", "venousTotal", "M460 120 H545", "venous");
  }
  svg.append(edges, labels);
  const legend = node("g", {class: "cv-width-legend", transform: "translate(200 610)"});
  for (const [index, flow] of [1, 5, 20].entries()) {
    const lx = index * 110;
    const line = node("line", {x1: lx, x2: lx + 30, y1: 0, y2: 0});
    legend.append(line, node("text", {x: lx + 37, y: 4}, `${flow} l/min`)); legendLines.push({line, flow});
  }
  svg.append(legend); host.replaceChildren(svg);
  if (host.scrollWidth > host.clientWidth) host.scrollLeft = (host.scrollWidth - host.clientWidth) / 2;
  return {update(input, widthScale = 1) {
    const values = schematicValues(input);
    for (const edge of edgeElements) {
      const flow = values[edge.flow], flowing = Number.isFinite(flow) && Math.abs(flow) > 1e-6;
      edge.line.setAttribute("stroke-width", flowWidth(flow, widthScale));
      edge.line.setAttribute("stroke-opacity", flowing ? 0.88 : 0.18);
      edge.line.setAttribute("data-flow", Number.isFinite(flow) ? flow : "");
      edge.line.setAttribute("marker-end", flow < 0 ? "none" : `url(#${id}-${edge.color}-arrow)`);
      edge.line.setAttribute("marker-start", flow < 0 ? `url(#${id}-${edge.color}-arrow)` : "none");
      edge.title.textContent = `${edge.name}: ${number(flow)} l/min`;
    }
    for (const item of readouts) {
      const value = values[item.key];
      const display = item.unit === "l/min" && Number.isFinite(value) && Math.abs(value) < 0.01 ? value.toFixed(3) : number(value);
      item.value.textContent = `${display} ${item.unit}`;
    }
    for (const item of legendLines) item.line.setAttribute("stroke-width", flowWidth(item.flow, widthScale));
  }};
}
