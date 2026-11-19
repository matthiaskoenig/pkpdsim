// Theme-aware inline SVG, using the .core-diagram / .cd-* vocabulary of
// sbmlsim-core's explanatory diagrams. Geometry is fixed; simulated pressures
// and flows update the nodes and edges without rebuilding the diagram.
const SVG = "http://www.w3.org/2000/svg";
const SYSTEMIC = ["Qhe", "Qbr", "Qmu", "Qsk", "Qad", "Qbo", "Qgo"];
const SPLANCHNIC = ["Qst", "Qsi", "Qla", "Qsp", "Qpa"];

export function flowWidth(flow) {
  if (!Number.isFinite(flow)) return 0.7;
  return 0.7 + 1.7 * Math.sqrt(Math.abs(flow));
}

export function schematicValues(row) {
  const sum = keys => keys.reduce((total, key) => total + (row[key] ?? 0), 0);
  return {...row, systemic: sum(SYSTEMIC), splanchnic: sum(SPLANCHNIC)};
}

function node(tag, attributes = {}, text) {
  const element = document.createElementNS(SVG, tag);
  for (const [key, value] of Object.entries(attributes)) element.setAttribute(key, value);
  if (text !== undefined) element.textContent = text;
  return element;
}

const number = value => Number.isFinite(value) ? value.toFixed(2) : "-";

export function createSchematic(host, dynamic = false) {
  const offset = dynamic ? 160 : 0;
  const id = dynamic ? "heart" : "pressure";
  const svg = node("svg", {
    viewBox: `0 0 1000 ${dynamic ? 690 : 530}`,
    class: "core-diagram cv-schematic",
    role: "img",
    "aria-labelledby": `${id}-title ${id}-description`,
  });
  svg.append(node("title", {id: `${id}-title`}, dynamic ? "Closed-loop cardiovascular circulation" : "Whole-body pressure-driven circulation"));
  svg.append(node("desc", {id: `${id}-description`},
    "Directed blood flow through systemic tissues, kidney, portal organs and the dual hepatic supply. " +
    (dynamic ? "The return loop includes both ventricles, four valves and the pulmonary circulation. " : "Arterial and venous pressures are prescribed. ") +
    "Edge width follows the square root of instantaneous flow on one shared scale. Labels give pressures in millimeters of mercury and flows in liters per minute."));
  const defs = node("defs");
  for (const color of ["arterial", "venous", "portal", "hepatic", "pulmonary"]) {
    const marker = node("marker", {id: `${id}-${color}-arrow`, viewBox: "0 0 10 10", refX: 9, refY: 5,
      markerWidth: 9, markerHeight: 9, orient: "auto-start-reverse", markerUnits: "userSpaceOnUse"});
    marker.append(node("path", {d: "M0 0 L10 5 L0 10 z", class: `cv-arrow cv-${color}`}));
    defs.append(marker);
  }
  svg.append(defs);
  const edges = node("g", {class: "cv-edges"});
  const nodes = node("g", {class: "cv-nodes"});
  const valves = node("g", {class: "cv-valves"});
  const edgeElements = [];
  const nodeElements = [];
  const y = value => value + offset;

  function edge(name, flow, path, color, label) {
    const group = node("g", {"data-edge": name});
    group.append(node("path", {d: path, class: "cv-edge-track"}));
    const line = node("path", {d: path, class: `cv-edge cv-${color}`, "marker-end": `url(#${id}-${color}-arrow)`});
    const title = node("title", {}, name);
    line.append(title);
    group.append(line);
    let value;
    if (label) {
      value = node("text", {x: label[0], y: label[1], class: "cd-sub cv-flow-value"});
      group.append(value);
    }
    edges.append(group);
    edgeElements.push({name, flow, line, title, value});
  }

  function compartment(name, title, x, cy, key, unit = "mmHg", kind = "", width = 132, volume) {
    const group = node("g", {class: `cd-step cd-module cv-compartment ${kind}`, "data-node": name});
    group.append(node("rect", {x: x - width / 2, y: cy - 37, width, height: 74, rx: 7}));
    group.append(node("text", {x, y: cy - 13, class: "cv-node-title"}, title));
    const value = node("text", {x, y: cy + 12, class: "cv-node-value"});
    group.append(value);
    const sub = node("text", {x, y: cy + 30, class: "cd-sub"});
    if (volume) group.append(sub);
    nodes.append(group);
    nodeElements.push({value, sub, key, unit, volume});
  }

  function valve(name, cx, cy, direction) {
    const group = node("g", {class: "cv-valve"});
    const rotation = direction === "up" ? 180 : direction === "left" ? 90 : 0;
    group.append(node("path", {d: "M-8 -7 L0 3 L8 -7 M-8 8 H8", transform: `translate(${cx} ${cy}) rotate(${rotation})`}));
    group.append(node("text", {x: direction === "left" ? cx : cx > 500 ? cx - 80 : cx + 80,
      y: cy + (direction === "left" ? 56 : 5), class: "cd-sub cv-valve-label"}, name));
    valves.append(group);
  }

  edge("Artery to systemic tissues", "systemic", `M154 ${y(194)} H225 V${y(80)} H395`, "arterial", [305, y(65)]);
  edge("Systemic tissue return", "systemic", `M585 ${y(80)} H850 V${y(194)} H844`, "venous", [745, y(65)]);
  edge("Renal arterial supply", "Qki", `M154 ${y(210)} H410`, "arterial", [285, y(195)]);
  edge("Renal venous return", "Qki", `M570 ${y(210)} H844`, "venous", [705, y(195)]);
  edge("Splanchnic arterial supply", "splanchnic", `M80 ${y(247)} V${y(400)} H210`, "arterial", [142, y(385)]);
  edge("Portal organ drainage", "splanchnic", `M390 ${y(400)} H480`, "portal", [435, y(385)]);
  edge("Portal vein to sinusoids", "Qpo", `M640 ${y(400)} H675 V${y(370)} H710`, "portal", [690, y(431)]);
  edge("Hepatic arterial supply", "Qha", `M154 ${y(230)} H200 V${y(350)} H710`, "hepatic", [455, y(335)]);
  edge("Hepatic venous return", "Qh", `M870 ${y(350)} H920 V${y(247)}`, "hepatic", [922, y(385)]);

  compartment("ar", dynamic ? "Systemic artery" : "Arterial boundary", 80, y(210), "CV_P_ar", "mmHg", "cv-arterial-node", 148);
  compartment("ve", dynamic ? "Systemic vein" : "Venous boundary", 920, y(210), "CV_P_ve", "mmHg", "cv-venous-node", 148);
  compartment("systemic", "Systemic tissues", 490, y(80), "systemic", "l/min", "", 190);
  compartment("kidney", "Kidney", 490, y(210), "Qki", "l/min", "", 160);
  compartment("gut", "Portal organs", 300, y(400), "splanchnic", "l/min", "", 180);
  compartment("po", "Portal vein", 560, y(400), "CV_P_po", "mmHg", "cv-portal-node", 160);
  compartment("li", "Hepatic sinusoids", 790, y(350), "CV_P_li", "mmHg", "cv-hepatic-node", 160);

  if (dynamic) {
    edge("Tricuspid valve", "CV_Q_tv", `M920 ${y(173)} V107`, "venous", [963, 222]);
    edge("Pulmonary valve", "CV_Q_pv", "M854 70 H792", "pulmonary", [815, 51]);
    edge("Pulmonary artery to lung", "Qlu", "M644 70 H566", "pulmonary", [605, 51]);
    edge("Lung to pulmonary vein", "Qlu", "M434 70 H372", "pulmonary", [395, 51]);
    edge("Mitral valve", "CV_Q_mv", "M224 70 H146", "venous", [185, 51]);
    edge("Aortic valve", "CO", `M80 107 V${y(173)}`, "arterial", [35, 222]);
    compartment("rv", "Right ventricle", 920, 70, "CV_P_rv", "mmHg", "cv-venous-node", 132, "CV_V_rv");
    compartment("pa", "Pulmonary artery", 710, 70, "CV_P_pa", "mmHg", "", 160);
    compartment("lu", "Lung", 500, 70, "Qlu", "l/min", "", 132);
    compartment("pv", "Pulmonary vein", 290, 70, "CV_P_pv", "mmHg", "", 160);
    compartment("lv", "Left ventricle", 80, 70, "CV_P_lv", "mmHg", "cv-arterial-node", 132, "CV_V_lv");
    valve("Tricuspid", 920, 165, "up");
    valve("Pulmonary", 815, 70, "left");
    valve("Mitral", 185, 70, "left");
    valve("Aortic", 80, 165, "down");
    svg.append(node("text", {x: 500, y: 145, class: "cd-title"}, "PULMONARY CIRCULATION AND HEART"));
  }
  svg.append(node("text", {x: 500, y: y(28), class: "cd-title"}, "SYSTEMIC AND PORTAL CIRCULATION"));
  svg.append(edges, nodes, valves);
  const legend = node("g", {class: "cv-width-legend"});
  const legendY = dynamic ? 657 : 497;
  legend.append(node("text", {x: 20, y: legendY, class: "cd-label", "text-anchor": "start"}, "One width scale:"));
  for (const [index, flow] of [0, 1, 5, 20, 100].entries()) {
    const x = 220 + index * 120;
    legend.append(node("line", {x1: x - 22, x2: x + 22, y1: legendY - 5, y2: legendY - 5, class: "cv-legend-line", "stroke-width": flowWidth(flow)}));
    legend.append(node("text", {x, y: legendY + 20, class: "cd-sub"}, String(flow)));
  }
  legend.append(node("text", {x: 900, y: legendY + 20, class: "cd-sub"}, "l/min"));
  svg.append(legend);
  host.replaceChildren(svg);
  return {
    update(input) {
      const values = schematicValues(input);
      for (const edge of edgeElements) {
        const flow = values[edge.flow];
        const flowing = Number.isFinite(flow) && Math.abs(flow) > 1e-6;
        edge.line.setAttribute("stroke-width", flowWidth(flow));
        edge.line.setAttribute("stroke-opacity", flowing ? 0.88 : 0.18);
        edge.line.setAttribute("data-flow", Number.isFinite(flow) ? flow : "");
        edge.title.textContent = `${edge.name}: ${number(flow)} l/min`;
        if (edge.value) edge.value.textContent = number(flow);
      }
      for (const item of nodeElements) {
        item.value.textContent = `${number(values[item.key])} ${item.unit}`;
        if (item.volume) item.sub.textContent = `${number(values[item.volume])} ml`;
      }
    },
  };
}
