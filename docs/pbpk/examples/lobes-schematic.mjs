// The clickable liver of the explorer: one block per lobe whose area is the
// share of the liver it holds, with the hepatic artery and the portal vein
// entering at the hilum and the hepatic vein leaving at the top. Geometry is
// computed from the lobe table; colours and widths come from the simulated
// values and no colour is written inline, so the light and the dark theme of
// the documentation both work (lobes.css).

const SVG = "http://www.w3.org/2000/svg";

/** An element of the diagram, with its attributes and its text. */
function node(tag, attributes = {}, text) {
  const element = document.createElementNS(SVG, tag);
  for (const [name, value] of Object.entries(attributes)) element.setAttribute(name, String(value));
  if (text !== undefined) element.textContent = text;
  return element;
}

/** The width of an edge of a flow: the same scale for every edge of the diagram. */
export function flowWidth(flow) {
  if (!Number.isFinite(flow)) return 0.7;
  return 0.7 + 1.7 * Math.sqrt(Math.abs(flow));
}

/** The columns of the lobes of a set, from the lobe table alone.
 *
 * The lobes stand in columns, the largest first, and the squarify of the
 * reference shares decides how many share a column, which keeps a block of
 * 4 % of the liver readable next to one of 30 %. The columns are the topology
 * of the drawing and are computed once per lobe table, never from the shares
 * of the moment: a lobe keeps its column and its neighbours in it, so a lobe
 * which atrophies while the one beside it grows only changes its size.
 */
function topology(lobes) {
  const cached = TOPOLOGY.get(lobes);
  if (cached) return cached;
  const entries = lobes.map(lobe => ({token: lobe.token, share: Math.max(lobe.mass_fraction ?? 0, 0)})).filter(entry => entry.share > 0);
  const total = entries.reduce((sum, entry) => sum + entry.share, 0) || 1;
  for (const entry of entries) entry.share /= total;
  const rest = entries.sort((a, b) => b.share - a.share);
  const rows = [];
  let free = {width: 1, height: 1}, remaining = 1;
  while (rest.length > 0) {
    const horizontal = true;
    const row = [rest.shift()];
    const worst = () => {
      const share = row.reduce((sum, entry) => sum + entry.share, 0);
      const length = horizontal ? free.height : free.width;
      const thickness = (share / remaining) * (horizontal ? free.width : free.height);
      return Math.max(...row.map(entry => {
        const side = (entry.share / share) * length;
        return Math.max(side / thickness, thickness / side);
      }));
    };
    while (rest.length > 0) {
      const before = worst();
      row.push(rest[0]);
      if (worst() > before) { row.pop(); break; }
      rest.shift();
    }
    const share = row.reduce((sum, entry) => sum + entry.share, 0);
    const thickness = (share / remaining) * (horizontal ? free.width : free.height);
    rows.push({horizontal, tokens: row.map(entry => entry.token)});
    remaining -= share;
    if (horizontal) free = {width: free.width - thickness, height: free.height};
    else free = {width: free.width, height: free.height - thickness};
    if (free.width <= 1e-9 || free.height <= 1e-9) break;
  }
  TOPOLOGY.set(lobes, rows);
  return rows;
}

const TOPOLOGY = new WeakMap();

/** The blocks of the lobes of a liver: the rows of its table, sized by share.
 *
 * The area of a block is the share of the lobe, so a reader sees the anatomy
 * of the lobe table and, once the remodelling moves the volumes, the atrophy
 * and the hypertrophy of a course as a change of area; the rows come from the
 * table, so a block only grows or shrinks and never moves to another side of
 * the liver. A lobe with no mass left, which is a resected one, has no block:
 * the page lists it beside the liver instead, so that the area of a block is
 * the mass of a lobe which is there and nothing else.
 */
export function blockLayout(lobes, {width, height, shares, scale = 1} = {}) {
  const area = Math.min(Math.max(Number.isFinite(scale) ? scale : 1, 0), 1);
  if (area <= 0) return [];
  // the blocks hold the area the liver has of the frame, so the whole liver of
  // the reference individual fills it and a liver which lost 70 % of its mass
  // holds 30 % of the area: a side of the laid out box is the square root of
  // the share. The box is centred across the frame and sits on its bottom
  // edge, where the hilum is, so the vessels keep their length.
  const side = Math.sqrt(area);
  const box = {width: width * side, height: height * side};
  const shift = {x: (width - box.width) / 2, y: height - box.height};
  const named = Object.fromEntries(lobes.map(lobe => [lobe.token, lobe]));
  const shareOf = token => Math.max(shares?.[token] ?? named[token]?.mass_fraction ?? 0, 0);
  const present = lobes.map(lobe => lobe.token).filter(token => shareOf(token) > 0);
  if (present.length === 0) return [];
  const total = present.reduce((sum, token) => sum + shareOf(token), 0) || 1;
  const blocks = [];
  let x = 0, y = 0, free = {width: box.width, height: box.height}, remaining = 1;
  for (const row of topology(lobes)) {
    const members = row.tokens.filter(token => shareOf(token) > 0);
    if (members.length === 0) continue;
    const horizontal = row.horizontal;
    const share = members.reduce((sum, token) => sum + shareOf(token) / total, 0);
    const thickness = (share / remaining) * (horizontal ? free.width : free.height);
    let offset = 0;
    for (const token of members) {
      const part = shareOf(token) / total;
      const length = (part / share) * (horizontal ? free.height : free.width);
      blocks.push({
        token,
        name: named[token]?.name,
        share: part,
        x: shift.x + (horizontal ? x : x + offset),
        y: shift.y + (horizontal ? y + offset : y),
        width: horizontal ? thickness : length,
        height: horizontal ? length : thickness,
      });
      offset += length;
    }
    remaining -= share;
    if (horizontal) { x += thickness; free = {width: free.width - thickness, height: free.height}; }
    else { y += thickness; free = {width: free.width, height: free.height - thickness}; }
    if (free.width <= 1e-9 || free.height <= 1e-9) break;
  }
  return blocks;
}

/** The two strokes which cut a ligated vessel, across the middle of its curve.
 *
 * The point and the tangent at the middle of the cubic the vessel is drawn
 * with, so the cut sits on the vessel whatever its bend; a vessel of no length
 * has no tangent and gets no cut.
 */
export function cutMark(p0, p1, p2, p3, {length = 11, gap = 4, at = 0.5} = {}) {
  const u = 1 - at;
  const weights = [u ** 3, 3 * u * u * at, 3 * u * at * at, at ** 3];
  const point = [0, 1].map(i => p0[i] * weights[0] + p1[i] * weights[1] + p2[i] * weights[2] + p3[i] * weights[3]);
  const mid = point;
  const tangent = [0, 1].map(i => 3 * u * u * (p1[i] - p0[i]) + 6 * u * at * (p2[i] - p1[i]) + 3 * at * at * (p3[i] - p2[i]));
  const size = Math.hypot(tangent[0], tangent[1]);
  if (!(size > 1e-9)) return [];
  const along = [tangent[0] / size, tangent[1] / size];
  const across = [-along[1], along[0]];
  return [-gap / 2, gap / 2].map(offset => {
    const centre = [mid[0] + along[0] * offset, mid[1] + along[1] * offset];
    return [
      [centre[0] - across[0] * length / 2, centre[1] - across[1] * length / 2],
      [centre[0] + across[0] * length / 2, centre[1] + across[1] * length / 2],
    ];
  });
}

const GEOMETRY = {width: 460, height: 260, margin: 16, hilum: 44, vein: 40, label: 18};
const HILUM = "hepatic artery and portal vein";
const VEIN = "hepatic vein";
// the shadings of a state: a hatch of the lobe's own colour, so the hue says
// which lobe it is and the shading says what was done to it
const SHADES = {
  portal: {angle: 45, step: 7, width: 3.2, opacity: 0.75},
  arterial: {angle: 135, step: 9, width: 2.2, opacity: 0.6},
};
let instances = 0;

/** The clickable liver of a species: the blocks of its lobes and its vessels.
 *
 * Geometry is rebuilt when the shares change (the remodelling moves them) and
 * values only update the labels and the widths otherwise.
 */
export function createLiver(host, lobes, {onPick} = {}) {
  const svg = node("svg", {
    class: "core-diagram liver",
    viewBox: `0 0 ${GEOMETRY.width} ${GEOMETRY.vein + GEOMETRY.height + GEOMETRY.hilum}`,
    role: "group",
    "aria-label": "the lobes of the liver, each a block of its share of the mass",
  });
  svg.append(node("title", {}, "The lobes of the liver"));
  svg.append(node("desc", {}, "One block per lobe, its area the share of the liver mass it holds; click a block to ligate or resect the lobe. The vessels at the bottom are the hepatic artery and the portal vein, the ones at the top the hepatic vein; the width of an edge grows with its flow and a ligated branch is cut across."));
  const uid = `liver${++instances}`;
  const defs = node("defs", {});
  const vessels = node("g", {class: "liver-vessels"});
  const blocks = node("g", {class: "liver-blocks"});
  // the outline of the whole liver of this individual: what a resection took
  // away is the area between it and the blocks
  const frame = node("rect", {
    class: "liver-frame", x: GEOMETRY.margin, y: GEOMETRY.vein + GEOMETRY.margin,
    width: GEOMETRY.width - 2 * GEOMETRY.margin, height: GEOMETRY.height - 2 * GEOMETRY.margin,
    rx: 8, opacity: "0",
  });
  const cuts = node("g", {class: "liver-cuts"});
  svg.append(defs, frame, vessels, blocks, cuts);
  host.replaceChildren(svg);

  const shapes = new Map();
  const edges = new Map();
  let drawn = "";
  let caption = null;

  // the blocks hang below the band which carries the hepatic vein
  const TOP = GEOMETRY.vein + GEOMETRY.margin;

  const draw = (shares, scale, colours) => {
    const layout = blockLayout(lobes, {
      width: GEOMETRY.width - 2 * GEOMETRY.margin,
      height: GEOMETRY.height - 2 * GEOMETRY.margin,
      shares,
      scale,
    });
    const key = layout.map(block => `${block.token}:${block.x.toFixed(2)},${block.y.toFixed(2)},${block.width.toFixed(2)},${block.height.toFixed(2)}`).join("|");
    if (key === drawn) return layout;
    drawn = key;
    shapes.clear();
    edges.clear();
    blocks.replaceChildren();
    vessels.replaceChildren();
    // one hatch per lobe and state, in the colour of that lobe: the pattern
    // cannot read the colour of the box which refers to it
    defs.replaceChildren(...lobes.flatMap(lobe => Object.entries(SHADES).map(([state, shade]) => {
      const pattern = node("pattern", {
        id: `${uid}-${lobe.token}-${state}`, width: shade.step, height: shade.step,
        patternUnits: "userSpaceOnUse", patternTransform: `rotate(${shade.angle})`,
      });
      pattern.append(node("line", {
        x1: 0, y1: 0, x2: 0, y2: shade.step,
        style: `stroke:${colours?.[lobe.token] ?? "var(--lobe-other)"};stroke-width:${shade.width};stroke-opacity:${shade.opacity}`,
      }));
      return pattern;
    })));
    const hilum = {x: GEOMETRY.width / 2, y: GEOMETRY.vein + GEOMETRY.height + GEOMETRY.hilum / 2};
    const vein = {x: GEOMETRY.width / 2, y: GEOMETRY.vein / 2};
    for (const block of layout) {
      const group = node("g", {class: "liver-lobe", tabindex: 0, role: "button", "data-lobe": block.token});
      const rect = node("rect", {
        x: block.x + GEOMETRY.margin, y: block.y + TOP,
        width: block.width, height: block.height, rx: 6, class: "liver-block",
        style: `fill:${colours?.[block.token] ?? "var(--lobe-other)"};stroke:${colours?.[block.token] ?? "var(--lobe-other)"}`,
      });
      const shade = node("rect", {
        x: block.x + GEOMETRY.margin, y: block.y + TOP,
        width: block.width, height: block.height, rx: 6, class: "liver-shade", fill: "none",
      });
      const label = node("text", {
        x: block.x + GEOMETRY.margin + block.width / 2,
        y: block.y + TOP + block.height / 2 - 2,
        class: "cd-label liver-token", "text-anchor": "middle",
      }, block.token);
      const value = node("text", {
        x: block.x + GEOMETRY.margin + block.width / 2,
        y: block.y + TOP + block.height / 2 + GEOMETRY.label - 4,
        class: "cd-sub liver-value", "text-anchor": "middle",
      }, "");
      const state = node("text", {
        x: block.x + GEOMETRY.margin + block.width / 2,
        y: block.y + TOP + block.height - 6,
        class: "cd-sub liver-state", "text-anchor": "middle",
      }, "");
      group.append(rect, shade, label, value, state);
      group.addEventListener("click", () => onPick?.(block.token));
      group.addEventListener("keydown", event => {
        if (event.key === "Enter" || event.key === " ") { event.preventDefault(); onPick?.(block.token); }
      });
      blocks.append(group);
      shapes.set(block.token, {group, rect, shade, label, value, state, block});
      const middle = block.x + GEOMETRY.margin + block.width / 2;
      const entry = {x: middle, y: TOP + GEOMETRY.height - 2 * GEOMETRY.margin + 2};
      const exit = {x: middle, y: TOP};
      // the artery and the portal branch of a lobe run side by side into the
      // bottom of the liver, so that both are visible and a closed one is read
      // off its own path; the hepatic vein leaves the top of the block
      for (const [branch, offset] of [["ha", -1], ["po", 1]]) {
        const start = hilum.x + offset * 9;
        const end = entry.x + offset * Math.min(6, block.width / 4);
        const points = [[start, hilum.y], [start, entry.y + 20], [end, entry.y + 24], [end, entry.y]];
        const path = node("path", {
          class: `liver-edge liver-edge-${branch}`,
          d: `M ${start} ${hilum.y} C ${start} ${entry.y + 20}, ${end} ${entry.y + 24}, ${end} ${entry.y}`,
        });
        vessels.append(path);
        edges.set(`${branch}_${block.token}`, {path, points});
      }
      const out = [[vein.x, vein.y], [vein.x, exit.y - 22], [exit.x, exit.y - 18], [exit.x, exit.y]];
      const outflow = node("path", {
        class: "liver-edge liver-edge-hv",
        d: `M ${vein.x} ${vein.y} C ${vein.x} ${exit.y - 22}, ${exit.x} ${exit.y - 18}, ${exit.x} ${exit.y}`,
      });
      vessels.append(outflow);
      edges.set(`hv_${block.token}`, {path: outflow, points: out});
    }
    caption = node("text", {x: hilum.x, y: hilum.y + 18, class: "cd-sub", "text-anchor": "middle"}, HILUM);
    vessels.append(caption, node("text", {x: vein.x, y: vein.y - 10, class: "cd-sub", "text-anchor": "middle"}, VEIN));
    return layout;
  };

  return {
    svg,
    /** Draw the anatomy of the shares and show the values of a sample. */
    update({shares, scale = 1, colours = {}, states = {}, values = {}, flows = {}, flowLabel = ""} = {}) {
      draw(shares, scale, colours);
      // the liver of this individual, drawn to the area it has of the whole
      frame.setAttribute("opacity", scale < 0.999 ? "1" : "0");
      for (const [token, shape] of shapes) {
        const state = states[token] ?? "intact";
        shape.group.setAttribute("data-state", state);
        shape.group.setAttribute("aria-label", `${shape.block.name}: ${STATE_LABEL[state] ?? state}`);
        shape.shade.setAttribute("fill", state in SHADES ? `url(#${uid}-${token}-${state})` : "none");
        // a small block carries the token alone: the colour and the legend
        // say what was done to it, and two labels in 30 px would collide
        const shown = values[token];
        const room = shape.block.height;
        shape.value.textContent = room < 34 ? "" : (shown ?? "");
        shape.state.textContent = room < 58 || state === "intact" ? "" : (STATE_SHORT[state] ?? state);
      }
      const cut = [];
      for (const [key, edge] of edges) {
        const flow = flows[key];
        edge.path.setAttribute("stroke-width", flowWidth(flow));
        edge.path.setAttribute("data-closed", flow !== undefined && Math.abs(flow) < 1e-12 ? "true" : "false");
        // a ligated branch is cut across, which is what was done to it; a
        // branch which merely carries nothing is dashed and is not cut
        const branch = key.slice(0, key.indexOf("_"));
        const token = key.slice(key.indexOf("_") + 1);
        const ligated = (states[token] === "portal" && branch === "po") || (states[token] === "arterial" && branch === "ha");
        if (!ligated) continue;
        // on the stretch where the branch runs alone, not in the bundle at the hilum
        for (const [from, to] of cutMark(...edge.points, {at: 0.78})) {
          cut.push(node("line", {
            class: `liver-cut liver-edge-${branch}`, "data-branch": key,
            x1: from[0], y1: from[1], x2: to[0], y2: to[1],
          }));
        }
      }
      cuts.replaceChildren(...cut);
      if (caption) caption.textContent = flowLabel ? `${HILUM}, ${flowLabel} in all` : HILUM;
    },
  };
}

const STATE_LABEL = {intact: "intact", portal: "portal vein ligated", arterial: "hepatic artery ligated", resected: "resected"};
const STATE_SHORT = {portal: "PVL", arterial: "HAL", resected: "resected"};
