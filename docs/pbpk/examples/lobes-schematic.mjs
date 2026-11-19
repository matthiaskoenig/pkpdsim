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

/** The blocks of the lobes of a liver: a squarified treemap of their shares.
 *
 * The area of a block is the share of the lobe, so a reader sees the anatomy
 * of the lobe table and, once the remodelling moves the volumes, the atrophy
 * and the hypertrophy of a course as a change of area. A lobe with no mass
 * left, which is a resected one, has no block: the page lists it beside the
 * liver instead, so that the area of a block is the mass of a lobe which is
 * there and nothing else.
 */
export function blockLayout(lobes, {width, height, shares} = {}) {
  const entries = lobes
    .map(lobe => ({
      token: lobe.token,
      name: lobe.name,
      share: Math.max(shares?.[lobe.token] ?? lobe.mass_fraction, 0),
    }))
    .filter(entry => entry.share > 0);
  if (entries.length === 0) return [];
  const total = entries.reduce((sum, entry) => sum + entry.share, 0) || 1;
  for (const entry of entries) entry.share /= total;
  const blocks = [];
  // the rows of a treemap: the largest lobes first, each row as square as it
  // can be, which keeps a block of 4 % of the liver readable next to one of 30 %
  const rest = [...entries].sort((a, b) => b.share - a.share);
  let x = 0, y = 0, free = {width, height}, remaining = 1;
  while (rest.length > 0) {
    const horizontal = free.width >= free.height;
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
    let offset = 0;
    for (const entry of row) {
      const side = (entry.share / share) * (horizontal ? free.height : free.width);
      blocks.push({
        token: entry.token,
        name: entry.name,
        share: entry.share,
        x: horizontal ? x : x + offset,
        y: horizontal ? y + offset : y,
        width: horizontal ? thickness : side,
        height: horizontal ? side : thickness,
      });
      offset += side;
    }
    remaining -= share;
    if (horizontal) { x += thickness; free = {width: free.width - thickness, height: free.height}; }
    else { y += thickness; free = {width: free.width, height: free.height - thickness}; }
    if (free.width <= 1e-9 || free.height <= 1e-9) break;
  }
  return blocks;
}

const GEOMETRY = {width: 460, height: 260, margin: 16, hilum: 44, label: 18};
const HILUM = "hepatic artery and portal vein";

/** The clickable liver of a species: the blocks of its lobes and its vessels.
 *
 * Geometry is rebuilt when the shares change (the remodelling moves them) and
 * values only update the labels and the widths otherwise.
 */
export function createLiver(host, lobes, {onPick} = {}) {
  const svg = node("svg", {
    class: "core-diagram liver",
    viewBox: `0 0 ${GEOMETRY.width} ${GEOMETRY.height + GEOMETRY.hilum}`,
    role: "group",
    "aria-label": "the lobes of the liver, each a block of its share of the mass",
  });
  svg.append(node("title", {}, "The lobes of the liver"));
  svg.append(node("desc", {}, "One block per lobe, its area the share of the liver mass it holds; click a block to ligate or resect the lobe. The vessels at the bottom are the hepatic artery and the portal vein, the one at the top the hepatic vein; the width of an edge grows with its flow."));
  const vessels = node("g", {class: "liver-vessels"});
  const blocks = node("g", {class: "liver-blocks"});
  svg.append(vessels, blocks);
  host.replaceChildren(svg);

  const shapes = new Map();
  const edges = new Map();
  let drawn = "";
  let caption = null;

  const draw = shares => {
    const layout = blockLayout(lobes, {
      width: GEOMETRY.width - 2 * GEOMETRY.margin,
      height: GEOMETRY.height - 2 * GEOMETRY.margin,
      shares,
    });
    const key = layout.map(block => `${block.token}:${block.x.toFixed(2)},${block.y.toFixed(2)},${block.width.toFixed(2)},${block.height.toFixed(2)}`).join("|");
    if (key === drawn) return layout;
    drawn = key;
    shapes.clear();
    edges.clear();
    blocks.replaceChildren();
    vessels.replaceChildren();
    const hilum = {x: GEOMETRY.width / 2, y: GEOMETRY.height + GEOMETRY.hilum / 2};
    for (const block of layout) {
      const group = node("g", {class: "liver-lobe", tabindex: 0, role: "button", "data-lobe": block.token});
      const rect = node("rect", {
        x: block.x + GEOMETRY.margin, y: block.y + GEOMETRY.margin,
        width: block.width, height: block.height, rx: 6, class: "liver-block",
      });
      const label = node("text", {
        x: block.x + GEOMETRY.margin + block.width / 2,
        y: block.y + GEOMETRY.margin + block.height / 2 - 2,
        class: "cd-label liver-token", "text-anchor": "middle",
      }, block.token);
      const value = node("text", {
        x: block.x + GEOMETRY.margin + block.width / 2,
        y: block.y + GEOMETRY.margin + block.height / 2 + GEOMETRY.label - 4,
        class: "cd-sub liver-value", "text-anchor": "middle",
      }, "");
      const state = node("text", {
        x: block.x + GEOMETRY.margin + block.width / 2,
        y: block.y + GEOMETRY.margin + block.height - 6,
        class: "cd-sub liver-state", "text-anchor": "middle",
      }, "");
      group.append(rect, label, value, state);
      group.addEventListener("click", () => onPick?.(block.token));
      group.addEventListener("keydown", event => {
        if (event.key === "Enter" || event.key === " ") { event.preventDefault(); onPick?.(block.token); }
      });
      blocks.append(group);
      shapes.set(block.token, {group, rect, label, value, state, block});
      const entry = {x: block.x + GEOMETRY.margin + block.width / 2, y: GEOMETRY.height - GEOMETRY.margin + 2};
      // the artery and the portal branch of a lobe run side by side, so that
      // both are visible and a closed one is read off its own path
      for (const [branch, offset] of [["ha", -1], ["po", 1]]) {
        const start = hilum.x + offset * 9;
        const end = entry.x + offset * Math.min(6, block.width / 4);
        const path = node("path", {
          class: `liver-edge liver-edge-${branch}`,
          d: `M ${start} ${hilum.y} C ${start} ${entry.y + 20}, ${end} ${entry.y + 24}, ${end} ${entry.y}`,
        });
        vessels.append(path);
        edges.set(`${branch}_${block.token}`, path);
      }
    }
    caption = node("text", {x: hilum.x, y: hilum.y + 18, class: "cd-sub", "text-anchor": "middle"}, HILUM);
    vessels.append(caption);
    return layout;
  };

  return {
    svg,
    /** Draw the anatomy of the shares and show the values of a sample. */
    update({shares, states = {}, values = {}, flows = {}, flowLabel = ""} = {}) {
      draw(shares);
      for (const [token, shape] of shapes) {
        const state = states[token] ?? "intact";
        shape.group.setAttribute("data-state", state);
        shape.group.setAttribute("aria-label", `${shape.block.name}: ${STATE_LABEL[state] ?? state}`);
        // a small block carries the token alone: the colour and the legend
        // say what was done to it, and two labels in 30 px would collide
        const shown = values[token];
        const room = shape.block.height;
        shape.value.textContent = room < 34 ? "" : (shown ?? "");
        shape.state.textContent = room < 58 || state === "intact" ? "" : (STATE_SHORT[state] ?? state);
      }
      for (const [key, path] of edges) {
        const flow = flows[key];
        path.setAttribute("stroke-width", flowWidth(flow));
        path.setAttribute("data-closed", flow !== undefined && Math.abs(flow) < 1e-12 ? "true" : "false");
      }
      if (caption) caption.textContent = flowLabel ? `${HILUM}, ${flowLabel} in all` : HILUM;
    },
  };
}

const STATE_LABEL = {intact: "intact", portal: "portal vein ligated", arterial: "hepatic artery ligated", resected: "resected"};
const STATE_SHORT = {portal: "PVL", arterial: "HAL", resected: "resected"};
