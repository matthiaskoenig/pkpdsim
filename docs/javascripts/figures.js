// Fit the frames of the interactive figures (`iframe.pkpd-figure`) to the height
// their page posts (`pkpdsim.pbpk.visual.interactive.html_page`).
window.addEventListener("message", (event) => {
  const data = event.data;
  if (!data || data.pkpdFigure !== true || typeof data.height !== "number") return;
  for (const frame of document.querySelectorAll("iframe.pkpd-figure")) {
    if (frame.contentWindow === event.source) {
      frame.style.aspectRatio = "auto";
      frame.style.height = `${Math.ceil(data.height)}px`;
    }
  }
});
