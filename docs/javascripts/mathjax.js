// The equations of the documentation: the generic mode of pymdownx.arithmatex
// writes them as `<script type="math/tex">`, which MathJax 3 typesets. The
// configuration follows the one of Material for MkDocs; `subscriptionCenter`
// of the instant navigation re-typesets a page which was swapped in.
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex",
  },
};

document.addEventListener("DOMContentLoaded", () => {
  if (window.MathJax && window.MathJax.typesetPromise) {
    window.MathJax.typesetPromise();
  }
});
