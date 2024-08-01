/* custom.js */

// Function to force MathJax to re-render the page
function renderMathJax() {
  if (typeof MathJax !== 'undefined') {
    MathJax.typesetPromise();
  }
}

// Ensure MathJax rendering is triggered at the beginning
document.addEventListener("DOMContentLoaded", renderMathJax);
