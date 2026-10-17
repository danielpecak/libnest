# Documentation Build Guide

## Local Development

### Prerequisites

From the repository root, in a virtual environment:

```bash
pip install -e ".[docs]"     # libnest + Sphinx toolchain
```

### Building Documentation

```bash
cd docs
make html                                       # quick build
make html SPHINXOPTS="-W --keep-going -n"       # strict build, exactly as CI runs it
```

Open `docs/_build/html/index.html` in your browser. Building needs network access
(intersphinx downloads the NumPy/SciPy/... inventories).

### Cleaning Build Files

```bash
make clean
```

## Figures and generated text

Every figure in the docs is generated from code at build time — no image files are
committed (only the logos in `_static/`). The same holds for tables and text that come
from the code or from repository metadata.

- `source/plots/plot_<name>.py` — one script per figure; it takes the output path as its
  only argument. Its module docstring is the figure's caption in the gallery. Shared
  helpers (`output_path`, `savefig`) and the categorical colors (`SERIES`, `LINESTYLES`:
  use them in this order, with a legend, for figures with several curves) live in
  `source/plots/_common.py`.
- `source/text/*.py` — generated rst fragments and files in `_generated/` (not committed):
  `functionals.inc` (table of the registered functionals, from `libnest.edf`),
  `gallery.inc` (every figure in `FILES` with its caption and script) and `citation.bib`
  (from `CITATION.cff`).
- `source/Makefile` — builds `_static/<name>.png` for every name in `FILES` (plus the
  favicon, `EXTRA`) and the `_generated/` files. Everything is rebuilt when its script or
  any file of the `libnest` package changes. Scripts always import the libnest of this
  checkout (the Makefile sets `PYTHONPATH`).
- `make html` runs `make -C source` first, so figures are always current.

To add a figure: create `source/plots/plot_<name>.py` with a one-sentence docstring (copy
`plot_pairing_vs_rho.py`), add `<name>.png` to `FILES` in `source/Makefile` (this also puts
it in the gallery), and reference `_static/<name>.png` from an `.rst` page.

## GitHub Pages Deployment

`.github/workflows/documentation.yml` builds the docs on every push and pull request with
`-W --keep-going -n`: **any Sphinx warning or unresolved cross-reference fails the build.**
On pushes to `main` the result is deployed to the `sphinx` branch, which GitHub Pages
serves at https://danielpecak.github.io/libnest/. A final `make linkcheck` step reports
broken external links without failing the job.

## Troubleshooting

#### A cross-reference does not resolve

Within libnest use a leading dot, e.g. :func:`.energy_per_nucleon`; for external
libraries use the full name, e.g. :func:`numpy.gradient`. Type names in docstrings must be
real types (`str`, `float`, `numpy.ndarray`), not words like `string`.

#### Missing LaTeX/Math Rendering

The project uses `sphinxcontrib-katex` for math rendering (not MathJax).

#### `make linkcheck` reports a publisher link as broken

APS, ScienceDirect and World Scientific answer automated requests with HTTP 403; they are
listed in `linkcheck_ignore` in `conf.py`. Add other such domains there.

## Dependencies

Declared in the `docs` extra of `pyproject.toml` (mirrored in `docs/requirements.txt`):
- `sphinx>=7.0` - Documentation generator
- `sphinx_rtd_theme` - ReadTheDocs theme
- `sphinxcontrib-bibtex` - Bibliography support
- `sphinxcontrib-katex` - Math equations
- `pillow` - Image processing
- `pyyaml` - Reads `CITATION.cff` for the generated BibTeX entry

## References

- [Sphinx Documentation](https://www.sphinx-doc.org/)
- [ReadTheDocs Theme](https://sphinx-rtd-theme.readthedocs.io/)
- [sphinxcontrib-bibtex](https://sphinxcontrib-bibtex.readthedocs.io/)
