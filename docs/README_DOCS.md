# Documentation Build Guide

## Local Development

### Prerequisites

```bash
# Install documentation dependencies
pip install -r requirements.txt
```

### Building Documentation

```bash
# Build HTML documentation
make html

# View documentation
# Open docs/_build/html/index.html in your browser
```

### Cleaning Build Files

```bash
make clean
```

## GitHub Pages Deployment

Documentation is automatically built and deployed to GitHub Pages when changes are pushed to the `main` branch.

### Workflow Details

- **File**: `.github/workflows/documentation.yml`
- **Python Version**: 3.11 (to avoid `imghdr` compatibility issues with Python 3.13)
- **Sphinx Version**: ≥7.0 (modern version with better Python 3.11+ support)
- **Deploy Branch**: `sphinx`

### Troubleshooting

#### Error: "No module named 'imghdr'"

This error occurs with:
- Python 3.13+ (where `imghdr` was removed)
- Older Sphinx versions (<7.0)

**Solution**: The workflow now uses Python 3.11 and Sphinx ≥7.0

#### Missing LaTeX/Math Rendering

The project uses `sphinxcontrib-katex` for math rendering (not MathJax).
Ensure it's installed in requirements.

#### Build Warnings

Common warnings can be suppressed in `conf.py`:
```python
suppress_warnings = ['epub.unknown_project_files']
```

## Dependencies

See `requirements.txt` for full list:
- `sphinx>=7.0` - Documentation generator
- `sphinx_rtd_theme` - ReadTheDocs theme
- `sphinxcontrib-bibtex` - Bibliography support
- `sphinxcontrib-katex` - Math equations
- `pillow` - Image processing

## References

- [Sphinx Documentation](https://www.sphinx-doc.org/)
- [ReadTheDocs Theme](https://sphinx-rtd-theme.readthedocs.io/)
- [sphinxcontrib-bibtex](https://sphinxcontrib-bibtex.readthedocs.io/)
