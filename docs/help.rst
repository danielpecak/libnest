Citing and Contributing
=======================

Citing
------
If you use libNeST in your research, please cite it. The entry below is generated from
`CITATION.cff <https://github.com/danielpecak/libnest/blob/main/CITATION.cff>`_ (GitHub
shows the same metadata as "Cite this repository"):

.. literalinclude:: _generated/citation.bib
   :language: bibtex

Please also cite the papers behind the parts of the functional you use:

* BSk31: :cite:`chamel2016further`; BSkG3: :cite:`grams2023skyrme`; BSkG4:
  :cite:`grams2025skyrme` (the reference of every registered functional is listed in the
  table on the page :doc:`functionals`);
* the reference pairing gaps: :cite:`cao2006screening`, and the interpolation scheme to
  asymmetric matter: :cite:`chamel2009pairing` (BSk), :cite:`grams2023skyrme` (BSkG3),
  :cite:`grams2025skyrme` (BSkG4) — see :doc:`pairing`;
* the pairing strength obtained from the gaps: :cite:`chamel2010effective`.

Contributing
------------
Bug reports and pull requests are welcome at
`github.com/danielpecak/libnest <https://github.com/danielpecak/libnest>`_.

1. Fork the repository and create a feature branch.
2. Install in editable mode with the test tools: ``pip install -e ".[test]"``.
3. Run the test suite with ``pytest`` and make sure it passes.
4. Open a pull request.
