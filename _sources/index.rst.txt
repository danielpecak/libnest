libNeST Documentation
=====================
libNeST stands for Library for Neutron Star Toolkit. It contains a variety
of tools that help with handling the complicated physics of neutron stars.
The aim is to enable users to work with and learn about the physics of nuclear 
matter inside a neutron star without worrying about the complexity of Skyrme 
functionals or data formats.

.. seealso::
  The following projects are strongly interconnected:

      * `WSLDA Toolkit <https://wslda.fizyka.pw.edu.pl/>`_
      * `WBSK Toolkit <https://wbsk.fizyka.pw.edu.pl/>`_
      * `py-WDATA <https://pypi.org/project/wdata/>`_




Contributors
============

- **Monika Marek** (2023)

  Developed module :mod:`.tools`.

- **Aleksandra Bochenek** (2022)

  Developed :mod:`.bsk` by implementing and testing equations for the Brussels-Montreal family
  of density functionals.

Acknowledgments
===============
Funding
-------
- **Sonata 20 (2025-)**

  Grant no. 2024/55/D/ST2/01516

- **Sonatina 5 (2021-2024)**

  Grant no. 2021/40/C/ST2/00072 supported the development of this project.

  .. image:: _static/logo-ncn-en_large.png
     :width: 350px
     :alt: Logo of Polish National Science Center


Computational grants
--------------------
**LUMI Consortium**

We acknowledge the Polish high-performance computing infrastructure PLGrid for 
awarding access to the LUMI supercomputer, owned by the EuroHPC Joint Undertaking, 
hosted by CSC (Finland) and the LUMI consortium through grants PLL/2022/03/016433 
and PLL/2023/04/016476.

.. figure:: _static/lumi_logo_medium.png
   :width: 150px
   :alt: Logo of LUMI

.. image:: _static/plgrid_small.png
   :width: 70px
   :alt: Logo of PLGrid






.. toctree::
   :caption: Documentation
   :maxdepth: 1
   :hidden:

   instalation
   physics
   help
   bibliography

.. toctree::
   :caption: Modules
   :maxdepth: 1

   units
   definitions
   bsk
   io
   pasta
   plots
   tools

.. toctree::
   :caption: Examples
   :maxdepth: 1

   tutorial
   plotting
   inner-crust


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
