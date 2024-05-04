#######
Physics
#######
Here you can find a brief introduction to the physics of neutron stars.

.. contents:: Table of Contents
   :depth: 10

*************************
Structure of Neutron Star
*************************

Crust
=====
In the crust we can show certain densities that are relevant. The major division
is to **outer crust** and **inner crust**.

Pairing S
---------
.. image:: _static/pairing_vs_kf.png
  :width: 48 %
.. image:: _static/pairing_vs_rho.png
  :width: 48 %

Pairing P
---------
Give some references.

Pasta phase
-----------

******************
Density Functional
******************

Skyrme force
============
The first paper about vortex pinning within fully dynamical approach in neutron
matter was :cite:`wlazlowski2016vortex`. A standard SLy4 type of functional was
used (see :cite:`bulgac2016induced`) of the following form:

.. math::

  v_{12}^{(2)}
  = \phantom{\frac{1}{2}}t_0 (1+x_0 \hat P_\sigma) \delta(\bm r)

  + \frac{1}{2} t_1 (1+x_1 \hat P_\sigma) (\bm k'^2\delta(\bm r) + \delta(\bm r) \bm k^2)

  + \phantom{\frac{1}{2}}t_2 (1+x_2 \hat P_\sigma) \bm k' \cdot \delta(\bm r) \bm k

  + \frac{1}{6} t_3 (1+x_3 \hat P_\sigma) \rho^\gamma(\bm r) \delta(\bm r)

  + i W_0 (\sigma_1 + \sigma_2) \cdot (\bm k'\times \delta(\bm r) \bm k),

where the consecutive terms have a physical meaning and interpretation (see Sec. 2
in :cite:`chabanat1998skyrme`):

 * :math:`t_0`: contact (zero range) interaction
 * :math:`t_1, t_2`: momentum exchange (finite range term) and effective mass term (which is proportional to the kinetic density :math:`\tau`)
 * :math:`t_3`: for :math:`\gamma \rightarrow 1` it is the three body term averaged over the third particle; in general a density-dependent term
 * :math:`W_0`: spin-orbit coupling


It turns out, that similarly as for the term :math:`t_3` where in general the coefficient
:math:`\gamma\neq1`, the terms :math:`t_1,t_2` can be generalized to density
dependent terms :math:`t_4,t_5`, respectively :cite:`chamel2009further`. Increasing
the number of parameters leads to the better fit of very large amount of data.


Why BSk functional?
===================
BSk density functional has been designed to describe accurately nuclear matter
in neutron stars. Hence, the fitting was done not only using experimental nuclear
data (as it is usually done) but also theoretical results for many-body calculations
(for example pairing gaps :cite:`cao2006screening`).

Experimental data
 * all atomic masses with :math:`Z,N \ge 8` from the Atomic Mass Evaluation (2353)
 * nuclear charge radii
 * symmetry energy :math:`29 \mathrm{MeV} \le J \le 32 \mathrm{MeV}`
 * incompressibility :math:`K_V = (240 \pm 10) \mathrm{MeV}`

N-body calculations using realistic forces
 * equation of state of pure neutron matter
 * :math:`{}^1S_0` pairing gaps in nuclear matter
 * effective masses in nuclear matter
 * stability against spin and spin-isospin fluctuations



History of development of BSk
=============================
There is a wide variety of different BSk functionals. Below the history of updates
and improvements of the model is shown:

 * fit to realistic :math:`{}^1S_0` pairing gaps (no self-energy)  (BSk16--17)
   :cite:`chamel2008further,chamel2009pairing,chamel2010effective`
 * removal of spurious spin-isospin instabilities :cite:`chamel2009further,chamel2010spin` (BSk18)
 * fit to realistic neutron-matter equations of state :cite:`goriely2010further` (BSk19--21)
 * fit to different symmetry energies :cite:`goriely2013further` (BSk22--26)
 * optimal fit of the 2012 AME :cite:`goriely2013hartree` (BSk27*)
 * genealized spin-orbit coupling :cite:`goriely2015further`  (BSk28--29)
 * fit to realistic :math:`{}^1S_0` pairing gaps with self-energy :cite:`chamel2016further` (BSk30--32)



BSk functional
==============
Based on :cite:`chamel2009further,chamel2016further` the Skyrme-type force reads:

.. math::

    v_{12}^{(2)}

    = \phantom{\frac{1}{2}}t_0 (1+x_0 \hat P_\sigma) \delta(\bm r)

    + \frac{1}{2} t_1 (1+x_1 \hat P_\sigma) (\bm k'^2\delta(\bm r) + \delta(\bm r) \bm k^2)

    + \phantom{\frac{1}{2}}t_2 (1+x_2 \hat P_\sigma) \bm k' \cdot \delta(\bm r) \bm k

    + \frac{1}{6} t_3 (1+x_3 \hat P_\sigma) \rho^\alpha(\bm r) \delta(\bm r)

    + \frac{1}{2} t_4 (1+x_4 \hat P_\sigma) \left[ \bm k'^2 \rho^\beta(\bm r) \delta(\bm r) + \delta(\bm r) \rho^\beta(\bm r) \bm k^2 \right]

    + \phantom{\frac{1}{2}}t_5 (1+x_5 \hat P_\sigma) \bm k' \rho^\gamma(\bm r) \delta(\bm r) \bm k

    + i W_0 (\sigma_1 + \sigma_2) \cdot (\bm k'\times \delta(\bm r) \bm k)

    + \sum_q f^\pm_q  \left(  v^{\pi,q}[\rho_n(\bm r),\rho_p(\bm r)]  + \kappa_q |\nabla\rho(\bm r)|^2  \right) \delta(\bm r).

The core of this BSk functional is the Skyrme functional showed above. Here, in
addition to :math:`t_1,t_2` terms, there are similar expressions that depend
additionally on the power of density :math:`\rho`.


Energy
------
.. image:: _static/energySM.png
  :width: 48 %
.. image:: _static/energySMZoom.png
  :width: 48 %


Pure neutron matter
^^^^^^^^^^^^^^^^^^^
Neutron stars are built from pure neutron matter - this is very good approximation.


Symmetric matter
^^^^^^^^^^^^^^^^
Atomic nuclei (especially those of lower mass) have similar number protons and neutrons.
Therefore, it's convenient to consider symmetric matter. We know that
for larger number of nucleons, the system tends to be more neutron-rich
to remove Coulomb interaction. In contrast to pure neutron matter, there is no
symmetric nuclear matter in the nature. However this concept is very useful
for theoretical considerations.


Pairing models
==============
References for pairing models:

 * :cite:`chamel2009pairing`
 * :cite:`chamel2010effective`
 * :cite:`chamel2013pairing`
 * :cite:`chamel2016further`


.. todo::
  Expand this section: Summarize and fill in pairing models, physical meaning
  and papers for reference.





.. bibliography::
   :filter: docname in docnames
