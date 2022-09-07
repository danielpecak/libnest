Module:BSk
==========
Formulas:



.. math::
	\varepsilon_{\Delta\rho}&(\rho_n,\vec\nabla\rho_n,\rho_p,\vec\nabla\rho_p)	= \\
	 &+\frac{3}{16} t_1 \left[
	  \left(          1 + \frac{1}{2}x_1\right)         \left( \nabla(\rho_n+\rho_p)   \right)^2
	 -\left(\frac{1}{2} +           x_1 \right)  \sum_q \left( \nabla\rho_q \right)^2
	 \right]  \\
	  &- \frac{1}{16}t_2 \left[
	     \left(           1 + \frac{1}{2} x_2\right)          \left( \nabla(\rho_n+\rho_p)  \right)^2
	   + \left( \frac{1}{2} +             x_2\right) \sum_{q} \left( \nabla\rho_q\right)^2  \right]  \\
	  &+\frac{3}{16} t_4  (\rho_n+\rho_p)^\beta \left[
	  \left(1 + \frac{1}{2}x_4\right)           \left( \nabla(\rho_n+\rho_p)   \right)^2
	 - \left(\frac{1}{2} + x_4 \right)  \sum_q  \left( \nabla\rho_q \right)^2
	 \right]  \\
	 &+ \frac{\beta}{8}t_4 (\rho_n+\rho_p)^{\beta-1} \left[
	 \left(1 + \frac{1}{2}x_4\right) (\rho_n+\rho_p)\left(\nabla(\rho_n+\rho_p)\right)^2
	 - \left(\frac{1}{2} + x_4 \right)
	 \nabla(\rho_n+\rho_p)\cdot\left(\sum_q \rho_q\nabla\rho_q\right)
	 \right]  \\
	 &- \frac{1}{16}t_5 (\rho_n+\rho_p)^\gamma \left[
	   \left(          1 + \frac{1}{2} x_5\right) \left( \nabla(\rho_n+\rho_p)\right)^2
	  +\left(\frac{1}{2} +             x_5\right) \sum_{q} \left( \nabla\rho_q\right)^2
	 \right].

.. math::
	\varepsilon_\rho(\rho_n,\rho_p)
	&= C^\rho_0(\rho_n+\rho_p) \left( \rho_n + \rho_p \right)^2
	 + C^\rho_1(\rho_n+\rho_p) \left( \rho_n - \rho_p \right)^2
	 \\
	\varepsilon_\tau(\rho_n,\tau_n,{\bm j}_n,\rho_p,\tau_p,{\bm j}_p)
	&= C^\tau_0(\rho_n+\rho_p) \left[ (\rho_n + \rho_p)(\tau_n + \tau_p) - ({\bm j}_n + {\bm j}_p)^2 \right] \nonumber \\
	&+ C^\tau_1(\rho_n+\rho_p) \left[ (\rho_n - \rho_p)(\tau_n - \tau_p) - ({\bm j}_n - {\bm j}_p)^2 \right].


Pairing
~~~~~~~





Effective mass
~~~~~~~~~~~~~~


.. todo::
   Describe properties

.. todo::
   Give references

.. automodule:: libnest.bsk
    :members:

References
----------
.. bibliography::
