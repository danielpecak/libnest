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


.. math::
	 C^0_\rho(\rho) &= \frac{3}{8} t_0 + \frac{3}{48} t_3 \rho^\alpha
	 \nonumber \\
	 C^1_\rho(\rho) &= -\frac{1}{4} t_0 \left(\frac{1}{2} + x_0\right) - \frac{1}{24} t_3 \left(\frac{1}{2} + x_3 \right) \rho^\alpha
	 \nonumber \\
	 C^0_\tau(\rho) &=  \frac{3}{16} t_1 + \frac{1}{4} t_2 \left(\frac{5}{4} + x_2\right) + \frac{3}{16} t_4 \rho^\beta + \frac{1}{4} t_5 \left(\frac{5}{4} + x_5 \right) \rho^\gamma
	 \nonumber \\
	 C^1_\tau(\rho) &= -\frac{1}{8} t_1 \left( \frac{1}{2} + x_1 \right) + \frac{1}{8} t_2 \left(\frac{1}{2} + x_2\right) - \frac{1}{8} t_4 \rho^\beta \left( \frac{1}{2} + x_4 \right)  + \frac{1}{8} t_5 \left(\frac{1}{2} + x_5 \right) \rho^\gamma

Pairing
~~~~~~~
.. math::
	\varepsilon_\pi(\rho_n,\vec\nabla\rho_n,\tilde{\rho}_n,\rho_p,\vec\nabla\rho_p,\tilde{\rho}_p)
	&=\frac{1}{4} f^\pm_n \left( v^{\pi n}(\rho_n,\rho_p) + \kappa_n|\nabla\rho_n|^2 \right) \tilde{\rho_n}^2 \nonumber \\
	&+\frac{1}{4} f^\pm_p \left( v^{\pi p}(\rho_n,\rho_p) + \kappa_p|\nabla\rho_p|^2 \right) \tilde{\rho_p}^2,

.. math::
	 \Delta_q^U(\rho_n,\rho_p) = 
	    \Delta_{SM}^U(\rho_n+\rho_p) \left( 1 - |\eta| \right)
	\pm \Delta_{NM}^U(\rho_n) \eta \frac{\rho_q}{\rho_n+\rho_p},


.. math::
	\Delta_{NM}^U(k_F< 1.38 ) = \frac{3.37968 k_F^2}{k_F^2+0.556092^2} \frac{(k_F-1.38236)^2}{(k_F-1.38236)^2+0.327517^2},

.. math::
	 \Delta_{SM}^U(k_F<1.31) =  \frac{11.5586 k_F^2}{k_F^2 + 0.489932^2}\frac{(k_F - 1.3142)^2}{(k_F - 1.3142)^2 + 0.906146^2}.




Effective mass
~~~~~~~~~~~~~~
.. math::
	B_q  = \frac{\hbar^2}{2 M^*_q} 
	 = \frac{\hbar^2}{2 M_q}  + C^\tau_0\rho  + C^\tau_1 (\rho_{q} - \rho_{q'}) 


.. todo::
   Describe properties
   
.. todo::
   Give references

.. automodule:: libnest.bsk
    :members:    


