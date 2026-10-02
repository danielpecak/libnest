Module: BSk
===========
The full reference below follows the order of the source code. Here, please find
the functions grouped by topic:




Pairing
~~~~~~~
.. image:: _static/pairing_vs_kf.png
  :width: 48 %
.. image:: _static/pairing_vs_rho.png
  :width: 48 %


* :func:`.neutron_pairing_field`
* :func:`.symmetric_pairing_field`
* :func:`.neutron_ref_pairing_field`
* :func:`.proton_ref_pairing_field`


Fermi energy, minigap, chemical potential (in :mod:`.definitions`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  * :func:`.eF_n`
  * :func:`.E_minigap_rho_n`
  * :func:`.E_minigap_delta_n`
  * :func:`.mu_q`


Statistics
~~~~~~~~~~
  * :func:`.energy_per_nucleon`
  * :func:`.pressure_n`
  * :func:`.epsilon_derivative_n`
  * :func:`.pressure_derivative_n`
  * :func:`.speed_of_sound_n`


Effective mass
~~~~~~~~~~~~~~
  * :func:`.effMn`
  * :func:`.effMp`
  * :func:`.isoscalarM`
  * :func:`.isovectorM`

Mean field
~~~~~~~~~~
  * :func:`.U_q`
  * :func:`.B_q`

Density energy functional
~~~~~~~~~~~~~~~~~~~~~~~~~
See also:
  * :func:`.epsilon_np`
  * :func:`.epsilon_rho_np`
  * :func:`.epsilon_tau_np`
  * :func:`.epsilon_delta_rho_np`
  * :func:`.epsilon_pi_np`
  * :func:`.v_pi`

Auxiliary
~~~~~~~~~
  .. image:: _static/Cs_vs_rho.png
    :width: 60 %
    :align: right

  * :func:`.C0_rho`
  * :func:`.C1_rho`
  * :func:`.C0_tau`
  * :func:`.C1_tau`
  * :func:`.I`
  * :func:`.Lambda`







.. automodule:: libnest.bsk
    :members:
