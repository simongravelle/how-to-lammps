How to measure surface heat of adsorption
=========================================

Here a thin layer of water is in contact with a NaCl(100) solid surface.
Grand canonical Monte Carlo simulation is performed to equilibrate the
content in water.

This simulation was taken from my `article`_ with Alex and Christian.

.. _video: https://youtu.be/05DgPNfjReY
.. _article: https://doi.org/10.1063/5.0099646

The input file is described here, and is available as well, together with all
required data and parameter files in this `folder`_.

.. _folder: files/


.. code-block:: python

     variable rdm equal 544645
     variable pre equal 1
     variable tem equal 297.15
     variable dis index 0.5
     variable tfac equal 5.0/3.0
     variable mu equal -12.275


The value of the variable mu, which is the imposed chemical potential, can be calibrated
by performing gcmc simulation in absence of solid (just vapor water). The calibration
consists in performing simulations for varying value of mu, and measuring the
equilibrated pressure and water density.

.. code-block:: python

     units real
     atom_style full
     bond_style harmonic
     angle_style harmonic
     boundary p p p
     pair_style lj/cut/tip4p/long 1 2 1 1 0.105 10.0
     kspace_style pppm/tip4p 1.0e-4
     pair_modify mix arithmetic tail yes


Here a four-point water models is used (tip4p/epsilon). tip4p/epsilon is
generally a good choice as it reproduces very accurately many properties
of water under ambient conditions.

.. code-block:: python

     read_data gcmc.data
     include PARM.lammps


.. code-block:: python

     group H2O type 1 2
     group NaCl type 3 4


.. code-block:: python

     molecule h2omol H2OTip4p.txt
     fix myshk H2O shake 1.0e-4 200 0 b 1 a 1 mol h2omol
     fix mynvt all nvt temp ${tem} ${tem} 100
     fix tether NaCl spring/self 10.0
     fix mygcmc H2O gcmc 100 50 0 0 ${rdm} ${tem} ${mu} ${dis} mol h2omol tfac_insert ${tfac} group H2O shake myshk full_energy
     timestep 2.0

Here the NaCl block is maintained into position with spring/self, in order to
avoid a vertical drift of the whole system.

The first number and second numbers in the gcmc fix (respectively 100 and 50) are
important, they set the ratio between GCMC attempts and molecular dynamics (MD) moves. Here every
100 steps, 50 GCMC moves are attempted, and in between those 100 step, MD is performed.


.. code-block:: python

     dump mydmp all atom 100 dump.lammpstrj
     thermo 100
     variable oxygen atom "type==1"
     group oxygen dynamic all var oxygen
     variable nO equal count(oxygen)
     variable myetot equal etotal
     fix at1 all ave/time 100 10 1000 v_nO file numbermolecule.dat
     fix at2 all ave/time 100 10 1000 v_myetot file totalenergy.dat
     variable ins_att equal f_mygcmc[3]
     variable ins_suc equal f_mygcmc[4]
     variable del_att equal f_mygcmc[5]
     variable del_suc equal f_mygcmc[6]
     fix at3 all ave/time 100 10 1000 v_ins_att v_ins_suc v_del_att v_del_suc file gcmc.dat
     run 100000
     write_data gcmc.data

The system total energy and number of molecule are needed for the calculation of the
adsorption heat.
