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

     variable rdm equal 544
     variable pre equal 1
     variable tem equal 297.15
     variable dis index 0.5
     variable tfac equal 5.0/3.0
     variable mu equal -12.60


.. code-block:: python

     units real
     atom_style full
     bond_style harmonic
     angle_style harmonic
     boundary p p p
     pair_style lj/cut/tip4p/long 1 2 1 1 0.105 10.0
     kspace_style pppm/tip4p 1.0e-4
     pair_modify mix arithmetic tail yes


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
     fix mygcmc H2O gcmc 50 100 0 0 ${rdm} ${tem} ${mu} ${dis} mol h2omol tfac_insert ${tfac} group H2O shake myshk full_energy
     timestep 2.0


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
     run 50000
     write_data gcmc.data
