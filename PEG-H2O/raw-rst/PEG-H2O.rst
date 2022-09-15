Pulling a PEG molecule in water
===============================

Here a single PEG molecule is immersed in water, and forces are applied to its 
ends in order to extend it.

Creating a box of water
-----------------------

.. raw:: html

   <div>Some stuff <pre>some <b>bold</b> text</pre>...</div>

   <p style="color:red;">I am red</p>

.. code-block:: python

    # LAMMPS input script
    units real
    atom_style full
    bond_style harmonic
    angle_style charmm
    dihedral_style charmm
    pair_style lj/cut/tip4p/long 1 2 1 1 0.1546 12.0
    kspace_style pppm/tip4p 1.0e-4

Note: The atom_style `full` is required for charged molecules, and the pair_style `lj/cut/tip4p/long`
is a Lennard-Jones (cut) - Coulomb (long) pair style specifically adapted to 4 points water models,
which is what we want here. 

.. code-block::

    region box block -40 40 -15 15 -15 15
    create_box 7 box &
    bond/types 6 &
    angle/types 9 &
    dihedral/types 14 &
    extra/bond/per/atom 2 &
    extra/angle/per/atom 1 &
    extra/special/per/atom 2

.. code-block::

    include ../PARM.lammps

.. code-block::

    molecule h2omol H2OTip4p.txt
    create_atoms 0 random 500 456415 NULL mol h2omol 454756
    group H2O type 1 2
    delete_atoms overlap 2 H2O H2O mol yes

.. code-block::

    fix myshk H2O shake 1.0e-4 200 0 b 1 a 1 mol h2omol
    fix mynpt all npt temp 300 300 100 iso 1 1 1000

.. code-block::

    dump mydmp all atom 100 dump.lammpstrj
    variable mytemp equal temp
    variable myvol equal vol
    fix myat1 all ave/time 10 10 100 v_mytemp file temperature.dat
    fix myat2 all ave/time 10 10 100 v_myvol file volume.dat

.. code-block::

    timestep 1.0
    thermo 1000
    run 20000

.. code-block::

    write_data H2O.data
