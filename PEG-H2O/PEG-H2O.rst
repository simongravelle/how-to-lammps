Pulling a PEG molecule in water
===============================

Creating a box of water
-----------------------

.. code-block::units real
    atom_style full
    bond_style harmonic
    angle_style charmm
    dihedral_style charmm
    pair_style lj/cut/tip4p/long 1 2 1 1 0.1546 12.0
    kspace_style pppm/tip4p 1.0e-4

The atom_style `full` is required for charged molecules.