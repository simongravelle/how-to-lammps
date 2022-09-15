Pulling a PEG molecule in water
===============================

Here a single PEG molecule is immersed in water, and forces are applied to its 
ends in order to extend it.

Creating a box of water
-----------------------

.. code-block::

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