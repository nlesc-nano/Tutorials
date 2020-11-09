.. _fitting:

Forcefield Optimization
=======================
The goal here is to fit a classical potential energy surface (PES) to an
*ab-initio* PES by optimizing the classical forcefield parameters.
This forcefield optimization is conducted using the Addaptive Rate Monte
Carlo (ARMC) method described by S. Cosseddu *et al* in
*J. Chem. Theory Comput.*, **2017**, *13*, 297–308.

You need to start from an *ab-initio* Molecular Dynamics (MD) trajectory.
In this tutorial, we will start from the DFT .... 

Then, we are going to build a .yaml file with the ARMC settings containing the following sections.

The param block
---------------
The :attr:`param<param.block>` key in the .yaml input contains all user-specified
to-be optimized parameters, i.e. the charge (Coulomb potential), and the epsilon and sigma parameters (Lennard-Jones potential).

There are three critical (and two optional) components to the ``"param"`` block:

    * The key of each block (charge_, epsilon_ & sigma_).
    * The ``"keys"`` sub-block, which points to the section path in the CP2K settings (*e.g.* `['input', 'force_eval', 'mm', 'forcefield', 'charge'] <https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/CHARGE.html>`_).
    * The sub-blocks containing either singular atoms_ or `atom pairs <https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_ATOMS>`_.

Besides the three above-mentioned mandatory components, one can (optionally) supply the unit of the parameter and/or constrain its value to a certain range. When supplying units, it is the responsibility of the user to ensure the units are supported by CP2K. Furthermore, parameter constraints are, as of the moment, limited to specifying minimum and/or maximum values (e.g. 0 < Cs < 2). Aditionally, one can add a "frozen" component, where the specified parameters will be treated as constants rather than to-be optimized variables.

Here is the ``"param"`` block for the fitting of our double perovskite nanocrystal capped with acetate ligands.

.. code:: yaml

    param:
        charge:
            param: charge
            keys: [input, force_eval, mm, forcefield, charge]
            Cs: 1.0
            Ag: 1.0
            In: 3.0
            Cl: -1.0
            C2O3:  0.247973
            O2D2: -0.5739865
            constraints:
                - '0 < Cs < 1.5'
                - '0.2 < Ag < 1.7'
                - '1.2 < In < 3.5'
                - '-1.5 < Cl < 0'
        lennard_jones:
            - param: epsilon
              unit: kjmol
              Cs Cs: 0.1882
              Cs Pb: 0.7227
              Pb Pb: 2.7740
            - unit: nm
              param: sigma
              constraints: 'Cs Cs == Pb Pb'
              Cs Cs: 0.60
              Cs Pb: 0.50
              Pb Pb: 0.60


