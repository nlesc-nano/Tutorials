.. _fitting:

Forcefield Optimization
=======================
The goal here is to obtain the classical forcefield parameters for a double perovksite nanocrystal (NC) capped by carboxylate ligands.
To do that, we will start from the *ab-initio* Molecular Dynamics (MD) trajectory of a 3nm sided cubic Cs2AgInCl6 NC capped by 50% of acetate molecules (see tutorial for Quantum Dot construction).

Before starting the fitting we invite you to read the documentation relative to Addaptive Rate Monte Carlo (ARMC) in Auto-FOX.

First, let's build the .yaml file with our ARMC settings.
The file is composed of the following sections.

The param block
---------------
.. code:: yaml

    param:
        charge:
            param: charge
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
              frozen:
                  guess: uff
            - param: sigma
              unit: nm
              Cs Cs:  0.385
              Cs Ag:  0.310
              Cs In:  0.370
              Cs Cl:  0.290
              Ag Ag:  0.635
              Ag In:  0.390
              Ag Cl:  0.215
              In In:  0.685
              In Cl:  0.215
              Cl Cl:  0.295
              C2O3 Cs: 0.295
              C2O3 Ag: 0.255
              C2O3 In: 0.240
              C2O3 Cl: 0.300
              O2D2 Cs: 0.245
              O2D2 Ag: 0.190
              O2D2 In: 0.195
              O2D2 Cl: 0.275
              constraints:
                  - 'Cs Cs   > 0.365'
                  - 'Cs Ag   > 0.290'
                  - 'Cs In   > 0.350'
                  - 'Cs Cl   > 0.270'
                  - 'Ag Ag   > 0.615'
                  - 'Ag In   > 0.370'
                  - 'Ag Cl   > 0.195'
                  - 'In In   > 0.665'
                  - 'In Cl   > 0.195'
                  - 'Cl Cl   > 0.275'
                  - 'C2O3 Cs > 0.275'
                  - 'C2O3 Ag > 0.235'
                  - 'C2O3 In > 0.220'
                  - 'C2O3 Cl > 0.280'
                  - 'O2D2 Cs > 0.225'
                  - 'O2D2 Ag > 0.170'
                  - 'O2D2 In > 0.175'
                  - 'O2D2 Cl > 0.255'
              frozen:
                 C331 Cs: 0.295
                 C331 Ag: 0.255
                 C331 In: 0.240
                 C331 Cl: 0.300
                 HGA3 Cs: 0.245
                 HGA3 Ag: 0.215
                 HGA3 In: 0.260
                 HGA3 Cl: 0.210

The ``"param"`` key contains all user-specified features for the to-be optimized parameters for the Coulomb potential (the charge_) and the Lennard-Jones potential
(epsilon_ & sigma_).
Let's have a look at the relative sub-blocks:

1.  Coulomb potential

    .. code:: yaml

            param: charge
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

    The initial parameters for the charges are simply:
    * for the nanocrystal core ions (Cs, Ag, In, Cl), the most stable oxidation state;
    * for the anchoring group of the ligand (COO group of the acetate, i.e. C2O3 and O2D2), the charges are choosen in order to have an overall charge neutral system.
    In this case, the core ions charges are constrained to a certain range in order to keep the correct oxidation state (for example In is constrained to values higher
    than 1 to keep the oxidation number +3 it has in a double perovskite). However the constraints component is optional.

Let's move to the :code:`lennard_jones` block.

2.  Lennard-Jones potential
    This sub-block is divided itself in two components.

    .. code:: yaml

            - param: epsilon
              unit: kjmol
              frozen:
                  guess: uff
    In our fitting the epsilon parameters treated as constants rather than to-be optimized variables (all frozen) and all the values are guessed using the `uff <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo.html#parameter-guessing>`_ procedure. Specifying the epsilon parameters (even without optimizing them) helps in a more accurate
    fitting procedure.

    .. code:: yaml

            - param: sigma
              unit: nm
              Cs Cs:  0.385
              Cs Ag:  0.310
              Cs In:  0.370
              Cs Cl:  0.290
              Ag Ag:  0.635
              Ag In:  0.390
              Ag Cl:  0.215
              In In:  0.685
              In Cl:  0.215
              Cl Cl:  0.295
              C2O3 Cs: 0.295
              C2O3 Ag: 0.255
              C2O3 In: 0.240
              C2O3 Cl: 0.300
              O2D2 Cs: 0.245
              O2D2 Ag: 0.190
              O2D2 In: 0.195
              O2D2 Cl: 0.275
              constraints:
                  - 'Cs Cs   > 0.365'
                  - 'Cs Ag   > 0.290'
                  - 'Cs In   > 0.350'
                  - 'Cs Cl   > 0.270'
                  - 'Ag Ag   > 0.615'
                  - 'Ag In   > 0.370'
                  - 'Ag Cl   > 0.195'
                  - 'In In   > 0.665'
                  - 'In Cl   > 0.195'
                  - 'Cl Cl   > 0.275'
                  - 'C2O3 Cs > 0.275'
                  - 'C2O3 Ag > 0.235'
                  - 'C2O3 In > 0.220'
                  - 'C2O3 Cl > 0.280'
                  - 'O2D2 Cs > 0.225'
                  - 'O2D2 Ag > 0.170'
                  - 'O2D2 In > 0.175'
                  - 'O2D2 Cl > 0.255'
              frozen:
                 C331 Cs: 0.295
                 C331 Ag: 0.255
                 C331 In: 0.240
                 C331 Cl: 0.300
                 HGA3 Cs: 0.245
                 HGA3 Ag: 0.215
                 HGA3 In: 0.260
                 HGA3 Cl: 0.210

    Here we need to optimize the sigma parameters for the all pair interactions of interest (provided with the corresponding `atom pairs <https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_ATOMS>`_): 
    the ion-ion interactions inside the nanocrystal core (eg. Cs-Cs) and the acetate anchoring group-core ions interactions (eg. O2D2-Cs).
    The initial parameters for these pairs are obtained from the DFT trajectory by mean of a small python script:

    .. code:: python

        >>> import pandas as pd
        >>> from FOX import MultiMolecule, example_xyz, estimate_lj

        >>> xyz_file: str = 'qmworks-cp2k-pos-1_RUN.xyz' # path of DFT trajectory
        >>> atom_subset = ['Cs', 'Ag', 'In', 'Cl', 'C', 'O'] # core ions and acetate anchoring group

        >>> mol = MultiMolecule.from_xyz(xyz_file)
        >>> rdf: pd.DataFrame = mol.init_rdf(atom_subset=atom_subset)
        >>> param: pd.DataFrame = estimate_lj(rdf)

        >>> print(param)
    
    The script provides the sigma values in Angstrom so we divided them by 10 to obtain the corresponding values in nm.
    In order to avoid atoms getting too close one from each other, we limited the sigma parameters with a miminal value (choosen to be 0.02nm lower than the initial value).



.. _charge: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/CHARGE.html#list_CHARGE
.. _epsilon: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_EPSILON
.. _sigma: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_SIGMA
