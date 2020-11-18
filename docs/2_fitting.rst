.. _fitting:

Forcefield Optimization
=======================
The goal here is to obtain the classical forcefield parameters for a lead perovksite nanocrystal (NC) capped by carboxylate ligands. While these parameters are commonly available in literature for the ligands, in the QD field we need to construct our own parameters for a proper description of:

    * the ion-ion interactions inside the nanocrystal core;
    * the ligand anchoring group-core ions interactions at the nanocrystal surface.

To do that, we will start from an *ab-initio* Molecular Dynamics (MD) trajectory (NVT, 300K, 5ps) of a 2.5 nm sided cubic CsPbBr_3 NC capped by 50% of acetate molecules (see `tutorial <https://nanotutorials.readthedocs.io/en/latest/1_build_qd.html>`_ for the Quantum Dot construction).

Before starting the fitting we invite you to read the `documentation <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo.html>`_ relative to Adaptive Rate Monte Carlo (ARMC) in Auto-FOX.

First, let's have a look at the .yaml file containing our ARMC settings.

.. code:: yaml

    param:
        charge:
            param: charge
            Cs: 0.367090
            Pb: 0.823531
            Br: -0.383415
            C2O3: 0.247973
            O2D2: -0.284449
            constraints:
                - '0 < Cs < 1.5'
                - '0 < Pb < 2'
                - '-1.5 < Br < 0'

        lennard_jones:
            - param: epsilon
              unit: kjmol
              frozen:
                  guess: uff
            - param: sigma
              unit: nm
              Cs Cs:  0.400
              Cs Pb:  0.355
              Cs Br:  0.300
              Pb Pb:  0.520
              Pb Br:  0.250
              Br Br:  0.320
              C2O3 Cs: 0.295
              C2O3 Pb: 0.265
              C2O3 Br: 0.305
              O2D2 Cs: 0.250
              O2D2 Pb: 0.210
              O2D2 Br: 0.280
              constraints:
                  - 'Cs Cs   > 0.380'
                  - 'Cs Pb   > 0.335'
                  - 'Cs Br   > 0.280'
                  - 'Pb Pb   > 0.500'
                  - 'Pb Br   > 0.230'
                  - 'Br Br   > 0.300'
                  - 'C2O3 Cs > 0.275'
                  - 'C2O3 Pb > 0.245'
                  - 'C2O3 Br > 0.285'
                  - 'O2D2 Cs > 0.230'
                  - 'O2D2 Pb > 0.190'
                  - 'O2D2 Br > 0.260'
              frozen:
                  C331 Cs: 0.295
                  C331 Pb: 0.265
                  C331 Br: 0.305
                  HGA3 Cs: 0.255
                  HGA3 Pb: 0.270
                  HGA3 Br: 0.235

    psf:
        rtf_file: acetate.rtf
        ligand_atoms: [C, O, H]

    pes:
        rdf:
            func: FOX.MultiMolecule.init_rdf
            kwargs:
                atom_subset: [Cs, Pb, Br, O2D2]

    job:
        molecule: last5000.xyz

        geometry_opt:
            template: qmflows.templates.geometry.specific.cp2k_mm
            settings:
                prm: acetate.prm
        md:
            template: qmflows.templates.md.specific.cp2k_mm
            settings:
                prm: acetate.prm


Now, let's see in detail the contents of each section of our input file.

The param block
---------------
The ``"param"`` key contains all user-specified features concerning the to-be optimized parameters for the Coulomb potential (the charge_)
and the Lennard-Jones potential (epsilon_ & sigma_). Let's have a look at the relative sub-blocks:

1.  Coulomb potential

    .. code:: yaml

        charge:
            param: charge
            Cs: 0.367090
            Pb: 0.823531
            Br: -0.383415
            C2O3: 0.247973
            O2D2: -0.284449
            constraints:
                - '0 < Cs < 1.5'
                - '0 < Pb < 2'
                - '-1.5 < Br < 0'

    Here, the to-be optimized charges are those of the nanocrystal core ions (Cs, Pb, Br) and of the ligand anchoring group atoms (carboxylate group of the acetate, i.e. C2O3 and O2D2). Their initial values are obtained:

    * For the nanocrystal core ions, from a previous fitting procedure. You can simply use the most stable oxidation state of each ion if you don't have a better starting point.
    * For the anchoring group of the ligand, by adjusting the charges found in the .rtf file of the ligand to have an overall charge neutral system.
    In this case, the core ions charges are constrained to a certain range in order to keep the correct oxidation state (for example cations constrained to values higher than 0).

Let's move to the :code:`lennard_jones` block.

2.  Lennard-Jones potential

    This sub-block is divided itself in two components: epsilon_ and sigma_. Let's have a look at them:

    .. code:: yaml

            - param: epsilon
              unit: kjmol
              frozen:
                  guess: uff
    In our fitting the epsilon parameters treated as constants rather than to-be optimized variables (all frozen) and all the values are guessed using
    the `uff <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo.html#parameter-guessing>`_ procedure. Specifying the epsilon parameters (even without optimizing them) helps achieving a more accurate fitting.

    .. code:: yaml

            - param: sigma
              unit: nm
              Cs Cs:  0.400
              Cs Pb:  0.355
              Cs Br:  0.300
              Pb Pb:  0.520
              Pb Br:  0.250
              Br Br:  0.320
              C2O3 Cs: 0.295
              C2O3 Pb: 0.265
              C2O3 Br: 0.305
              O2D2 Cs: 0.250
              O2D2 Pb: 0.210
              O2D2 Br: 0.280
              constraints:
                  - 'Cs Cs   > 0.380'
                  - 'Cs Pb   > 0.335'
                  - 'Cs Br   > 0.280'
                  - 'Pb Pb   > 0.500'
                  - 'Pb Br   > 0.230'
                  - 'Br Br   > 0.300'
                  - 'C2O3 Cs > 0.275'
                  - 'C2O3 Pb > 0.245'
                  - 'C2O3 Br > 0.285'
                  - 'O2D2 Cs > 0.230'
                  - 'O2D2 Pb > 0.190'
                  - 'O2D2 Br > 0.260'
              frozen:
                  C331 Cs: 0.295
                  C331 Pb: 0.265
                  C331 Br: 0.305
                  HGA3 Cs: 0.255
                  HGA3 Pb: 0.270
                  HGA3 Br: 0.235

    Here we need to optimize the sigma parameters for the all pair interactions of interest (provided with the corresponding `atom pairs <https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_ATOMS>`_): 
    the ion-ion interactions inside the nanocrystal core (eg. Cs-Cs) and the acetate anchoring group-core ions interactions (eg. O2D2-Cs).
    The initial parameters for these pairs are obtained from the DFT trajectory by means of a small python script:

    .. code:: python

        >>> import pandas as pd
        >>> from FOX import MultiMolecule, example_xyz, estimate_lj

        >>> xyz_file: str = 'last5000.xyz' # path of DFT trajectory
        >>> atom_subset = ['Cs', 'Pb', 'Br', 'C', 'O', 'H'] # core ions and acetate atoms

        >>> mol = MultiMolecule.from_xyz(xyz_file)
        >>> rdf: pd.DataFrame = mol.init_rdf(atom_subset=atom_subset)
        >>> param: pd.DataFrame = estimate_lj(rdf)

        >>> print(param)


The script provides the sigma values in Angstrom so we divided them by 10 to obtain the corresponding values in nm.
In order to avoid atoms getting too close one from each other, we constrained the sigma parameters to be higher than a miminal value (choosen to be exactly 0.02 nm lower than the initial value).
Finally, in the ``"frozen"`` subsection,  we specified the sigma values for the acetate methyl group - core ions interactions (eg. C331 Cs) as frozen components 
(so without optimizing them). Similarly to the to-be optimized sigmas, the corresponding frozen values are taken from the output of the python script shown above. 
Once again, this specification results in a smoother fitting procedure.
    
The psf block
-------------
    .. code:: yaml

    psf:
           rtf_file: acetate.rtf
           ligand_atoms: [C, O, H]

This 

.. _charge: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/CHARGE.html#list_CHARGE
.. _epsilon: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_EPSILON
.. _sigma: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_SIGMA
