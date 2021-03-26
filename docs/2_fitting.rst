.. _fitting:

Forcefield Optimization
=======================
The goal of this tutorial is to obtain the classical forcefield parameters for a lead perovksite core and for the nanocrystal (NC) obtained by capping the fitted core with carboxylate ligands. In fact, the construction of a forcefield for a Quantum Dot (QD) requires parameters for a proper description of:

    * the ion-ion interactions inside the nanocrystal core;
    * the ligand anchoring group-core ions interactions at the nanocrystal surface.
    
The third "category" of parameters, accounting for the organic ligands, are commonly available in literature and we thus won't need to fit them.



To fit the parameters, we will start from an *ab-initio* Molecular Dynamics (MD) trajectory (NVT, 300K, 5ps) of a 2.3 nm sided cubic CsPbBr_3 core (see `tutorial <https://nanotutorials.readthedocs.io/en/latest/1_build_qd.html>`_ for the Quantum Dot construction using **CAT**).

Before starting the fitting we invite you to read the `documentation <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo.html>`_ relative to Adaptive Rate Monte Carlo (ARMC) in Auto-FOX.

The inorganic core
---------------

First, let's have a look at the .yaml file containing our ARMC settings for the fitting of the forcefield parameters of our CsPbBr_3 core from its previously calculated quantum mechanic, Density Functional Theory (DFT) trajectory.

    .. code:: yaml
    
        param:
        charge:
        param: charge
        Cs: 0.4174
        Pb: 0.8348
        Br: -0.4174
        constraints:
            - '0 < Cs < 1.5'
            - '0 < Pb < 2'
            - '-1.5 < Br < 0'
            - 'Cs == -1 * Br'
            - 'Pb == -2 * Br'

        lennard_jones:
            - param: epsilon
              unit: kjmol
              frozen: 
                  guess: uff
            - param: sigma
              unit: nm
              Cs Cs:  0.453
              Cs Pb:  0.367
              Br Cs:  0.363
              Pb Pb:  0.610
              Br Pb:  0.298
              Br Br:  0.369
              constraints:
                  - 'Cs Cs   > 0.433'
                  - 'Cs Pb   > 0.347'
                  - 'Br Cs   > 0.343'
                  - 'Pb Pb   > 0.590'
                  - 'Br Pb   > 0.278'
                  - 'Br Br   > 0.349'

        pes:
            rdf:
                func: FOX.MultiMolecule.init_rdf
                kwargs:
                    atom_subset: [Cs, Pb, Br]
    
        job:
            molecule: 2.3nm_cspbbr3_NVT_300K-pos-1.xyz
        
            md_settings:
                template: qmflows.templates.md.specific.cp2k_mm
                settings:
                    input:
                        global:
                            print_level: LOW
                        force_eval:
                            mm:
                              poisson:
                                 periodic: xyz
                                 ewald:
                                   ewald_type: spme
                                   gmax: '62 62 62'
                                   o_spline: 4
                            subsys:
                                cell:
                                    abc: '[angstrom] 100.0 100.0 100.0'
                                    periodic: xyz
        
                        motion:
                            print:
                                restart:
                                   each:
                                      md: 10
                                trajectory:
                                   each:
                                      md: 10
                                velocities:
                                   each:
                                      md: 10
                                forces:
                                   each:
                                      md: 10
                            md:
                                ensemble: NVT
                                temperature: 300.0
                                timestep: 2.5
                                steps: 10000
                                thermostat:
                                    type: csvr
                                    csvr:
                                        timecon: 10000
        
        monte_carlo:
            type: FOX.armc.ARMC
            iter_len: 50000
            sub_iter_len: 10
            logfile: armc.log
            hdf5_file: armc.hdf5
            path: ./
            folder: MM_MD_workdir
            keep_files: True

Now, let's see in detail the contents of each section of our input file.

The param block
---------------
The ``"param"`` key contains all user-specified features concerning the to-be optimized parameters for the Coulomb potential (the charge_)
and the Lennard-Jones potential (epsilon_ & sigma_). Let's have a look at the relative sub-blocks:

1.  **Coulomb potential**

    .. code:: yaml
    
        param:
        charge:
        param: charge
        Cs: 0.4174
        Pb: 0.8348
        Br: -0.4174
        constraints:
            - '0 < Cs < 1.5'
            - '0 < Pb < 2'
            - '-1.5 < Br < 0'
            - 'Cs == -1 * Br'
            - 'Pb == -2 * Br'

    Here, the to-be optimized charges are those of the nanocrystal core ions (Cs, Pb, Br). Their initial values are usually obtained from their DFT trajectory. You can simply use the most stable oxidation state of each ion if you don't have a better starting point.
    In this case, the core ions charges are constrained to a certain range in order to keep the correct oxidation state (for example cations constrained to values higher than 0), as well as the prerequisite of the overall neutrality of the system. Additional constraints are added to ensure that the ions correctly balance each other in case of the detachment of a neutral species, i.e. CsBr and PbBr_2, from the surface of the core.

Let's move to the :code:`lennard_jones` block.

2.  **Lennard-Jones potential**

    This sub-block is divided in two further components: epsilon_ and sigma_. Let's have a look at them:

    .. code:: yaml

            - param: epsilon
              unit: kjmol
              frozen:
                  guess: uff
                  
    In our fitting the epsilon parameters treated as constants rather than to-be optimized variables (all frozen) and all the values are guessed using
    the `uff <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo.html#parameter-guessing>`_ procedure, as specified by their so-called ``"frozen"`` subsection. Specifying the epsilon parameters (even without optimizing them) helps achieving a more accurate fitting.

    .. code:: yaml

            - param: sigma
              unit: nm
              Cs Cs:  0.453
              Cs Pb:  0.367
              Br Cs:  0.363
              Pb Pb:  0.610
              Br Pb:  0.298
              Br Br:  0.369
              constraints:
                  - 'Cs Cs   > 0.433'
                  - 'Cs Pb   > 0.347'
                  - 'Br Cs   > 0.343'
                  - 'Pb Pb   > 0.590'
                  - 'Br Pb   > 0.278'
                  - 'Br Br   > 0.349'
                  
    Here we need to optimize the sigma parameters for the pair interactions of interest (provided with the corresponding `atom pairs <https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_ATOMS>`_), i.e.
    the ion-ion interactions inside the nanocrystal core (eg. Cs-Cs).
    The initial parameters for these pairs are obtained from the DFT trajectory by means of a small python script:

    .. code:: python

        >>> import pandas as pd
        >>> from FOX import MultiMolecule, example_xyz, estimate_lj

        >>> xyz_file: str = '2.3nm_cspbbr3_NVT_300K-pos-1.xyz' # path of DFT trajectory
        >>> atom_subset = ['Cs', 'Pb', 'Br'] # core ions

        >>> mol = MultiMolecule.from_xyz(xyz_file)
        >>> rdf: pd.DataFrame = mol.init_rdf(atom_subset=atom_subset)
        >>> param: pd.DataFrame = estimate_lj(rdf)

        >>> print(param)


The script provides the sigma values in Angstrom so we divided them by 10 to obtain the corresponding values in nm.
In order to avoid atoms getting too close one from each other, we constrained the sigma parameters to be higher than a minimal value (choosen to be exactly 0.02 nm lower than the initial value).

The pes block
-------------
The `pes <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo_args.html?highlight=rtf#pes>`_ block contains the setting and descriptors aimed at the construction of the Potential Energy Surface (PES) of the atoms we aim to fit, specified in the kwargs_ subsection. We chose to calculate their radial distribution function (rdf_).

    .. code:: yaml
    
        pes:
            rdf:
                func: FOX.MultiMolecule.init_rdf
                kwargs:
                    atom_subset: [Cs, Pb, Br]
                    

The job block
-------------
The `job <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo_args.html?highlight=job#job>`_ section is divided into two subsections:

    * ``molecule``, containing the reference .xyz file with the reference QM rdf;
    * ``md_settings``, specifying the the settings of the calculation we want to perform (in our case the MD simulations).     

    .. code:: yaml
    
        job:
            molecule: 2.3nm_cspbbr3_NVT_300K-pos-1.xyz
        
            md_settings:
                template: qmflows.templates.md.specific.cp2k_mm
                settings:
                    input:
                        global:
                            print_level: LOW
                        force_eval:
                            mm:
                              poisson:
                                 periodic: xyz
                                 ewald:
                                   ewald_type: spme
                                   gmax: '62 62 62'
                                   o_spline: 4
                            subsys:
                                cell:
                                    abc: '[angstrom] 100.0 100.0 100.0'
                                    periodic: xyz
        
                        motion:
                            print:
                                restart:
                                   each:
                                      md: 10
                                trajectory:
                                   each:
                                      md: 10
                                velocities:
                                   each:
                                      md: 10
                                forces:
                                   each:
                                      md: 10
                            md:
                                ensemble: NVT
                                temperature: 300.0
                                timestep: 2.5
                                steps: 10000
                                thermostat:
                                    type: csvr
                                    csvr:
                                        timecon: 10000


This section containts the actual parameters that will figure in the CP2K input file: for further inquiries on the keywords, we invite you to refer to the relative `documentation <https://manual.cp2k.org/cp2k-7_1-branch/CP2K_INPUT.html>`_. These parameters can be tailored according to need: for example, in our case, we tailored the value of ``gmax`` on the dimension of our cubic cell (whose periodic parameters are thus provided as ``abc``) and we chose which properties - the trajectory, velocities and forces - to print over each MD run depending on the future calculations we aimed to perform. Moreover, we performed NVT MD simulations on systems at room temperature and, in the absence of organic molecules, we opted for 2.5 fs integration timesteps. 

The monte_carlo block
-------------
The `monte_carlo <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo_args.html?highlight=md_settings#monte-carlo>`_ block contains all the settings required to operate the Monte Carlo procedure (in our case, we are making use of the `Adaptive Rate Monte Carlo <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo.html#addaptive-rate-monte-carlo>`_ algorithm), including the total number of iterations and sub_iterations in the procedure, the name and path of the logfile containing the summary of the performed jobs and their respective errors calculated through a comparison with our chosen PES descriptor (rdf), the paths of the working directory and whether or not the directories containing the single MD jobs are being kept in the main working directory (``keep_files: True`` or ``False``).

    .. code:: yaml
    
        monte_carlo:
            type: FOX.armc.ARMC
            iter_len: 50000
            sub_iter_len: 10
            logfile: armc.log
            hdf5_file: armc.hdf5
            path: ./
            folder: MM_MD_workdir
            keep_files: True

We will thus perform the fitting procedure by opening our conda environment containing **Auto-FOX** and computing the command prompt ``init_armc settings.yaml``.

The nanocrystal
---------------
Once we obtain reliable parameters (i.e. when the comparison between our reference function, the MM radial distribution function calculated with the fitted parameters, and the QM-computed radial distribution function displays a very low error), we can use these parameters as a starting point to build a new .yaml input for the fitting of the forcefield parameters of the NC obtained by capping the fitted CsPbBr_3 core with acetate ligands. Let's have a brief look at the new input file.

    .. code:: yaml
    
        param:
            charge:
                param: charge
                Cs: 0.4
                Pb: 0.8
                Br: -0.4
                C2O3: 0.25
                O2D2: -0.275
                constraints:
                    - '0 < Cs < 1.5'
                    - '0 < Pb < 2'
                    - '-1.5 < Br < 0'
                    - 'Cs == -1 * $LIGAND'
                    - 'Pb == -2 * $LIGAND'
                    - 'Cs == -1 * Br'
                    - 'Pb == -2 * Br'
        
            lennard_jones:
                - param: epsilon
                  unit: kjmol
                  frozen:
                      guess: uff
                - param: sigma
                  unit: nm
                  Cs Cs:  0.433
                  Cs Pb:  0.362
                  Br Cs:  0.389
                  Pb Pb:  0.636
                  Br Pb:  0.316
                  Br Br:  0.369
                  C2O3 Cs: 0.437
                  C2O3 Pb: 0.348
                  Br C2O3: 0.383
                  Cs O2D2: 0.331
                  O2D2 Pb: 0.264
                  Br O2D2: 0.369
                  constraints:
                      - 'Cs Cs   > 0.523'
                      - 'Cs Pb   > 0.342'
                      - 'Br Cs   > 0.369'
                      - 'Pb Pb   > 0.616'
                      - 'Br Pb   > 0.296'
                      - 'Br Br   > 0.349'
                      - 'C2O3 Cs > 0.417'
                      - 'C2O3 Pb > 0.328'
                      - 'Br C2O3 > 0.363'
                      - 'Cs O2D2 > 0.311'
                      - 'O2D2 Pb > 0.244'
                      - 'Br O2D2 > 0.349'
                  frozen:
                      C331 Cs: 0.295
                      C331 Pb: 0.265
                      Br C331: 0.305
                      Cs HGA3: 0.255
                      HGA3 Pb: 0.270
                      Br HGA3: 0.235
        
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
        
            md_settings:
                template: qmflows.templates.md.specific.cp2k_mm
                settings:
                    prm: acetate.prm
                    input:
                        global:
                            print_level: LOW
                        force_eval:
                            mm:
                              poisson:
                                 periodic: xyz
                                 ewald:
                                   ewald_type: spme
                                   gmax: '62 62 62'
                                   o_spline: 4
                            subsys:
                                cell:
                                    abc: '[angstrom] 100.0 100.0 100.0'
                                    periodic: xyz
                        motion:
                            print:
                                cell:
                                   each:
                                      md: 10
                                restart:
                                   each:
                                      md: 10
                                trajectory:
                                   each:
                                      md: 10
                                velocities:
                                   each:
                                      md: 10
                                forces:
                                   each:
                                      md: 10
                            md:
                                ensemble: NVT
                                temperature: 300.0
                                timestep: 1
                                steps: 10000
                                thermostat:
                                    type: csvr
                                    csvr:
                                        timecon: 10000
                                print:
                                    energy:
                                        each:
                                           md: 10
        
        monte_carlo:
            type: FOX.armc.ARMC
            iter_len: 50000
            sub_iter_len: 10
            logfile: armc.log
            hdf5_file: armc.hdf5
            path: ./
            folder: MM_MD_workdir
            keep_files: True

The yaml code above shows a clear resemblance to the one used for the core, except for a few key differences. We hereby provide a brief comparison of their features.

The param block
---------------

    .. code:: yaml
    
        param:
            charge:
                param: charge
                Cs: 0.4
                Pb: 0.8
                Br: -0.4
                C2O3: 0.25
                O2D2: -0.275
                constraints:
                    - '0 < Cs < 1.5'
                    - '0 < Pb < 2'
                    - '-1.5 < Br < 0'
                    - 'Cs == -1 * $LIGAND'
                    - 'Pb == -2 * $LIGAND'
                    - 'Cs == -1 * Br'
                    - 'Pb == -2 * Br'
                    
          
Here, the Coulomb potential sub-block shows both the charges of the nanocrystal core ions (Cs, Pb, Br) and those of the ligand anchoring group atoms (in this specific case, the carboxylate group of the acetate, i.e. C2O3 and O2D2). Their initial values are usually obtained: 
    
    * For the nanocrystal core ions, from the approximated results of the previous fitting procedure used for the inorganic core or by their most stable oxidation state, in absence of more accurate parameters.
    * For the anchoring group of the ligand, by adjusting the charges (found both in the .yaml input and in the .rtf file of the ligand) to achieve the overall charge neutrality of the system. More specifically, the total charge of the ligand needs to equal the charge of the atom it replaces: in this specific case, our ligand is an acetate group, and it thus needs to balance the charge of the Br atom (-0.4). We will provide an example of this procedure in the following section.
    
    .. code:: yaml    
    
            lennard_jones:
                - param: epsilon
                  unit: kjmol
                  frozen:
                      guess: uff
                - param: sigma
                  unit: nm
                  Cs Cs:  0.553
                  Cs Pb:  0.367
                  Br Cs:  0.363
                  Pb Pb:  0.610
                  Br Pb:  0.298
                  Br Br:  0.379
                  C2O3 Cs: 0.437
                  C2O3 Pb: 0.348
                  Br C2O3: 0.383
                  Cs O2D2: 0.331
                  O2D2 Pb: 0.264
                  Br O2D2: 0.369
                  constraints:
                      - 'Cs Cs   > 0.523'
                      - 'Cs Pb   > 0.337'
                      - 'Br Cs   > 0.333'
                      - 'Pb Pb   > 0.580'
                      - 'Br Pb   > 0.268'
                      - 'Br Br   > 0.349'
                      - 'C2O3 Cs > 0.407'
                      - 'C2O3 Pb > 0.318'
                      - 'Br C2O3 > 0.353'
                      - 'Cs O2D2 > 0.301'
                      - 'O2D2 Pb > 0.234'
                      - 'Br O2D2 > 0.339'
                  frozen:
                      C331 Cs: 0.295
                      C331 Pb: 0.265
                      Br C331: 0.305
                      Cs HGA3: 0.255
                      HGA3 Pb: 0.270
                      Br HGA3: 0.235

    In the :code:`lennard_jones` block we will need to optimize the sigma parameters for all the `atom pair <https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_ATOMS>`_ interactions of interest, including both the ion-ion interactions inside the nanocrystal core (eg. Cs-Cs) and the acetate anchoring group-core ions interactions (eg. O2D2-Cs). In addition, the sigmas between the ions in the inorganic core and the ligand atoms which are not in the anchoring group are treated as frozen (non-optimized, constant parameters): their values are thus inserted in the ``"frozen"`` subsection. The initial parameters for these pairs are obtained from the DFT trajectory by means of a small python script:

    .. code:: python

        >>> import pandas as pd
        >>> from FOX import MultiMolecule, example_xyz, estimate_lj

        >>> xyz_file: str = 'last5000.xyz' # path of DFT trajectory
        >>> atom_subset = ['Cs', 'Pb', 'Br', 'C', 'O', 'H'] # core ions and acetate atoms

        >>> mol = MultiMolecule.from_xyz(xyz_file)
        >>> rdf: pd.DataFrame = mol.init_rdf(atom_subset=atom_subset)
        >>> param: pd.DataFrame = estimate_lj(rdf)

        >>> print(param)


In this case, the output of this python script provides both the sigma values for both to the to-be optimized sigmas and the frozen components. Once again, in order to avoid atoms getting too close one from each other, we constrained the sigma parameters to be 0.02 nm lower than their estimated value: resulting in a smoother fitting procedure.

The psf block
-------------

The `psf <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo_args.html?highlight=psf#psf>`_ section contains the settings required for the construction of the protein structure files. In our case the required data is the name of the .rtf file and a list identifying the atoms of the ligands.

    .. code:: yaml

        psf:
               rtf_file: acetate.rtf
               ligand_atoms: [C, O, H]
           
The .rtf file is used for assigning atom types and charges to ligands. In fact, any information on the ligand which isn't contained in the .yaml input is read from its .rtf file. Let's see an example of its structure in detail for our acetate ligands:

::

  harmm RTF built by MATCH
  *
    22     0
  MASS   122 C2O3  12.01100 C
  MASS   123 C331  12.01100 C
  MASS   124 HGA3  1.008000 H
  MASS   125 O2D2  15.99900 O
  
  AUTO ANGLES DIHE
  
  RESI  LIG   -1.000000
  GROUP
  ATOM C    C331  -0.370000
  ATOM C2   C2O3   0.288746
  ATOM O    O2D2  -0.328684
  ATOM O5   O2D2   0.288746
  ATOM H6   HGA3   0.090000
  ATOM H7   HGA3   0.090000
  ATOM H    HGA3   0.090000
  BOND C2   C
  BOND C    H
  BOND C    H6
  BOND C    H7
  BOND C2   O
  BOND C2   O5
  IMPR C2   C    O    O5
  PATCH FIRST NONE LAST NONE
  
  END

As we can see, this file contains a block indicating the masses of the ligand atoms and one containing their charges. The line ``RESI LIG -1.000000`` highlights the total charge on each ligand, which is the sum of the charges of its constituent atoms (i.e. -0.37 + 0.288746 + (-0.328684) + 0.288746 + 3*0.09 = -1).
Since any information on the ligand which isn't contained in the .yaml input is read from its .rtf file, we can modulate the charge for our anchoring group (``C2O3`` and ``O2D2``) in our yaml input, and they will be overwritten. More specifically, the total charge on each acetate molecule needs to balance the charge we indicated for Br atoms (i.e. ``Br  -0.4``), so that the charge of the system is kept neutral during the replacement. This means that the sum of the charges needs to be adjusted to satisfy the relationship: -0.37 + C2O3 + 2O2D2 + 3*0.09 = -0.4. We have thus chosen the values ``C2O3  0.25`` and ``O2D2  -0.275`` in the .yaml input because they satisfied these requirements mantaining the correct proportions between the atoms in the anchoring group.

The job block
-------------

        job:
            molecule: last5000.xyz
        
            md_settings:
                template: qmflows.templates.md.specific.cp2k_mm
                settings:
                    prm: acetate.prm
                    input:
                        global:
                            print_level: LOW
                        force_eval:
                            mm:
                              poisson:
                                 periodic: xyz
                                 ewald:
                                   ewald_type: spme
                                   gmax: '62 62 62'
                                   o_spline: 4
                            subsys:
                                cell:
                                    abc: '[angstrom] 100.0 100.0 100.0'
                                    periodic: xyz
                        motion:
                            print:
                                cell:
                                   each:
                                      md: 10
                                restart:
                                   each:
                                      md: 10
                                trajectory:
                                   each:
                                      md: 10
                                velocities:
                                   each:
                                      md: 10
                                forces:
                                   each:
                                      md: 10
                            md:
                                ensemble: NVT
                                temperature: 300.0
                                timestep: 1
                                steps: 10000
                                thermostat:
                                    type: csvr
                                    csvr:
                                        timecon: 10000
                                print:
                                    energy:
                                        each:
                                           md: 10

The main differences with the previous `job <https://auto-fox.readthedocs.io/en/latest/4_monte_carlo_args.html?highlight=job#job>`_ section are:

1. The presence of the ``settings.prm`` subsection, containing the homonymous file for the ligand;
2. The choice of a 1 fs timestep in the MDs, which is motivated by the need of an appropriate description of the vibration of the organic bonds in the ligands.

The remainder of the sections are structured in a parallel fashion to the previous input. We will once again perform the fitting procedure by opening our conda environment containing **Auto-FOX**, FOX, and computing the command prompt ``init_armc settings.yaml``.

.. _charge: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/CHARGE.html#list_CHARGE
.. _epsilon: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_EPSILON
.. _sigma: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES.html#list_SIGMA
.. _kwargs: https://auto-fox.readthedocs.io/en/latest/4_monte_carlo_args.html?highlight=rtf#pes.block.kwargs
.. _rdf: https://auto-fox.readthedocs.io/en/latest/1_rdf.html?highlight=init_rdf#radial-angular-distribution-function
