.. _simulation_box:

Preparing a Simulation Box
==========================

The goal of this tutorial is to prepare a 11.5 nm simulation box for Classical Molecular Dynamics (CMD) simulations. We have chosen to simulate the products of a synthesis procedure to obtain monodisperse CsPbBr\ :sub:`3`\ NCs, described by L. Protesescu *et al* in *Nano Letters*, **2015**, *15(6)*, 3692-3696.
The box specifically contains:
- One CsPbBr\ :sub:`3`\ core capped by 20% oleate (OA) and 20% oleylammonium (OLA) ligands;
- 2293 octadecene (ODA) molecules, used as solvent for the reaction;
- 287 oleylamine (OLAM) molecules, used as a reagent in the synthesis;
- 75 ionic oleate-oleylammonium couples (by-products obtained from the reaction).
    
Installation Requirements
-------------------------

This tutorial requires the download and use of the following programs:
- the **CAT**, **data-CAT** and **nano-CAT** packages for the construction of the NC model and of the organic molecules in the box. The relative `documentation <https://cat.readthedocs.io/en/latest/0_documentation.html#cat-documentation>`__ is hereby provided for more information on the keywords. 
- The **Packmol** package (see the following `link <http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml>`__ for insight on the installation and the keywords) for the construction of the coordinate file of the entire simulation box.
- The **Auto-FOX** `package <https://auto-fox.readthedocs.io/en/latest/includeme.html>`__, required for the construction of the Protein Structure File (.psf) containing the molecular-level information required to apply any force field to our simulation box;
- Ultimately, the **VMD** `software <https://www.ks.uiuc.edu/Research/vmd/>`__ has been employed as a molecular visualization program for the construction of the topology (.top) and the Gromacs (.gro) files to run the MD simulations.

Getting started - The molecules
-------------------------------

The first step towards creating the simulation box involves the creation of the .xyz files containing the geometries of the molecules. 
We'll achieve this by running a small .yaml script with **CAT**. First of all, we need to create an .xyz file according to our needs and to move the newly created .xyz file inside our working directory (see the `General Overview <https://cat.readthedocs.io/en/latest/1_get_started.html#default-settings>`_ for further information). We will then create a ``input_settings.yaml`` `input file <https://cat.readthedocs.io/en/latest/includeme.html#input-files>`_ in the working directory and customize it with the desired settings.
We hereby provide a .yaml input example for the construction of a cesium oleate molecule from scratch. Let's take a look:

.. code:: yaml

    path: null
    input_cores:
        - CsCl.xyz:
            guess_bonds: False
    input_ligands:
        - CCCCCCCC/C=C\CCCCCCCC(=O)[O-]
    optional:
    core:
        dirname: core
        anchor: Cl
        allignment: sphere
    ligand:
        dirname: ligand
        optimize: True
        split: True 
    qd:
        dirname: qd
        construct_qd: True
        optimize: False
            
The `path <https://cat.readthedocs.io/en/latest/2_path.html#path>`_, `input_cores & input_ligands <https://cat.readthedocs.io/en/latest/3_input_core_ligand.html#input-cores-input-ligands>`_ and  sections, together with the meaning of the `optional <https://cat.readthedocs.io/en/latest/4_optional.html#optional>`_ keywords and their relative `arguments <https://cat.readthedocs.io/en/latest/4_optional.html#arguments>`_, can be easily found inside the **CAT** `documentation <https://cat.readthedocs.io/en/latest/0_documentation.html#cat-documentation>`_. Let's take a look at them in detail:

1. ``path``: The `path <https://cat.readthedocs.io/en/latest/2_path.html#path>`_ section, as suggested, contains the path to the so-called working directory - i.e. where all the files are stored.
2. ``input_cores``: This section requires a little more insight: **CAT** is originally born as a program to build nanocrystals consisting of cores and ligands. When we build the .xyz files of our molecules, we thus treat them as if we're replacing the superficial ions of a "minimal" core made of two atoms.
This section contains the coordinates of the "minimal", biatomic core, specified by our .xyz file (``CsCl.xyz``). The ``guess_bonds: False`` keyword tells **CAT** that, since our core is ionic, it is not necessary the bonds and bond orders from the content of the .xyz file (i.e. it is not required to generate the internal coordinates of the system).
3. ``input_ligands``: This section contains information on both the structure and the chemistry of the ligand. This information is stored in its `SMILES <https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Description>`_ (Simplified molecular-input line-entry system) string, specifically ``CCCCCCCCC=CCCCCCCCC(=O)[O-]`` for oleate.
4. ``optional``: The `optional <https://cat.readthedocs.io/en/latest/4_optional.html#optional>`_ section contains three fairly similar subsections: ``core``, ``ligand``, ``qd``. The subsections contain keywords with several specifications, such as:

- the directories where inorganic cores/ligands/qd (or, in this case, the .xyz file containing our new molecule) will be stored (``optional.*.dirname``);
- whether or not their optimization is required (``optional.ligand.optimize`` and ``optional.*.optimize``);
- the dummy atom that needs to be replaced with the chosen ligand (``optional.*.anchor``);
- how the to-be attached ligands should be alligned with the core (``optional.*.allignment``). In our case, since we're building a molecule instead of a nanocrystal (NC), this key is mandatorily set to``optional.core.allignment: sphere``.
- whether or not to remove protons from the ligand (``optional.ligand.split``). Specifically, since the SMILES string we are using in the input (i.e. ``CCCCCCCCC=CCCCCCCCC(=O)[O-]``) refers to the anionic ligand, we will opt for ``optional.ligand.split: False``, so no protons have been removed from the ligand anchoring group. Conversely, if the SMILES is provided in the neutral form, then ``optional.ligand.split: True``, meaning that a proton is cleaved from the functional group (in this case carboxylate) to ensure that the ligand is still added in its anionic form. Note that the latter form is preferrable when the ligand present more than one functional group.  

In all cases, the ``*`` in the keywords accounts for the name of the subsection it refers to (i.e ``core``, ``ligand``, ``qd``).

We are finally ready to run CAT with the following command: ``init_cat input_settings.yaml``
After running **CAT** the .xyz file corresponding to our cesium oleate molecule can be found in the specified directory, 'qd'. Don't worry, the directory will be created from scratch if it does not yet exist!
After having renamed the .xyz file, we can just remove the counterion Cs from the molecule and we'll have the complete .xyz file for OA.

In a parallel fashion, the same script can be used to build the remaining .xyz files as follows:
- OLA molecules can be built by replacing the ligand SMILES string to that of oleylammonium (``CCCCCCCCC=CCCCCCCCC[NH3+]``) and the core anchor to ``Cl``, in order to obtain oleylammonium bromide. The Cl atom will then be replaced from the resulting .xyz file;
- OLAM molecules can be obtained in a similar fashion by replacing the ligand SMILES string to that of oleylammine (``CCCCCCCCC=CCCCCCCCCN``), by setting ``optional.ligand.split: False`` and the core anchor to ``Cl``, in order to obtain an .xyz containing oleylamine and bromide. The Cl atom will then be replaced from the resulting .xyz file to obtain OLAM.
- The QD can be built in a very similar fashion using this script. We invite you to refer to the `tutorial <https://nanotutorials.readthedocs.io/en/latest/1_build_qd.html>`__ for the step-by-step construction of the structure from scratch.

All of the remaining molecules (such as the CsCl.xyz and the .xyz file for ODA) can be built 
