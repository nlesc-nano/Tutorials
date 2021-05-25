.. _simulation_box:

Preparing a Simulation Box
==========================

The goal of this tutorial is to prepare a 11.5 nm simulation box for Classical Molecular Dynamics (CMD) simulations. We have chosen to simulate the products of a synthesis procedure to obtain monodisperse CsPbBr\ :sub:`3`\  NCs, described by L. Protesescu *et al* in *Nano Letters*, **2015**, *15(6)*, 3692-3696.
The box specifically contains:

- One CsPbBr\ :sub:`3`\ core capped by 20% oleate (OA) and 20% oleylammonium (OLA) ligands;
- 75 ionic oleate-oleylammonium couples (by-products obtained from the reaction);
- 287 oleylamine (OLAM) molecules, used as a reagent in the synthesis;
- 2293 octadecene (ODA) molecules, used as solvent for the reaction;

    
Installation Requirements
-------------------------

This tutorial requires the download and use of the following programs:

- The **CAT**, **data-CAT** and **nano-CAT** packages for the construction of the NC model and of the organic molecules in the box. The relative `documentation <https://cat.readthedocs.io/en/latest/0_documentation.html#cat-documentation>`__ is hereby provided for more information on the keywords. 
- The **Packmol** package (see the following `link <http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml>`__ for insight on the installation and the keywords) for the construction of the coordinate file of the entire simulation box.
- The **Auto-FOX** `package <https://auto-fox.readthedocs.io/en/latest/includeme.html>`__, required for the construction of the Protein Structure File (.psf);
- Ultimately, the **VMD** `software <https://www.ks.uiuc.edu/Research/vmd/>`__ has been employed as a molecular visualization program for the construction of the topology (.top) and the Gromacs (.gro) files to run the MD simulations.

Other programs and/or packages will be mentioned over the course of this tutorial, but are not mentioned here as their installation is not required. In those cases, you can either refer to online pages performing the same function (e.g. online converters) or use the package you're most comfortable with (e.g. Molden, Gaussian etc. to obtain .xyz files).

Getting started - The molecules
-------------------------------
The first step towards creating the simulation box involves the creation of the files containing the geometries and the structural information regarding the molecules in our system. Let's see how to obtain the different files.

The .xyz files
^^^^^^^^^^^^^^
First of all, let's create the .xyz files containing the geometries of our molecules. 
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
            
The `path <https://cat.readthedocs.io/en/latest/2_path.html#path>`_, `input_cores & input_ligands <https://cat.readthedocs.io/en/latest/3_input_core_ligand.html#input-cores-input-ligands>`_ and  sections, together with the meaning of the `optional <https://cat.readthedocs.io/en/latest/4_optional.html#optional>`_ keywords and their relative `arguments <https://cat.readthedocs.io/en/latest/4_optional.html#arguments>`_, can be easily found inside the **CAT** `documentation <https://cat.readthedocs.io/en/latest/0_documentation.html#cat-documentation>`_. Let's take a look at some of them in detail:

1. ``path``: The `path <https://cat.readthedocs.io/en/latest/2_path.html#path>`_ section, as suggested, contains the path to the so-called working directory - i.e. where all the files are stored.
2. ``input_cores``: This section requires a little more insight: **CAT** is originally born as a program to build nanocrystals consisting of cores and ligands. When we build the .xyz files of our molecules, we thus treat them as if we're replacing the superficial ions of a "minimal" core made of two atoms. This section contains the coordinates of the "minimal", biatomic core, specified by our .xyz file (``CsCl.xyz``). The ``guess_bonds: False`` keyword tells **CAT** that, since our core is ionic, it is not necessary the bonds and bond orders from the content of the .xyz file (i.e. it is not required to generate the internal coordinates of the system).
3. ``input_ligands``: This section contains information on both the structure and the chemistry of the ligand. This information is stored in its canonical `SMILES <https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Description>`_ (Simplified molecular-input line-entry system) string, specifically ``CCCCCCCC/C=C\CCCCCCCC(=O)[O-]`` for oleate.
4. ``optional``: The `optional <https://cat.readthedocs.io/en/latest/4_optional.html#optional>`_ section contains three fairly similar subsections: ``core``, ``ligand``, ``qd``. The subsections contain keywords with several specifications, such as:

- the directories where inorganic cores/ligands/qd (or, in this case, the .xyz file containing our new molecule) will be stored (``optional.*.dirname``);
- the dummy atom that needs to be replaced with the chosen ligand (``optional.*.anchor``);
- how the to-be attached ligands should be alligned with the core (``optional.*.allignment``). In our case, since we're building a molecule instead of a nanocrystal (NC), this key is mandatorily set to``optional.core.allignment: sphere``.
- whether or not to remove protons from the ligand (``optional.ligand.split``). Specifically, since the SMILES string we are using in the input (i.e. ``CCCCCCCC/C=C\CCCCCCCC(=O)[O-]``) refers to the anionic ligand, we will opt for ``optional.ligand.split: False``, so no protons have been removed from the ligand anchoring group. Conversely, if the SMILES is provided in the neutral form, then ``optional.ligand.split: True``, meaning that a proton is cleaved from the functional group (in this case carboxylate) to ensure that the ligand is still added in its anionic form. Note that the latter form is preferrable when the ligand present more than one functional group.  

In all cases, the ``*`` in the keywords accounts for the name of the subsection it refers to (i.e ``core``, ``ligand``, ``qd``).

We are finally ready to run CAT with the following command: ``init_cat input_settings.yaml``
After running **CAT** the .xyz file corresponding to our cesium oleate molecule can be found in the specified directory, 'qd'. Don't worry, the directory will be created from scratch if it does not yet exist!
After having renamed the .xyz file, we can just remove the counterion Cs from the molecule and we'll have the complete .xyz file for OA.

In a parallel fashion, the same script can be used to build the remaining .xyz files as follows:
- OA+OLA molecules (i.e. our ionic oleate-oleylammonium couples) can be obtained by means of a similar script, which we hereby report and comment briefly:

.. code:: yaml

    path: null
    input_cores:
        - rbcl.xyz:
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
            split: False
    
        qd:
            dirname: qd
            construct_qd: True
            multi_ligand:
               ligands:
                 - CCCCCCCC/C=C\CCCCCCCC[NH3+]
               anchor:
                 - Rb
            optimize: False

the only difference from the previous script is the presence of the ``optional.qd.multiligand`` key and of its relative specifications. All the keys under this section are completely parallel to the ``optional.ligand`` key block: Rb atoms are being replaced by oleylammonium molecules.
**Please note** that, in order to work effectively, this block acccepts SMILES strings by assuming a ``split: True`` specification.
- OLAM can be obtained by replacing the ligand SMILES string in the first script to that of oleylammine (``CCCCCCCC/C=C\CCCCCCCCN``), by setting ``optional.ligand.split: False`` and the core anchor to ``Cl``, in order to obtain an .xyz containing oleylamine and bromide. The Cl atom will then be replaced from the resulting .xyz file to obtain OLAM.
- The QD can be built in a very similar fashion using this script. We invite you to refer to the `tutorial <https://nanotutorials.readthedocs.io/en/latest/1_build_qd.html>`__ for the step-by-step construction of the structure from scratch.

All of the remaining molecules (such as the CsCl.xyz and the .xyz file for ODA) can be built using any (commonly available) molecular structure processing program, such as `Molden <https://www3.cmbi.umcn.nl/molden/>`__.
To sum up, in our case, we have now successfully built **these** files (the names have been chosen to represent their chemical formula for simplicity):

- qd.xyz, containing our ligand-capped NC;
- oaola.xyz;
- olam.xyz;
- octadecene.xyz.

Other file extensions
^^^^^^^^^^^^^^^^^^^^^

Now that we've obtained our .xyz files, we need to convert it to other extensions in order to complete our 3D structure with additional, "missing" information. Let's see how to obtain the other files:

1. *.pdb file*: The Protein Data Bank (.pdb) extension provides a description of the atomic coordinates, secondary structure assignments and atomic connectivity of a molecule. An .xyz file can be easily converted to this format by means of `Open Babel <https://openbabel.org/docs/dev/Installation/install.html>`__, a commonly employed chemical format converter. You can follow this link and follow the instructions for the installation (or just look for any Open Babel-based format converters available online). Once the program is correctly installed, the .xyz files can be converted to the .pdb format by running this simple command for each _organic_ molecule (*NOTE that this step does **NOT** apply to our qd.xyz file): ``obabel -ixyz file.xyz -opdb file.pdb``.
To sum up, we will now have the following .pdb files:

    - oaola.pdb;
    - olam.pdb;
    - octadecene.pdb.
    
2. *.prm and .rtf files*: Each .pdb file we created now needs to be converted to the following formats:

    - The CHARMM forcefield Parameter (.prm) file, including all of the numerical constants needed to evaluate forces and energies;
    - The Residue Topology File (.rtf) This file defines the main groups (atoms, properties, bond and charge information) for our molecular structures.
    
These formats can be easily obtained from our .pdb files by inserting our .pdb files in `MATCH <https://openbabel.org/docs/dev/Installation/install.html>`__. This online server will convert our files into the three required formats, which we will download as a zipped directory. We will first of all need to rename the new files to match their molecular formulas (2 for each .pdb file, for a total of 6 new files).
The .rtf files are ready for our next step, so we can put them aside for the present moment. Let's instead focus on the .prm files. An example of a .prm file (here we chose OLAM) looks like this:

.. code:: yaml

    * prm file built by MATCH
    *
    
    BONDS
    C324   N3P3   200.00     1.4800
    HGP2   N3P3   403.00     1.0400
    C321   C324   222.50     1.5300
    C324   HGA2   284.50     1.1000
    C321   C321   222.50     1.5300
    C321   HGA2   309.00     1.1110
    C321   C331   222.50     1.5280
    C2D1   C321   365.00     1.5020
    C331   HGA3   322.00     1.1110
    C2D1   C2D1   440.00     1.3400
    C2D1   HGA4   360.50     1.1000
    
    ANGLES
    C321   C324   N3P3   67.70      110.00
    HGA2   C324   N3P3   45.00      107.50
    C324   N3P3   HGP2   30.00      109.50
    HGP2   N3P3   HGP2   44.00      109.50
    C321   C321   C324   58.35      110.50
    HGA2   C321   C324   26.50      110.10
    C321   C324   HGA2   26.50      111.80
    HGA2   C324   HGA2   35.50      109.00
    HGA2   C321   C321   26.50      110.10
    C321   C321   C321   58.35      113.60
    HGA2   C321   HGA2   35.50      109.00
    C321   C321   C331   58.00      115.00
    C321   C321   C2D1   32.00      112.20
    HGA2   C321   C331   34.60      110.10
    C321   C331   HGA3   34.60      110.10
    C2D1   C2D1   C321   48.00      123.50
    HGA4   C2D1   C321   40.00      116.00
    C2D1   C321   HGA2   45.00      111.50
    HGA3   C331   HGA3   35.50      108.40
    HGA4   C2D1   C2D1   52.00      119.50
    
    DIHEDRALS
    C321   C321   C324   N3P3   0.1950     3      0.00
    HGA2   C321   C324   N3P3   0.1950     3      0.00
    C321   C324   N3P3   HGP2   0.1000     3      0.00
    HGA2   C324   N3P3   HGP2   0.1000     3      0.00
    C321   C321   C321   C324   0.1950     3      0.00
    HGA2   C321   C321   C324   0.1950     3      0.00
    C321   C321   C324   HGA2   0.1950     3      0.00
    HGA2   C321   C324   HGA2   0.1950     3      0.00
    HGA2   C321   C321   C321   0.1950     3      0.00
    HGA2   C321   C321   HGA2   0.2200     3      0.00
    C321   C321   C321   C321   0.14975    3      180.00
    C321   C321   C321   C321   0.09458    4      0.00
    C321   C321   C321   C321   0.11251    5      0.00
    C321   C321   C321   C321   0.06450    2      0.00
    C321   C321   C321   C331   0.08133    3      180.00
    C321   C321   C321   C331   0.10824    4      0.00
    C321   C321   C321   C331   0.20391    5      0.00
    C321   C321   C321   C331   0.15051    2      0.00
    C321   C321   C321   C2D1   0.1700     2      0.0
    C321   C321   C321   C2D1   0.0500     3      180.0
    C321   C321   C321   C2D1   0.1400     1      180.0
    HGA2   C321   C321   C331   0.1800     3      0.00
    C321   C321   C331   HGA3   0.1600     3      0.00
    HGA2   C321   C321   C2D1   0.1950     3      0.00
    C321   C321   C2D1   C2D1   0.6000     1      180.00
    C321   C321   C2D1   HGA4   0.1200     3      0.00
    HGA2   C321   C331   HGA3   0.1600     3      0.00
    C321   C2D1   C2D1   C321   8.5000     2      180.00
    C321   C2D1   C2D1   C321   0.4500     1      180.00
    HGA4   C2D1   C2D1   C321   1.0000     2      180.00
    C2D1   C2D1   C321   HGA2   0.3000     3      180.00
    HGA4   C2D1   C321   HGA2   0.0000     3      0.00
    HGA4   C2D1   C2D1   HGA4   1.0000     2      180.00
    
    IMPROPER
    
    NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
    cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
    C324   0.0000     -0.0550    2.1750
    N3P3   0.0000     -0.2000    1.8500
    HGP2   0.0000     -0.0460    0.2245
    C321   0.0000     -0.0560    2.0100
    HGA2   0.0000     -0.0350    1.3400
    C331   0.0000     -0.0780    2.0500
    C2D1   0.0000     -0.0680    2.0900
    HGA3   0.0000     -0.0240    1.3400
    HGA4   0.0000     -0.0310    1.2500

The input for our MD simulation, however, requires only **one** .prm file, so we will need to merge all of our .prm files in a single, global one. We will achieve this by manually copying and pasting the lines of each individual .prm file into a "global" one section by section. Pay attention to this step: the .prm file won't be read correctly if lines are missing or repeated twice. Take your time with this step and check twice to make sure everything has been pasted appropriately!
Now that our .prm and .rtf files are ready, we are _finally_ ready to proceed to the next step!
    
Preparing the box
-----------------


    - The Protein Structure File (.psf), containing the molecular-level information required to apply any force field to our simulation box;
