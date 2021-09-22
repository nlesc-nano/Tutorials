.. _simulation_box:

Preparing a Simulation Box
==========================

The goal of this tutorial is to prepare a 11.5 nm simulation box for Classical Molecular Dynamics (CMD) simulations. We have chosen to simulate the products of a synthesis procedure to obtain monodisperse CsPbBr\ :sub:`3`\  NCs, described by M. Imran *et al* in *J. Am. Chem. Soc.*, **2018**, *140*, 2656−2664.
The box specifically contains:

- One CsPbBr\ :sub:`3`\ core capped by 20% oleate (OA) and 20% oleylammonium (OLA) ligands;
- 75 ionic oleate-oleylammonium couples (OA+OLA, by-products obtained from the reaction);
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
First of all, let's start by building the NC capped with two different ligands (OA and OLA). We invite you to refer to the `tutorial <https://nanotutorials.readthedocs.io/en/latest/1_build_qd.html>`__ for the step-by-step construction of the structure from scratch.
Let's take a look at the .yaml input:

.. code:: yaml

    path: null
    input_cores:
        - core.xyz:
            guess_bonds: False
    
    input_ligands:
        - CCCCCCCCC=CCCCCCCCC(=O)O
    
    optional:
        core:
            dirname: core
            anchor: Cl
    
        ligand:
            dirname: ligand
            optimize: True
            split: True
    
        qd:
            dirname: qd
            construct_qd: True
            multi_ligand:
               ligands:
                 - CCCCCCCCC=CCCCCCCCC[NH3+]
               anchor:
                 - Rb
            optimize: False
            
The `path <https://cat.readthedocs.io/en/latest/2_path.html#path>`_, `input_cores & input_ligands <https://cat.readthedocs.io/en/latest/3_input_core_ligand.html#input-cores-input-ligands>`_ and  sections, together with the meaning of the `optional <https://cat.readthedocs.io/en/latest/4_optional.html#optional>`_ keywords and their relative `arguments <https://cat.readthedocs.io/en/latest/4_optional.html#arguments>`_, can be easily found inside the **CAT** `documentation <https://cat.readthedocs.io/en/latest/0_documentation.html#cat-documentation>`_. Let's take a look at some of them:

1. ``input_cores``: This section contains the coordinates of the core, specified by our .xyz file (``core.xyz``). The ``guess_bonds: False`` keyword tells **CAT** that, since our core is ionic, it does not need to guess the bonds and bond orders from the content of the .xyz file.
2. ``input_ligands``: This section contains information on both the structure and the chemistry of the ligand. This information is stored in its canonical `SMILES <https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Description>`_ (Simplified molecular-input line-entry system) string (``CCCCCCCC/C=C\CCCCCCCC(=O)[O-]`` for oleate, ``CCCCCCCC/C=C\CCCCCCCC[NH3+]`` for oleylammonium);
3. ``optional``: The `optional <https://cat.readthedocs.io/en/latest/4_optional.html#optional>`_ section contains three fairly similar subsections: ``core``, ``ligand``, ``qd``. The subsections contain keywords with several specifications, such as:

- the dummy atom that needs to be replaced with the chosen ligand (``optional.core.anchor``);
- whether or not to remove protons from the ligand (``optional.ligand.split``). Specifically, since the SMILES string we are using in the input (i.e. ``CCCCCCCC/C=C\CCCCCCCC(=O)[O-]``) refers to the anionic ligand, we will opt for ``optional.ligand.split: False``, so no protons have been removed from the ligand anchoring group. Conversely, if the SMILES is provided in the neutral form, then the ``optional.ligand.split: True`` key should be used, meaning that a proton will be cleaved from the functional group (in this case carboxylate) to ensure that the ligand is still added in its anionic form. Note that the latter form is preferrable when the ligand presents more than one functional group;
- the ``optional.qd.multiligand`` block. All the keys under this section are completely parallel to the aforementioned ones: Rb atoms are now being replaced by oleylammonium molecules. **Please note** that, in order to work effectively, this block acccepts SMILES strings by assuming a ``split: True`` specification.

An important concept to remember here, which we will need in a while, is that **CAT** builds the .xyz file in the following order: all the core atoms in the exact order we gave in the ``core.xyz``, followed by a certain number of ligand molecules (depending on the chosen coverage). If the model comprises more than one ligand, we will first have all of the molecules of the first ligand, followed by those of the second ligand. In our specific case, the order of our .xyz file will therefore be: Cs, Pb, Br, OA, OLA.
We are finally ready to run CAT with the following command: ``init_cat input_settings.yaml``.
After running **CAT** the .xyz file corresponding to our NC can be found in the specified directory, 'qd'. Don't worry, the directory will be created from scratch if it does not yet exist. Remember to rename the file before using it!

In a parallel fashion, the same script can be used to build the .xyz file containing OA+OLA molecules (i.e. our ionic oleate-oleylammonium couples) with two main differences: we will use a RbCl molecule as our "minimal", biatomic core, specified by our .xyz file (``RbCl.xyz``). In addition, we'll use the ``optional.core.allignment: sphere`` key, which is mandatory on **CAT** when diatomic molecules are set as cores in the script. The .xyz files of the remaining molecules (i.e. the .xyz files for ODA and OLAM) can be built using any (commonly available) molecular structure processing program, such as `Molden <https://www3.cmbi.umcn.nl/molden/>`__.
We will now have successfully built the following files (the names have been chosen to represent their chemical formula for simplicity):

- qd.xyz (our ligand-capped NC);
- oaola.xyz;
- olam.xyz;
- oda.xyz.

Other file extensions
^^^^^^^^^^^^^^^^^^^^^
.. only:: builder_html

Now that we've obtained our .xyz files, we need to convert them to other extensions to ensure our 3D structures can be read and used by the softwares while building our simulation boxes. Let's see the other extensions and how to obtain them:

1. *.pdb file*: The Protein Data Bank (.pdb) extension provides a description of the atomic coordinates, secondary structure assignments and atomic connectivity of our molecules. An .xyz file can be easily converted to this format by means of `Open Babel <https://openbabel.org/docs/dev/Installation/install.html>`__, a commonly employed chemical format converter. You can follow this link for the installation instructions (or just look for any Open Babel-based format converters available online). Once the program is correctly installed, we can convert our .xyz files to the .pdb format by running this simple command (note that this step does only apply to our organic molecules, i.e. **NOT** to our qd.xyz file): ``obabel -ixyz file.xyz -opdb file.pdb``.
We will now have the following .pdb files:

- oaola.pdb;
- olam.pdb;
- oda.pdb.
    
2. *.prm and .rtf files*: Each .pdb file we created now needs to be converted to the following formats:

- The CHARMM forcefield Parameter (.prm) file, including all of the numerical constants needed to evaluate forces and energies;
- The Residue Topology File (.rtf) This file defines the main groups (atoms, properties, bond and charge information) of our molecular structures.
    
These formats can be easily obtained from our .pdb files by inserting our .pdb files in `MATCH <https://brooks.chem.lsa.umich.edu/index.php?matchserver=submit>`__. This online server will convert our files into the three required formats, which we will download as a zipped directory (the one we obtained for OLAM can be found :download:`here </_files/3.olam.zip>`. We will first of all need to rename the new files to match their molecular formulas (2 for each .pdb file, for a total of 6 new files in this example).
Let's put the .rtf files aside and focus on the .prm files. An example of a MATCH-built .prm file (here, once again, we chose OLAM) looks like this:

::

    * prm file built by MATCH
    *
    
    BONDS
    C321   C321   222.50     1.5300     
    C321   HGA2   309.00     1.1110     
    C321   C331   222.50     1.5280     
    C321   N321   263.00     1.4740     
    C2D1   C321   365.00     1.5020     
    C331   HGA3   322.00     1.1110     
    HPA2   N321   453.10     1.0140     
    C2D1   C2D1   440.00     1.3400     
    C2D1   HGA4   360.50     1.1000     
    
    ANGLES
    C321   C321   C321   58.35      113.60     
    HGA2   C321   C321   26.50      110.10     
    HGA2   C321   HGA2   35.50      109.00     
    C321   C321   C331   58.00      115.00     
    C321   C321   N321   32.00      112.20     
    C321   C321   C2D1   32.00      112.20     
    HGA2   C321   C331   34.60      110.10     
    C321   C331   HGA3   34.60      110.10     
    HGA2   C321   N321   32.40      109.50     
    C321   N321   HPA2   41.00      112.10     
    C2D1   C2D1   C321   48.00      123.50     
    HGA4   C2D1   C321   40.00      116.00     
    C2D1   C321   HGA2   45.00      111.50     
    HGA3   C331   HGA3   35.50      108.40     
    HPA2   N321   HPA2   29.50      105.85     
    HGA4   C2D1   C2D1   52.00      119.50     
    
    DIHEDRALS
    C321   C321   C321   C321   0.14975    3      180.00     
    C321   C321   C321   C321   0.09458    4      0.00       
    C321   C321   C321   C321   0.11251    5      0.00       
    C321   C321   C321   C321   0.06450    2      0.00       
    HGA2   C321   C321   C321   0.1950     3      0.00       
    HGA2   C321   C321   HGA2   0.2200     3      0.00       
    C321   C321   C321   C331   0.08133    3      180.00     
    C321   C321   C321   C331   0.10824    4      0.00       
    C321   C321   C321   C331   0.20391    5      0.00       
    C321   C321   C321   C331   0.15051    2      0.00       
    C321   C321   C321   N321   0.1700     2      0.0        
    C321   C321   C321   N321   0.0500     3      180.0      
    C321   C321   C321   N321   0.1400     1      180.0      
    C321   C321   C321   C2D1   0.1700     2      0.0        
    C321   C321   C321   C2D1   0.0500     3      180.0      
    C321   C321   C321   C2D1   0.1400     1      180.0      
    HGA2   C321   C321   C331   0.1800     3      0.00       
    C321   C321   C331   HGA3   0.1600     3      0.00       
    HGA2   C321   C321   N321   0.1950     3      0.00       
    C321   C321   N321   HPA2   0.1600     3      0.00       
    HGA2   C321   C321   C2D1   0.1950     3      0.00       
    C321   C321   C2D1   C2D1   0.6000     1      180.00     
    C321   C321   C2D1   HGA4   0.1200     3      0.00       
    HGA2   C321   C331   HGA3   0.1600     3      0.00       
    HGA2   C321   N321   HPA2   0.0100     3      0.00       
    C321   C2D1   C2D1   C321   8.5000     2      180.00     
    C321   C2D1   C2D1   C321   0.4500     1      180.00     
    HGA4   C2D1   C2D1   C321   1.0000     2      180.00     
    C2D1   C2D1   C321   HGA2   0.3000     3      180.00     
    HGA4   C2D1   C321   HGA2   0.0000     3      0.00       
    HGA4   C2D1   C2D1   HGA4   1.0000     2      180.00     
    
    IMPROPER
    
    NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
    cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
    C321   0.0000     -0.0560    2.0100     
    HGA2   0.0000     -0.0350    1.3400     
    C331   0.0000     -0.0780    2.0500     
    N321   0.0000     -0.0600    1.9900     
    C2D1   0.0000     -0.0680    2.0900     
    HGA3   0.0000     -0.0240    1.3400     
    HPA2   0.0000     -0.0100    0.8750     
    HGA4   0.0000     -0.0310    1.2500     

The input for our Molecular Dynamics (MD) simulation requires only **one** .prm file, so we will need to merge all of our .prm files into one. We will achieve this by manually copying and pasting the lines of each individual .prm file into a "global" one section by section (BONDS, ANGLES, DIHEDRALS etc). Pay attention to this step: the .prm file won't be read correctly if lines are missing or repeated twice. Take your time with this step and check twice to make sure everything has been pasted appropriately!
Now that our .prm and .rtf files are ready, we are finally ready to proceed to the next step!
    
Preparing the box
-----------------
Once all of our .xyz files are ready, we need to build our final .xyz file by randomly inserting all of our molecules into a pre-shaped box. An useful tool for this purpose is provided by the **Packmol** package - again, the following `link <http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml>`__ provides all the information we need for its installation and proper usage. In order to build our box, we will first of all need to move all of our .xyz files into our working directory. For simplicity, let's assume that the packmol.exe executable is located in the same directory. The box will then be built by running a small script, ``settings.inp``, on the program. Let's take a brief look at our settings.inp file:

.. code:: yaml

    tolerance 2.0
    
    filetype xyz
    
    structure qd.xyz
      number 1
      inside cube -80. -80. -80. 80.
      center
      fixed 0. 0. 0. 0. 0. 0.
    end structure
    
    structure oaola.xyz
      number 75
      inside cube -80. -80. -80. 80.
    end structure
    
    structure olam.xyz
      number 287
      inside cube -80. -80. -80. 80.
    end structure
    
    structure octadecene.xyz
      number 2293
      inside cube -80. -80. -80. 80.
    end structure
    
    output box.xyz

The used keywords can be very easily found in the relative  `User Guide <http://leandro.iqm.unicamp.br/m3g/packmol/userguide.shtml>`__. Here is a very brief explanation:

- The line ``tolerance 2.0`` specifies the tolerance required for the distances between molecules. Here, the value has been set at 2.0 Å, a common value for systems at room temperature and pressure;
- The ``filetype xyz`` key specifies the formats of the provided molecular inputs;
- Individual blocks containing several specifications for the molecules which will figure in the box, such as their .xyz file and the number of molecules of each type that will be placed inside the box. In our case, as specified by the ``inside cube -80. -80. -80. 80.`` key, we will be placing the molecules inside a cube with minimum coordinates (x,y,z) = (-80,-80,-80) and maximum coordinates (80,80,80): in other words, we will fill a cube of side 160.0 Å with our molecules. We set the coordinates between -80 and 80 (instead of, for example, 0 to 160) because, as specified by the keywords ``center`` and ``fixed 0. 0. 0. 0. 0. 0.``, we want to place our NC model in the center of our box.

Once our input is ready, we can simply run the following command: ``packmol < settings.inp``.
Once the script has run, the ``box.xyz`` output containing the box will be inside the working directory. 

Generating the .psf file
------------------------
We will now need to build the Protein Structure File (.psf) of our simulation box, containing the molecular-level information required to apply the force field to our system over the course of our MD trajectory (here is an `example <https://github.com/nlesc-nano/Tutorials/build_qd/docs/_files/3.box_ordered.psf.zip>`__ of what ours looks like. You can take a look at this `website <https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html>`__ to get an idea of its structure). Let's now take a peek at its first lines:

::

    PSF EXT
    
             2 !NTITLE
       REMARKS PSF file generated with Auto-FOX
       REMARKS https://github.com/nlesc-nano/Auto-FOX
    
    
        153333 !NATOM
             1 MOL1     1        COR      Cs       Cs      0.000000      132.905450        0
             2 MOL1     1        COR      Cs       Cs      0.000000      132.905450        0
             3 MOL1     1        COR      Cs       Cs      0.000000      132.905450        0
             ..........

As mentioned in the website, each line in a .psf file is structured according to the following fields:

- atom ID (the number of the atom in the .xyz file);
- segment name (the number associated to each molecule: for us ``1`` is the whole NC, ``2`` is the **first** OA molecule, ``3`` is the second OA etc.)
- residue ID (in our case, ``MOL1`` to ``MOL3`` are the atoms of the NC core, ``MOL4`` is OA, ``MOL5`` is OLA, ``MOL6`` is OLAM and ``MOL7`` is ODA);
- residue name (COR specifically refers to our NC, while LIG is associated to ligand molecules);
- the remaining fields: atom name (e.g. C, H), atom type (e.g. C324, HGP2), charge, mass, and an unused 0.

The .psf file for our .xyz molecule can be easily built using the **Auto-FOX** package by means of a straightforward python script. Let's take a look at it:

.. code:: python

    from scm.plams import Molecule
    from FOX import PSFContainer
    from FOX.io.read_psf import overlay_rtf_file
    from FOX.recipes import generate_psf2
    
    qd = Molecule('box.xyz') 
    ligands = ('CCCCCCCCC=CCCCCCCCC(=O)[O-]', 'CCCCCCCCC=CCCCCCCCC[NH3+]', 'CCCCCCCCC=CCCCCCCCCN', 'CCCCCCCCCCCCCCCCC=C')
    psf = generate_psf2(qd, *ligands, ret_failed_lig=True)
    psf.write('box.psf')
    
    segment_dict = {"MOL4": Molecule('oa.xyz'), "MOL5": Molecule('ola.xyz'),  "MOL6": Molecule('olam.xyz'),  "MOL7": Molecule('oda.xyz')}
    psf_new, argsort = psf.sort_values(["segment name", "residue ID"], return_argsort=True)
    qd.atoms = [qd.atoms[i] for i in argsort]
    qd.write('box_ordered.xyz')
    
    for mol in segment_dict.values():
        mol.guess_bonds()
        
    psf_new.generate_bonds(segment_dict=segment_dict)
    psf_new.generate_angles(segment_dict=segment_dict)
    psf_new.generate_dihedrals(segment_dict=segment_dict)
    psf_new.generate_impropers(segment_dict=segment_dict)
    
    overlay_rtf_file(psf_new, 'oa.rtf', list(range(2, 129)))
    overlay_rtf_file(psf_new, 'ola.rtf', list(range(129, 245)))
    overlay_rtf_file(psf_new, 'olam.rtf', list(range(245, 532)))
    overlay_rtf_file(psf_new, 'oda.rtf', list(range(532, 2825)))
    
    psf_new.write('box_ordered.psf')

We'll now provide a step-by-step explanation of the purpose of the most important blocks in the script.

.. code:: python

    from scm.plams import Molecule
    from FOX import PSFContainer
    from FOX.io.read_psf import overlay_rtf_file
    from FOX.recipes import generate_psf2
    
    qd = Molecule('box.xyz') 
    ligands = ('CCCCCCCCC=CCCCCCCCC(=O)[O-]', 'CCCCCCCCC=CCCCCCCCC[NH3+]', 'CCCCCCCCC=CCCCCCCCCN', 'CCCCCCCCCCCCCCCCC=C')
    psf = generate_psf2(qd, *ligands, ret_failed_lig=True)
    psf.write('box.psf')
    
This section includes the generation of the .psf file in the order provided by our .xyz input. The ``generate_psf2`` key is motivated by the fact that our NC is capped by multiple ligands. You can find a very exhaustive documentation for this section in the `FOX.recipes.psf <https://auto-fox.readthedocs.io/en/latest/7_recipes.html?highlight=generate_psf#FOX.recipes.generate_psf2>`__ section of the relative documentation.

.. code:: python
    
    segment_dict = {"MOL4": Molecule('oa.xyz'), "MOL5": Molecule('ola.xyz'),  "MOL6": Molecule('olam.xyz'),  "MOL7": Molecule('oda.xyz')}
    psf_new, argsort = psf.sort_values(["segment name", "residue ID"], return_argsort=True)
    qd.atoms = [qd.atoms[i] for i in argsort]
    qd.write('box_ordered.xyz')

Before using our newly generated .psf file, we need to remember that the atoms/molecules in box.xyz have been packed by **Packmol** in the order specified by our input (settings.inp). As we've mentioned earlier, in our qd.xyz file this order is Cs, Pb, Br, OA, OLA. The residueIDs for the NC will thus be in ascending order (``MOL1`` to ``MOL5``) in the .psf file. In addition, each OA+OLA molecule has got an OA and an OLA in its .xyz file, so their lines in the .psf file will alternate between two residueIDs, ``MOL4`` and ``MOL5``. The file will then look like this:

::

      2959 MOL4     56       LIG      C        C   -0.180000       12.010600        0
      2960 MOL4     56       LIG      C        C   -0.180000       12.010600        0
      ....
      3011 MOL4     56       LIG      H        H    0.090000        1.007980        0
      3012 MOL5     57       LIG      N        N   -0.300000       14.006850        0
      3013 MOL5     57       LIG      C        C    0.210000       12.010600        0
      ....
      3066 MOL5     57       LIG      H        H    0.330000        1.007980        0
      3067 MOL4     58       LIG      C        C   -0.180000       12.010600        0
      3068 MOL4     58       LIG      C        C   -0.180000       12.010600        0

In order to build an ordered .psf file, we thus need to reorder our .xyz file so that all the molecules - as well as their residueIDs - are provided in ascending order. 
To do so, we created a dictionary (``segment_dict``) connecting every residueID in our box.psf file to the matching .xyz file. After that, we proceeded to reorder our .psf file by means of the ``sort_values`` key (you can find it in the `PSFContainer <https://auto-fox.readthedocs.io/en/latest/8_psf.html?highlight=sort_values#module-FOX.io.read_psf>`__ section). Specifically, the ``["segment name", "residue ID"]`` segment establishes that the molecules are ordered according to their residueIDs (``MOL4`` and **then** ``MOL5``):

::

      2959 MOL4     56       LIG      C        C   -0.180000       12.010600        0
      2960 MOL4     56       LIG      C        C   -0.180000       12.010600        0
      ....
      3011 MOL4     56       LIG      H        H    0.090000        1.007980        0
      3012 MOL4     58       LIG      C        C   -0.180000       12.010600        0
      3013 MOL4     58       LIG      C        C   -0.180000       12.010600        0
      ....
      3064 MOL4     58       LIG      H        H    0.090000        1.007980        0
      3065 MOL5     57       LIG      N        N   -0.300000       14.006850        0
      3066 MOL5     57       LIG      C        C    0.210000       12.010600        0
      ....
      3121 MOL5     57       LIG      H        H    0.330000        1.007980        0
      
**and** that, at the same time, their segment names are then reset to match this new order (``56``, ``57`` and **then** ``58``), as in:

::

      2959 MOL4     56       LIG      C        C   -0.180000       12.010600        0
      2960 MOL4     56       LIG      C        C   -0.180000       12.010600        0
      ....
      3011 MOL4     56       LIG      H        H    0.090000        1.007980        0
      3012 MOL4     57       LIG      C        C   -0.180000       12.010600        0
      3013 MOL4     57       LIG      C        C   -0.180000       12.010600        0
      ....
      3064 MOL4     57       LIG      H        H    0.090000        1.007980        0
      3065 MOL5     58       LIG      N        N   -0.300000       14.006850        0
      3066 MOL5     58       LIG      C        C    0.210000       12.010600        0
      ....
      3121 MOL5     58       LIG      H        H    0.330000        1.007980        0

we then proceeded to order the atoms in our box.xyz file (``qd.atoms``) in the same order of this .psf file, and we saved our new .xyz file as ``box_ordered.xyz``. Let's move on to the next section:

.. code:: python
    
    for mol in segment_dict.values():
        mol.guess_bonds()
        
    psf_new.generate_bonds(segment_dict=segment_dict)
    psf_new.generate_angles(segment_dict=segment_dict)
    psf_new.generate_dihedrals(segment_dict=segment_dict)
    psf_new.generate_impropers(segment_dict=segment_dict)
    
The contents of this section are pretty self-explanatory: the MultiMolecule `guess_bond <https://auto-fox.readthedocs.io/en/latest/3_multimolecule.html?highlight=guess_bonds#FOX.MultiMolecule.guess_bonds>`__ instance was used to guess the bonds in the file based on their atom types and inter-atomic distances. The bonds, angles, dihedrals and improper angles were then generated in the ordered .psf file for each residueID in ``segment_dict``.

.. code:: python

    overlay_rtf_file(psf_new, 'oa.rtf', list(range(2, 129)))
    overlay_rtf_file(psf_new, 'ola.rtf', list(range(129, 245)))
    overlay_rtf_file(psf_new, 'olam.rtf', list(range(245, 532)))
    overlay_rtf_file(psf_new, 'oda.rtf', list(range(532, 2825)))
    
    psf_new.write('box_ordered.psf')
    
We're almost there! This section of the script, which is specific for the organic molecules in our structure, matches each atom name in the .psf to its corresponding atom type (for example ``C321`` and ``N3P3``) which is specified in its .rtf file. The resulting .psf file, which finally looks like this:

::

      2959 MOL4     56       LIG      C        C321   -0.180000       12.010600        0
      2960 MOL4     56       LIG      C        C321   -0.180000       12.010600        0
      ....
      3011 MOL4     56       LIG      H        HGA2    0.090000        1.007980        0
      3012 MOL4     57       LIG      C        C321   -0.180000       12.010600        0
      3013 MOL4     57       LIG      C        C321   -0.180000       12.010600        0
      ....
      3064 MOL4     57       LIG      H        HGA2    0.090000        1.007980        0
      3065 MOL5     58       LIG      N        N3P3   -0.300000       14.006850        0
      3066 MOL5     58       LIG      C        C324    0.210000       12.010600        0
      ....
      3121 MOL5     58       LIG      H        HGP2    0.330000        1.007980        0

is then saved by means of the ``write`` method as ``box_ordered.psf``, and it is now ready to be used with our previously ordered .xyz file.

Preparing the simulation
------------------------

We have now got all the files we need to start our MD simulation. In our specific case, we will run the simulations on **GROMACS**, so we will need the `.gro <https://manual.gromacs.org/documentation/2018/user-guide/file-formats.html#gro>`__ file (for the starting molecular structure) and the topology file (`.top <https://manual.gromacs.org/documentation/2018/user-guide/file-formats.html#top>`__) of our box. As we mentioned earlier, we will use the **VMD** `software <https://www.ks.uiuc.edu/Research/vmd/>`__ package for this purpose.
First of all, we will open our .psf file on **VMD** (click on File > New Molecule in the Main Window and then Load the .psf file). Once the file is correctly loaded, we can proceed to load the .xyz structure in our .psf file by right clicking on the loaded .psf and selecting Load Data Into Molecule and our .xyz file). This procedure is common to both formats.
Let's now see how to obtain the two separate file extensions:

- *.gro file*: This file can be very easily obtained by selecting File > Save Coordinates > File Type: gro. The resulting `file <https://github.com/nlesc-nano/Tutorials/tree/build_qd/docs/_files/3.box_ordered.gro.zip>`__ (hereby provided) is now ready to be used.
- *.top file*: This `file <https://github.com/nlesc-nano/Tutorials/tree/build_qd/docs/_files/3.box_ordered.top.zip>`__ (you can find it by following the previous link) can be obtained from the **VMD** command line. We will first need to move to the directory containing our .prm file. After that, we can just insert the following commands in the terminal: ``topo writegmxtop box_ordered.top box.prm`` (``box.prm`` being our `previously built <https://nanotutorials.readthedocs.io/en/latest/3_simulation_box.html#other-file-extensions>`__ "global" .prm file). The .top file will be generated in the same directory with the name we specified in the command line. As our very last step before running the simulation, we will need to perform a few small modifications to the file:

1. The ``[ atomtypes ]`` section is to be updated to include the inorganic atoms (Cs, Pb, Br), as well as their relative parameters (atomic number, mass, charge etc.) in the description;
2. In the ``[ nonbond_params ]`` section each couple of atoms is associated to a sigma and an epsilon. In our case, these parameters account for the description of the Lennard-Jones terms in our force field, and we will need to insert their corresponding values in the column. The section would then look like this:
    
::
    
    [ nonbond_params ]
    ;type1 type2 1 sigma epsilon
    Pb Pb    1  0.6248523340799998  2.773992
    Br Pb    1  0.31212  1.7068104799678259
    ....
        
3. The charges in the ``[ moleculetype ]`` section, containing all the information on the atoms and molecules figuring in the structure, need to be updated as well to coincide to those in our force field. In our case we updated those of the inorganic core (Cs, Pb, Br) as well as those belonging to the anchoring groups of the ligands (C2O3, O2D2 etc). Here's a snippet of what the section should look like:
    
::
    
    [ moleculetype ]
    ; Name      nrexcl
    molecule0     3
    
    [ atoms ]
    ; nr  type  resnr residue atom cgnr charge  mass
      1   Cs    1     COR     Cs   1    0.6976  132.9055
    ....

We now have all of the files required to run our **GROMACS** simulation!

:download:`here  <https://github.com/nlesc-nano/Tutorials/tree/build_qd/docs/_files/3.olam.zip>`__
