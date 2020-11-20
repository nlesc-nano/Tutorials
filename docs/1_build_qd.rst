.. _build_qd:

Build a Quantum Dot Model
=========================
Quantum Dots (QDs) are colloidal semiconductor nanocrystals, usually spanning 2-10 nm in diameter. Their optoelectronic properties arise from surface-dependent quantum effects, and are thus dependent on their size and shape. 

These semiconductors are characterized by an intrinsically low photoluminescence quantum yield (PLQY) due to the presence of midgap states known as *surface traps*. A commonly used approach to optimize the activity of these materials consists in covering the photoactive inorganic core with a shell of a wider band gap material, usually an organic ligand.

The properties of the Quantum Dot will then depend on the surface coverage of the ligand on the core and on the degree at which the traps are "filled".
Due to these premises, QDs are built by capping an appropriately designed crystalline core with a chosen organic ligand. 

The goal of this tutorial is to outline the steps to build a Quantum Dot from scratch. In this tutorial we will build a 4.2 nm sided cubic perovskite CsPbBr_3 QD capped by 80% of oleate molecules.

The inorganic core
---------------
For starters, we need to download the Crystallographic Information File (CIF) of the inorganic species we want to use as a framework (CsPbBr_3 in our case) for our QD. The CIF file provides a precise numerical description of the crystallographic structure, and it can be downloaded from several different databases and libraries.

Once the file is updated, it needs to be uploaded in a visualization program. Visualization programs usually show the unit cell of the crystal. This cell will be the starting point to build a supercell of the appropriate dimension. We generated a 8x8x8 cell from a cubic CsPbBr_3 unit cell, but the choice of the supercell dimension is usually a compromise between the computational cost of the follow-up calculations and the necessity to provide a realistic description of the QD of interest. 

In our specific case (i.e. cubic CsPbBr_3 capped by oleate molecules), the available experimental data from the synthesis of the QDs shows that the crystalline framework terminates with the Cs-Br layer. We thus manually deleted the external Pb-Br layers from the supercell so that its ending layers were the Cs-Br ones.
This is a fairly common procedure used to adapt the crystal framework to mimick the experimentally obtained inorganic cores.

It is now necessary to ensure that the newly built supercell is neutral. Calculating the charge of a supercell is fairly easy, since it can be done by counting its atoms and summing their charges. In our CsPbBr_3 supercell, for example, we used our visualization program to count:

- 512 Cs atoms, each carrying charge +1 in the crystal;
- 343 Pb atoms, each carrying charge +2 in the crystal;
- 1176 Br atoms, each carrying charge -1 in the crystal;

The charge of the supercell can then be obtained as:

512x(+1) + 343(+2) + 1176(-1) = 512 + 686 - 1176 = +12

The supercell therefore has an excess of 12 cations in the structure. Neutral charge can then be achieved by manually removing any charged atoms (cations in our specific case) from the surface of the inorganic core. 

The neutral charge requirement ensures that the QD is effectively stoichiometric, so that the band gap is clean and free of midgap states in principle. It is known that the removal of atoms on the corners and edges of the QD minimizes the distortion associated to the crystalline framework. Those atoms (Cs in our specific case) should then be deleted one by one until the charge of the supercell has been balanced. If possible, the removal of the atoms should be done in a symmetrical fashion (e.g. opposite corners, atoms at the same "coordinate" of the edge). In our case we need to remove 12 atoms, so we removed 8 from the corners of the supercell "cube" and 4 from the edges.
Once the core is neutral we are ready to save and export the resulting cartesian coordinates of the atoms in the supercell to an .xyz file.

Using dummies
---------------
Capping our inorganic scaffold with organic ligands means replacing a certain amount of its superficial ions with appropriately charged ligands (in our case, this would mean replacing Br- anions with oleate anions). Since neutrality has been mentioned to be a requirement, the replacement has to ensure that the total charge of the final QD is still neutral (i.e. if a certain number *n* of Br- anions are removed, *n* oleate anions need to be added). Moreover, the capping procedure depends on the surface coverage we want our Quantum Dot to have (in our case, we chose to cover 80% of the surface of our perovksite core).

All of these requirements can be fulfilled by an intermediate step: the initial replacement of the inorganic ion (Br- in our specific example) with a dummy ion of the same charge. We chose to use Cl- as our dummy ion. The replacement can be done by means of a small python script.
Please note that the script requires a **CAT** module. We invite you to read the relative `documentation <https://cat.readthedocs.io/en/latest/0_documentation.html#cat-documentation>`_ before continuing this tutorial.
Let's now have a look at the script:

.. code:: python

    >>> from scm.plams import Molecule
    >>> from CAT.recipes import replace_surface
    >>> mol = Molecule('cspbbr3_4.2nm.xyz')
    >>> mol_new = replace_surface(mol, symbol='Br', symbol_new='Cl', f=0.8, mode='uniform', displacement_factor=0.7)
    >>> mol_new.write('cspbbr3_4.2nm_80Cl.xyz')
    
The script is pretty self-explanatory: the .xyz coordinates are imported from our previously-built file (``'cspbbr3_4.2nm.xyz'``). A specifically built recipe, ``replace_surface``, is then able to recognize **only** the requested atoms on the surface of the supercell (``symbol='Br'``) from their chemical symbol in the .xyz file and replace them with dummies (``symbol_new='Cl'``) in a new .xyz file (``'cspbbr3_4.2nm_80Cl.xyz'``). We specifically chose to replace Br atoms with Cl atoms because we aim to use oleate as a ligand, but this setup can be varied if needed.

The file also specifies the `fraction <https://cat.readthedocs.io/en/latest/4_optional.html#optional.core.subset.f>`_ of the atoms that are being replaced (80% in our case, hence the ``f=0.8`` in the script), their `distribution <https://cat.readthedocs.io/en/latest/4_optional.html#optional.core.subset.mode>`_ and the displacement factor resulting from placing the new dummies on the surface of the core (i.e. 0.7).

CAT input: building the Quantum Dot
---------------
We are now ready to use **CAT** to build our Quantum Dot. We will first of all need to build our 'core' and 'ligand' directories inside our working directory (see the `General Overview <https://cat.readthedocs.io/en/latest/1_get_started.html#default-settings>`_ for further information).
Therefore, our newly built .xyz file needs to be moved into the 'core' directory.

We will then need to write the .yaml `input file <https://cat.readthedocs.io/en/latest/includeme.html#input-files>`_,containing all the desired settings, and to put it in the working directory.

Let's take a look at the keywords required for our .yaml file:

.. code:: yaml

    path: null

    input_cores:
        - cspbbr3_4.2nm_80Cl.xyz:
            guess_bonds: False

    input_ligands:
        - CCCCCCCCC=CCCCCCCCC(=O)[O-]

    optional:
        database:
            dirname: database
            read: False
            write: True
            overwrite: False
            thread_safe: False
            mol_format: xyz
            mongodb: False

        core:
            dirname: core
            anchor: Cl
            subset: null

        ligand:
            dirname: ligand
            optimize: True
            split: False
            anchor: null
            cosmo-rs: False

        qd:
            dirname: qd
            construct_qd: True
            optimize: False
            bulkiness: False
            activation_strain: False
            dissociate: False
        
The `path <https://cat.readthedocs.io/en/latest/2_path.html#path>`_, `input_cores & input_ligands <https://cat.readthedocs.io/en/latest/2_path.html#path>`_ and  sections, together with the meaning of the `optional <https://cat.readthedocs.io/en/latest/4_optional.html#optional>`_ keywords and their relative `arguments <https://cat.readthedocs.io/en/latest/4_optional.html#arguments>`_, can be easily found inside the **CAT** `documentation <https://cat.readthedocs.io/en/latest/0_documentation.html#cat-documentation>`_.

The sections are all fairly similar: their keywords contain several specifications, such as:

- the directories where our inorganic cores and ligands will be stored (``optional.database.dirname``);
- whether or not their optimization is required (``optional.ligand.optimize`` and ``optional.core.optimize``);
- the dummy atom that needs to be replaced with the chosen ligand (``optional.core.anchor``)

Please note that, in this specific case, we chose to opt for ``optional.ligand.split: False`` since the SMILES string we are using in the input (i.e. ``CCCCCCCCC=CCCCCCCCC(=O)[O-]``) refers to an ionic structure.

Now that all of the files are in their respective directories, we are finally ready to run CAT with the following command: ``init_cat input_settings.yaml``
After running **CAT** the new .xyz file, containing the coordinates of the desired Quantum Dot, will be exported to the directory we specified in ``optional.qd.dirname`` ( we named it 'qd'). Don't worry, the directory will be created from scratch if it does not yet exist!

Make sure to rename the .xyz file so that you know what it is, as its name is randomly generated by **CAT**. Once renamed, the .xyz file is ready to be used.
