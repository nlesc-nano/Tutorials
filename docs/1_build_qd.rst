.. _build_qd:

Build a Quantum Dot Model
=========================
QDs are colloidal semiconductor nanocrystals, usually spanning 2-10 nm in diameter. Their optoelectronic properties arise from surface-dependent quantum effects, and are thus dependent on their size and shape. 

These semiconductors are characterized by an intrinsically low photoluminescence quantum yield (PLQY) due to the presence of midgap states known as *surface traps*. A commonly used approach to optimize the activity of these materials consists in covering the photoactive inorganic core with a shell of a wider band gap material, usually an organic ligand.

The properties of the Quantum Dot will then depend on the surface coverage of the ligand on the core and on the degree at which the traps are "filled".
Due to these premises, QDs are built by capping an appropriately designed crystalline core with a chosen organic ligand. 

The goal of this tutorial is to outline the steps to build a Quantum Dot (QD) from scratch. In this tutorial we will build a 4.2 nm sided cubic CsPbBr_3 NC capped by 80% of oleate molecules.

The inorganic core
---------------
For starters, we need to build the inorganic framework from a precise numerical description of the crystallographic structure we want to use as a core (CsPbBr_3 in our case). This kind of structural information is stored in the Crystallographic Information File (CIF), which can be downloaded from several different databases and libraries.

Once the file is updated, it needs to be uploaded in a visualization program. The unit cell will be used as a starting point to build a supercell of the appropriate dimension - the choice is usually based on a compromise between the computational cost of the follow-up calculations and the necessity to provide a realistic description of the QD of interest. If the layers at the "sides" and the "corners" of the supercell contain less sterically demanding atoms (Pb in our specific example) it is common practice to manually delete them, so that the ending layers contain the bigger atoms (Cs in our case).

It is now necessary to ensure that the newly built supercell is neutral. This requirement ensures that the QD is effectively stoichiometric, so that the band gap is clean and free of midgap states in principle. Neutral charge can be achieved by manually removing any charged atoms from the surface of the inorganic core. It is known that the removal of atoms on the corners and edges of the QD minimizes the distortion associated to the crystalline framework. The atoms should then be deleted one by one until the charge of the supercell has been balanced. If possible, the removal of the atoms should be done in a symmetrical fashion (e.g. opposite corners, atoms at the same "coordinate" of the edge).

We are now ready to export the resulting cartesian coordinates of the atoms in the cell in an .xyz file.

Using dummies
---------------
Capping our inorganic scaffold with organic ligands means replacing a certain amount of its superficial ions with appropriately charged ligands (in our case, this would mean replacing Br- anions with oleate anions). Since neutrality has been mentioned to be a requirement, the replacement has to ensure that the total charge of the final QD is still neutral (i.e. if a certain number *n* of Br- anions are removed, *n* oleates need to be added). Moreover, a certain surface coverage needs to be ensured.

All of these requirements can be fulfilled by an intermediate step: the initial replacement of the inorganic ion with a dummy ion of the same charge. The replacement can be done by means of a small python script.
Please note that the script requires a **CAT** module. We invite you to read the relative `documentation <https://cat.readthedocs.io/en/latest/0_documentation.html#cat-documentation>`_ before continuing this tutorial.
Let's now have a look at the script:

.. code:: python

    >>> from scm.plams import Molecule
    >>> from CAT.recipes import replace_surface
    >>> mol = Molecule('cspbbr3_4.2nm.xyz')
    >>> mol_new = replace_surface(mol, symbol='Br', symbol_new='Cl', f=0.8, mode='uniform', displacement_factor=0.7)
    >>> mol_new.write('cspbbr3_4.2nm_80Cl.xyz')
    
The script is pretty self-explanatory: the .xyz coordinates are imported from our previously-built file (``'cspbbr3_4.2nm.xyz'``). A specifically built recipe, ``replace_surface``, is then able to recognize **only** the atoms on the surface of the supercell from their chemical symbol in the .xyz file and replace them with dummies in a new .xyz file (``'cspbbr3_4.2nm_80Cl.xyz'``). We specifically chose to replace Br atoms with Cl atoms because we aim to use oleate as a ligand, but this setup can be varied if needed.

The file also specifies the `fraction <https://cat.readthedocs.io/en/latest/4_optional.html#optional.core.subset.f>`_ of the atoms that are being replaced (80% in our case, hence the ``f=0.8`` in the script), their `distribution <https://cat.readthedocs.io/en/latest/4_optional.html#optional.core.subset.mode>`_ and the displacement factor resulting from the new dummies (i.e. 0.7).

CAT input: building the Quantum Dot
---------------
We are now ready to use **CAT** to build our Quantum Dot. We will need to build our 'core' and 'ligand' directories inside our working directory (see the `General Overview <https://cat.readthedocs.io/en/latest/1_get_started.html#default-settings>`_ for further information). Our newly built .xyz file needs then to be moved into the 'core' directory.
Let's now take a look at the .yaml file containing our settings:

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
        
The contents of the ``path``.....




After running **CAT** the .xyz file containing the QD will be exported to a newly built directory, named 'qd'. Once renamed, the .xyz file is ready to be used.
