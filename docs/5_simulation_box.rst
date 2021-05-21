.. _simulation_box:

Preparing a Simulation Box
==========================

The goal of this tutorial is to prepare a simulation box for Classical Molecular Dynamics (CMD) simulations. We have chosen to simulate the products of a synthesis procedure to obtain monodisperse CsPbBr\ :sub:`3`\ NCs, described by L. Protesescu *et al* in *Nano Letters*, **2015**, *15(6)*, 3692-3696.
The box specifically contains:
    * One CsPbBr\ :sub:`3`\ core capped with oleate (OA) and oleylammonium (OLA) ligands;
    * 2293 octadecene (ODA) molecules, used as solvent for the reaction;
    * 287 oleylamine (OLAM) molecules, used as a reagent in the synthesis;
    * 75 ionic oleate-oleylammonium couples (by-products obtained from the reaction).
    
Installation Requirements
-------------------------

This tutorial requires the download and use of the following programs:
    * the **CAT**, **data-CAT** and **nano-CAT** packages for the construction of the NC model and of the organic molecules in the box. The relative `documentation <https://cat.readthedocs.io/en/latest/0_documentation.html#cat-documentation>`__ is hereby provided for more information on the keywords. 
    * The **Packmol** package (see the following `link <http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml>`__ for insight on the installation and the keywords) for the construction of the coordinate file of the entire simulation box.
    * The **Auto-FOX** `package <https://auto-fox.readthedocs.io/en/latest/includeme.html>`__, required for the construction of the Protein Structure File (.psf) containing the molecular-level information required to apply any force field to our simulation box;
    * Ultimately, the **VMD** `software <https://www.ks.uiuc.edu/Research/vmd/>`__ has been employed as a molecular visualization program for the construction of the topology (.top) and the Gromacs (.gro) files to run the MD simulations.
