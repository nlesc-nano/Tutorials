.. _Dataset_and_Preparation:

Machine Learning Force Fields (MLFF): Comprehensive Dataset Generation Guide
============================================================================

This tutorial provides a comprehensive guide to generating a robust and diverse dataset for training **Machine Learning Force Fields (MLFF)** for **Quantum Dots (QDs)**. The workflow integrates several computational techniques to ensure extensive coverage of the QD system’s configurational and chemical space.

Workflow Overview
-----------------

- **Ab-initio Molecular Dynamics (AIMD) Simulations (DFT-based)**  
  The process begins with AIMD simulations based on **Density Functional Theory (DFT)**. This step generates realistic atomic configurations and force data by accurately modeling the **QD system's dynamic behavior** at the atomic scale. The AIMD simulations are run using the **CP2K** quantum chemistry package.

- **Enhanced Sampling via Principal Component Analysis (PCA)**  
  To **broaden** the configurational space explored by AIMD, **PCA** is applied to the molecular dynamics trajectory. This statistical method identifies dominant modes of structural variation, enabling the generation of new, **diverse configurations** that capture essential system dynamics.

- **High-Accuracy DFT Calculations on New Samples**  
  The newly generated configurations are further refined through **high-precision DFT calculations** using the **QMflows** package. **QMflows** is a library that automates the input generation for CP2K-based DFT calculations. This step **computes accurate energy and force data**, enriching the dataset for effective MLFF training.

- **Dataset Preparation for Machine Learning Models**  
  The final step involves organizing the collected data into a **machine-learning-friendly format** for further processing.

By systematically combining **AIMD simulations, PCA-enhanced sampling, and high-accuracy DFT calculations**, this workflow ensures the development of a high-quality dataset.

Step 1: Running Ab-initio Molecular Dynamics (AIMD) Simulation
--------------------------------------------------------------

Objective
~~~~~~~~~

Obtain initial atomic configurations and force data necessary for **MLFF development** by simulating realistic **QD system behavior**.

Simulation Setup
~~~~~~~~~~~~~~~~

- Start by **preparing a QD model** (see relevant tutorial) and **relax its geometry**.  
  - Typically, a **2–3 nm QD passivated with Cl atoms** to charge balance the system is sufficient.
- Perform an **AIMD simulation** using an **NVT ensemble** at **300 K** (or another desired temperature).  
- Set the simulation duration to **5–10 picoseconds** to explore the system's **potential energy landscape**.  
- Extract **2000 structural snapshots** uniformly throughout the simulation to capture diverse atomic configurations.  

CP2K Input Configuration
~~~~~~~~~~~~~~~~~~~~~~~~

Use the following **CP2K input configuration** in the ``&MOTION`` block to enable printing for each frame. Adjust the **printing frequency** according to your needs.

.. code-block:: bash

   &MOTION
     &PRINT
       &TRAJECTORY
         FORMAT XYZ
         UNIT angstrom
         &EACH
           MD 1
         &END EACH
       &END TRAJECTORY
       &FORCES
         &EACH
           MD 1
         &END EACH
       &END FORCES
     &END PRINT
   &END MOTION

Considerations for Quantum Dots (QDs)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- For **quantum dots like CdSe passivated with halogens**, atomic motions are slower.  
  - Use a **timestep of 2–4 fs** to effectively explore the configurational space.  
- **Equilibrate** the system with **1000 NVT steps** before running a **production simulation** of **2000 frames**.  
  - Example: An **8 ps** simulation with a **4 fs** timestep.  
- **Record positions and forces at every step** to ensure a **detailed dataset** for subsequent analysis.

Step 2: Enhancing AIMD Data with Principal Component Analysis (PCA)
-------------------------------------------------------------------

**Objective:**  
Broaden the sampled configuration space by applying PCA to the AIMD data and generating additional, diverse structures.

Principal Component Analysis (PCA) in Dataset Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Principal Component Analysis (PCA)** is a statistical technique used to reduce the dimensionality of data by transforming it into a new set of variables called **principal components**.  

In the context of dataset generation for **MLFF**, PCA helps identify the most significant variations in atomic configurations by analyzing:
- **Energies**
- **Atomic positions**
- **Forces**
- **RMSD (Root-Mean-Square Deviation)**
- **SOAP (Smooth Overlap of Atomic Positions) descriptors**

By projecting the data onto the principal components, PCA effectively reveals directions of **maximum variance**, enabling the creation of **diverse and representative configurations** that capture essential structural dynamics. This approach ensures that the generated dataset spans the most relevant regions of the chemical space.

You can download the script `generate_mlff_dataset.py` from:  
`https://github.com/nlesc-nano/MLFF_QD/tree/main/src/mlff_qd/preprocessing`  

Then run the script with the following command:

.. code-block:: bash

   generate_mlff_dataset.py input.yaml

The YAML input processes the AIMD trajectory by reading positions and forces. An example YAML configuration file is provided below.

Example YAML Configuration for the Script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml 

   pos_file: "mean_md-pos-1.xyz"
   frc_file: "mean_md-frc-1.xyz"
   scaling_factor: 0.4
   scaling_surf: 0.6
   scaling_core: 0.4
  max_random_displacement: 0.15
  surface_atom_types:
    - "In"
    - "P"
    - "Cl"
  clustering_method: "KMeans"
  num_clusters: 100
  num_samples_pca: 1200
  num_samples_pca_surface: 600
  num_samples_randomization: 200
  SOAP:
    species: ["In", "P", "Cl"]
    r_cut: 12.0
    n_max: 7
    l_max: 3
    sigma: 0.1
    periodic: False
    sparse: False

Structure Generation Breakdown
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **1200 Structures from PCA Sampling**  
  These structures are generated by perturbing configurations along the **principal components** derived from the **entire atomic system**.  
  - This enhances the dataset by exploring high-variance directions in the molecular dynamics trajectory.

- **600 Structures with Surface-Specific PCA Sampling**  
  Here, PCA is applied **specifically to surface atoms** (e.g., **Cs** and **Br** in QDs), which are **more dynamic** than core atoms.  
  - This approach ensures the **surface chemistry** is well-represented by applying **larger displacements** to surface atoms, reflecting their **natural mobility**.

- **200 Structures from Random Sampling**  
  Random displacements are applied **uniformly** across selected structures to introduce **additional diversity** and help avoid biases in the sampled configurations.

Detailed Explanation of YAML Input Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **pos_file**  
  Path to the `.xyz` file containing **atomic positions** from the AIMD simulation.  
  - This file serves as the input for PCA analysis.

- **frc_file**  
  Path to the `.xyz` file containing **corresponding atomic forces**.  
  - These forces are used to evaluate **structural dynamics**.

- **max_random_displacement**  
  The **maximum displacement** applied in the **random sampling step**.  

- **surface_atom_types**  
  A list of **atomic species** (e.g., `"In"`, `"Cl"`) considered as **surface atoms**.  
  - These atoms are **more prone to movement** and are treated differently during **PCA sampling**.

- **clustering_method**  
  The algorithm used for **clustering structures** in the **PCA space**.  
  - Here, `KMeans` is used to **group similar configurations** and **sample representative ones**.

- **num_clusters**  
  The **number of clusters** to create in **PCA space** for **diversity sampling**.  
  - Each cluster provides **representative structures**.

- **num_samples_pca**  
  Number of structures generated by **perturbing configurations along PCA components** applied to the **entire system (core + surface)**.

- **num_samples_pca_surface**  
  Number of structures generated by applying **PCA perturbations** specifically to **surface atoms**, allowing them **greater freedom to move**.

- **num_samples_randomization**  
  Number of **randomly perturbed structures** added to the dataset to increase **diversity**.

**SOAP** refers to **Smooth Overlap of Atomic Positions**:

- **species**: adjust according to your model.
- **r_cut**: a cutoff for the neighbouring environment.
- **n_max**: max number of radial basis functions (RBF).
- **l_max**: max degree of spherical harmonics.
- **sigma**: the width of smearing.


Output Files and Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Generated Structures:**  
  The script outputs a `dataset_2000.xyz` file containing:
  - 1200 PCA-sampled structures
  - 600 surface-PCA structures
  - 200 randomized structures

- **PCA Plots:**  
  Visualizations illustrate the distribution of the sampled structures in PCA space, providing insights into the configurational diversity achieved compared to the reference AIMD trajectory.

By combining **PCA-driven sampling, surface-specific perturbations, and randomization**, this approach ensures a well-balanced dataset that thoroughly explores the system's **chemical and configurational space**.

Step 3: High-Accuracy DFT Calculations on the Generated Structures
------------------------------------------------------------------

In this step, the **2000 structures** generated in the previous step using the **enhanced sampling process** will be computed at the **Density Functional Theory (DFT)** level of theory.  
This process enables the calculation of **energy** and **force** data for these new configurations, which will be **added to the starting AIMD dataset**.

Detailed Workflow
~~~~~~~~~~~~~~~~~

1. **Organize the Working Directory**

   - Create a new folder to run the **DFT calculations**.  
   - Copy the file `dataset_2000.xyz` (which contains the **2000 sampled structures**) into this folder.  
   - Copy the `train.yaml` configuration file into the same directory.  
     - This file will guide the **DFT calculation setup**.

   **Example Commands:**

   .. code-block:: bash

      mkdir DFT_Calculations
      cp dataset_2000.xyz DFT_Calculations/
      cp train.yaml DFT_Calculations/
      cd DFT_Calculations/

2. **Configure the `train.yaml` Input File**

   Open the `train.yaml` file and adjust the settings according to your **computational environment**.

Example YAML Configuration for the Script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `train.yaml` file contains various parameters necessary for **DFT calculations**, including:

- **DFT functional** (e.g., PBE, B3LYP)  
- **Basis set**  
- **Convergence criteria**  
- **HPC-specific settings** (e.g., number of cores, memory allocation)  

**Example Configuration:**

.. code-block:: yaml

   workflow:
       distribute_single_points

   project_name: PbSe_Cl
   calculate_guesses: "all"
   active_space: [100, 100]
   path_traj_xyz: “dataset_2000.xyz"
   path_hdf5: "CdSe_Cl.hdf5"
   scratch_path: "cp2k_chunks"
   workdir: "."
   blocks: 5

   job_scheduler:
       free_format: "
           #!/bin/bash \n
           #SBATCH --job-name=PbSe_cl_single_point_cal \n
           #SBATCH --time=24:00:00 \n
           #SBATCH --nodes 2 \n
           #SBATCH --ntasks-per-node=112 \n
           module load cp2k/2024.1\n"

   cp2k_general_settings:
       path_basis: “cp2k_basis"
       basis_file_name: "BASIS_MOLOPT"
       potential_file_name: "GTH_POTENTIALS"
       basis: "DZVP-MOLOPT-SR-GTH"
       potential: "GTH-PBE"
       cell_parameters: 49.0
       periodic: none
       executable: cp2k.popt
       wfn_restart_file_name: "scf.wfn”

   cp2k_settings_main:
       specific:
           template: "train_main"  

   cp2k_settings_guess:
       specific:
           template: "train_guess"  

Detailed Explanation of YAML Input Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **workflow:**  
  Defines the overall workflow. Here, ``"distribute_single_points"`` is used for **single-point energy and force calculations**.

- **project_name:**  
  Specifies the **name of the project**, e.g., ``"PbSe_Cl"``, which will be used for organizing output files.

- **calculate_guesses:**  
  Determines whether **initial wavefunction guesses** should be computed for each frame.  
  - ``"all"`` means **all frames** will undergo:
  
    1. **Orbital Transformation (OT) calculations** to obtain an efficient initial guess.
    2. A **main calculation**, which then **fully diagonalizes the Fock matrix**.

- **active_space:**  
  Defines the **number of active molecular orbitals** whose **coefficients and energies** will be stored in the **HDF5 file**.

- **path_traj_xyz:**  
  Path to the ``dataset_2000.xyz`` file containing **generated atomic structures**.  
  - This file **stacks all generated frames**, which will be computed using **DFT**.

- **path_hdf5:**  
  Path to an existing **HDF5 database**, which stores **DFT-derived properties**, including:

  - **Atomic structures (XYZ format)**
  - **Molecular Orbital (MO) coefficients**
  - **Other relevant electronic structure data**

- **scratch_path:**  
  Specifies the **temporary directory** where **DFT calculations** will be executed.

- **workdir:**  
  Specifies the **working directory**, where all **calculation results** will be stored.

- **blocks:**  
  Defines the **number of blocks** into which the **original dataset** will be split.  
  - This helps manage **computational efficiency** when running large-scale DFT calculations.

- **job_scheduler:**  
  Contains **HPC job submission settings**, including:

  .. code-block:: bash

     #SBATCH --job-name=PbSe_cl_single_point_cal
     #SBATCH --time=24:00:00
     #SBATCH --nodes=2
     #SBATCH --ntasks-per-node=112
     module load cp2k/2024.1

  - ``#SBATCH --job-name``: Specifies the **job name**.
  - ``#SBATCH --time``: Maximum **runtime allocation**.
  - ``#SBATCH --nodes``: Number of **compute nodes** requested.
  - ``#SBATCH --ntasks-per-node``: Number of **tasks per node**.
  - ``module load cp2k/2024.1``: Loads the **CP2K module** on the **HPC system**.

- **cp2k_general_settings:**  
  Contains **general CP2K input settings**, including:

  - ``basis_file_name``: Specifies the **MOLOPT basis set**.
  - ``potential_file_name``: Defines the **GTH pseudopotentials**.
  - ``cell_parameters``: Defines the **simulation box size** (e.g., **49.0 Å**).
  - ``periodic``: Specifies **boundary conditions** (``none`` for **isolated QDs**).
  - ``executable``: Points to the **CP2K binary** (e.g., ``cp2k.popt``).

- **cp2k_settings_main:**  
  Specifies the **main CP2K input template**, based on the **PBE functional**.  
  - Example: ``"train_pbe_main"``.

- **cp2k_settings_guess:**  
  Specifies the **wavefunction guess template**, also based on the **PBE functional**.  
  - Example: ``"train_pbe_guess"``.


3. **Run the QMflows Job Distribution Script**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   Use **`qmflows`** and **`nano-qmflows`** to automate the **DFT calculations**. Assuming both are already installed, launch the job distribution script:

   .. code-block:: bash

      distribute_jobs.py -i train.yaml

   **Process Explanation:**  
   - The script **splits the dataset** into **5 folders** (or more/less depending on settings) to parallelize calculations.  
   - In each folder, it generates the necessary **input files**, **Slurm job scripts**, and setup for the **DFT calculations**.  
   - The `input.yaml` generated is a **pre-processed YAML file** containing **all keywords, including default ones**, that will be used to generate the **CP2K input files**.  
   - It is always recommended to **check this file** before running the job to ensure all settings are correct.

4. **Submit the Jobs to the HPC Cluster**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   Navigate into each **generated folder** and submit the job to the HPC queue:

   .. code-block:: bash

      cd chunk_1/
      sbatch lauch.sh
      cd ../chunk_2/
      sbatch launch.sh
      # Repeat for all chunks

   *Tip:*  
   - You can **increase the number of chunks** to reduce the computational load per job and **speed up the calculations**, depending on the available **HPC resources**.

5. **Handling Interrupted Jobs**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   If any calculation is **interrupted** (e.g., due to **wall-time limits**), simply **rerun** the job distribution script.  
   **QMflows** efficiently manages restarts, ensuring **only incomplete calculations** are resumed:

   .. code-block:: bash

      sbatch launch.sh

Key Points to Consider
^^^^^^^^^^^^^^^^^^^^^^

- **Parallelization:** Adjust the **number of chunks** for optimal performance on your HPC system. More chunks with fewer structures can speed up computations.  
- **Resource Management:** Customize the Slurm scripts (`job.sh`) as needed for your **HPC environment**.  
- **Automatic Restart:** **QMflows** handles restarts smoothly, allowing you to **resume incomplete jobs** without manual intervention.

----

Step 4: Convert All DFT Structures to ML-Ready Format
-----------------------------------------------------

1. **Extract Structures from Chunk Folders**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   Once the **DFT calculations** in each chunk are completed, download the `extract.py` script and run it inside each folder:

   .. code-block:: bash

      python extract.py

   **Process Explanation:**  
   - The script scans the folder for **output files** containing **positions** and **forces**.  
   - It **identifies redundant structures** if some calculations have **failed** and required **restarts**.  
   - The most relevant output files are:

     * `positions_hartree_n.xyz`
     * `forces_hartree_n.xyz`  

     where `n` is the chunk number.

   **Merge all chunk outputs into a single file** (assuming **5 chunks** were generated):

   .. code-block:: bash

      cat positions_hartree_0.xyz positions_hartree_1.xyz positions_hartree_2.xyz \
          positions_hartree_3.xyz positions_hartree_4.xyz > positions_hartree_final.xyz 

      cat forces_hartree_0.xyz forces_hartree_1.xyz forces_hartree_2.xyz \
          forces_hartree_3.xyz forces_hartree_4.xyz > forces_hartree_final.xyz

2. **Merge Extracted Structures with MD Structures**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   Now, merge these **DFT-calculated structures** with those obtained from the **initial MD simulation**:

   .. code-block:: bash

      cat mean_md-pos-1.xyz positions_hartree_final.xyz > merged_positions.xyz 
      cat mean_md-frc-1.xyz forces_hartree_final.xyz > merged_forces.xyz

3. **Convert All DFT Structures to ML-Ready Format**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   Use the `compact_xyz.py` script to generate a **single XYZ file** ready for **ML training**.  
   The script ensures:
   - **All frames computed with DFT** are included.
   - The file contains **energies, positions, and forces**.
   - Unit conversion is performed:

     * **Energies** → **eV**
     * **Positions** → **Ångström** (already in this format)
     * **Forces** → **eV/Ångström** (preferred for ML training)

   Run the script with:

   .. code-block:: bash

      python compact_input.py --pos merged_positions.xyz --frc merged_forces.xy

  Use `consolidate.py` to pick random structures suitable for ML training:

    .. code-block:: bash

        python consolidate.py input.yaml

  An example of input YAML file:

   .. code-block:: bash

      dataset:
         input_file: "dataset_pos_frc_ev.xyz"
         output_prefix: "consolidated_dataset" 
         sizes: [500, 1000, 2000, 4000]
      # Subset counts (number of structures from each method)
         subset_counts:
            MD: 2533   
            PCA: 1200 
            PCA_Surface: 600 
            Random: 200
            contamination: 0.05 
      SOAP: 
         species: ["In", "P", "Cl"] 
         r_cut: 12.0 
         n_max: 7 
         l_max: 3 
         sigma: 0.1 
        periodic: False
        sparse: False

Detailed Explanation of YAML Input :
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- ``input_file``: specifies the input file name.
- ``output_prefix``: specifies the prefix of the output files
- ``sizes``: creates chunks of different sizes.

Subset counts:

- ``MD``: structures obtained from Molecular Dynamics (MD) simulation. Adjust according to your data.
- ``PCA``: structures obtained from Principal Component Analysis (PCA).
- ``PCA_Surface``: surface-focused structures from PCA sampling.
- ``Random``: randomly selected structures for additional diversity.
- ``contamination``: fraction of outliers removed by Isolation Forest. 

The output files contain:
     * `consolidated_dataset`: a chunk of dataset with the most diverse structures (preferred for ML training).
     * `MD_random_dataset`: random structures picked from MD data.
     * `random_dataset`: random structures from the whole dataset.

Choose the subset preferred for your method and convert according `xyz` file to `npz` format using: 

  .. code-block:: bash

       python xyztonpz.py consolidated_dataset_1000.xyz
