This repository contains a specialized suite of high-performance C engines and Python analysis pipelines developed for studying the collective behavior of **Active Brownian Particles (ABP)** and **Hard Sphere Models**. 

The codebase is designed to bridge the gap between microscopic particle dynamics and macroscopic phase behavior, specifically focusing on non-equilibrium steady states.

---

## 🔬 Scientific Context & Research Focus

The tools in this repository investigate the following physical phenomena:

* **Motility-Induced Phase Separation (MIPS):** Analyzing how purely repulsive self-propelled particles undergo phase separation without attractive forces.
* **Structural Correlation Functions:** Quantifying spatial order using the Static Structure Factor $S(k)$ and Radial Distribution Function $g(r)$.
* **Dynamics & Anomalous Diffusion:** Characterizing transport properties via Mean Square Displacement (MSD) and the extraction of translational diffusion coefficients.
* **Topological Order:** Identifying local symmetries (e.g., Hexatic order) and global polarization in dense active systems.

---

## 📂 Project Organization

### 1. High-Performance Simulation Engines (`/core_simulations`)
Written in C for computational efficiency, implementing Langevin dynamics and Brownian integrators.
* `ABP.c`: Core engine for Active Brownian Particle dynamics.
* `hard_sphere_model.c`: Simulation of hard-core interactions and collision dynamics.
* `mt19937-64.c`: Implementation of the 64-bit Mersenne Twister for high-quality pseudo-random number generation (crucial for stochastic simulations).

### 2. Dynamics & Transport Analysis (`/analysis/msd`)
* `calculate_msd_data_save.py`: Optimized algorithm for computing time-averaged Mean Square Displacement.
* `trans_diff_coeff.py`: Computes the effective diffusion constant $D_{eff}$ from the long-time limit of MSD.
* `msd_from_0_ref.py`: Investigates cage effects and relaxation times.

### 3. Structural & Pair Correlations (`/analysis/structure`)
* `structure_factor_2D.c`: High-speed C implementation for calculating the structure factor $S(k)$ in 2D Fourier space.
* `RDF_latest.c`: Computes the Pair Correlation Function $g(r)$ to identify fluid, hexatic, or crystalline phases.
* `exact_structure.py`: Theoretical comparison and verification against analytical models.

### 4. Order Parameters & Visualization (`/analysis/order` & `/plots`)
* `polar_order_snapshots.py`: Computes global polarization $\Phi(t) = \frac{1}{N} \left| \sum_{i=1}^N \hat{n}_i(t) \right|$.
* `hexatic_order_com.py`: Calculates the bond-orientational order parameter $\psi_6$.
* `video_creation.py`: Tools for rendering simulation snapshots into high-quality MP4/GIF animation
