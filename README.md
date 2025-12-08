# matGoFEM

**matGoFEM** is a MATLAB front-end for the **GoFEM** package, aimed at users who prefer working in MATLAB instead of Python.

The main focus of matGoFEM is to simplify:

- Working with **GoFEM meshes and models** (e.g. VTU outputs),
- Building and inspecting **3-D resistivity models**,
- Producing **2-D/3-D visualisations** (cross-sections, maps, isosurfaces),
- Preparing and exploring models in the context of **magnetotelluric (MT) workflows**,
- Visualising and comparing input and output datasets.

If you are already comfortable in MATLAB but not in Python, matGoFEM provides a lightweight way to interact with GoFEM results and integrate them into your existing MATLAB toolchain.

---

## Motivation

GoFEM is a powerful finite-element framework for large-scale 3-D EM modelling and inversion.  
The GoFEM core and the Python front-end **pyGoFEM** have been developed by **Alexander Grayver**, see the original Python repository:

➡️ **pyGoFEM (Python front-end for GoFEM)**:  
https://github.com/GoFEM/pyGoFEM/

Existing tooling such as **pyGoFEM** provides a rich Python interface, but not all users are equally familiar with Python or Jupyter.

**matGoFEM** is **strongly inspired by pyGoFEM** and has two main goals:

1. Provide a **MATLAB-friendly interface** to GoFEM meshes and model files.
2. Offer **ready-to-run examples** (Live Scripts and Jupyter notebooks) showing typical post-processing and visualisation workflows for GoFEM models.

matGoFEM is not meant to replace pyGoFEM itself; instead, it complements them for users and projects that are heavily MATLAB-centric.

---

## Relation to FFMT and MT2GoFEM

Many magnetotelluric-specific functions that are not (yet) included in this repository are part of the **FFMT** ecosystem:

➡️ **FFMT project page**: https://www.ffmt.xyz/

In the long term, several matGoFEM capabilities are intended to be:

- Integrated into a GUI within FFMT (working title: **MT2GoFEM**),
- Exposed here in this GitHub repository to encourage **community participation and collaboration**.

In other words, matGoFEM focuses on GoFEM-oriented MATLAB workflows, while FFMT provides a broader MT processing and visualisation framework. Where possible, the two will be kept consistent and interoperable.

---

## Repository structure

The repository is organised as follows:

- `src/`  
  Core MATLAB functions for reading, processing and visualising GoFEM models (e.g. VTU readers, model interpolation, cross-sections, isosurfaces, MT-oriented utilities).

- `examples/`  
  Self-contained examples. Each subfolder typically contains:
  - A **MATLAB Live Script** (`*.mlx`) showing an interactive workflow.
  - A **Jupyter notebook** (`*.ipynb`) exported from the Live Script, which can be viewed directly on GitHub or run in a Jupyter environment.

- `data/`  
  Example input data or small test models (if provided).

- `docs/`  
  Additional documentation, exported HTML/PDF versions of Live Scripts, and notes.

---

## What matGoFEM does (and does not) do

**matGoFEM _does_:**

- Read and work with GoFEM meshes / model files (e.g. VTU outputs).
- Provide utilities to:
  - Interpolate models onto regular grids,
  - Extract and plot 2-D/3-D resistivity cross-sections,
  - Build isosurfaces and simple volume visualisations,
  - Prepare and inspect models used in MT workflows,
  - Compare input and output datasets in a consistent way.

**matGoFEM _does not_:**

- Distribute or compile GoFEM itself.
- Implement the full deal.II / GoFEM API in MATLAB.

You still need a **working GoFEM installation** (or access to GoFEM outputs) to produce the mesh and model files that matGoFEM will process.

---

## Requirements

- **MATLAB** R2023b or newer (developed and tested with R2024b).
- No proprietary toolboxes are strictly required for the basic examples, although some visualisation routines may benefit from:
  - Image Processing Toolbox
  - Mapping Toolbox  
  (optional, depending on your use case).

matGoFEM is platform-agnostic and should work on **Linux, Windows, and macOS**, as long as MATLAB is supported.

---

## Installation

There is no special installation step. The typical workflow is:

1. **Clone or download** this repository:
   - Using Git:
     ```bash
     git clone https://github.com/<your-user>/matGoFEM.git
     ```
   - Or download the ZIP from GitHub and extract it.

2. **Add matGoFEM to the MATLAB path**:
   - Open MATLAB and set the current folder to the root of the repository.
   - Either:
     - Use the MATLAB Project file (`.prj`), or
     - Run:
       ```matlab
       addpath(genpath(pwd));
       ```

3. Verify that functions in `src/` can be called from the MATLAB command window.

---

## Getting started

A good way to start is to run the examples:

1. In MATLAB, open the **Live Script** in `examples/LiveScript01`:
   - `LiveScript01.mlx`  
   This example demonstrates how to import a GoFEM VTU model and produce 2-D and 3-D resistivity cross-section plots.

2. Explore upcoming examples:
   - `LiveScript02` – Creating isosurfaces from GoFEM VTU models *(coming soon)*.

Each example folder also contains a corresponding **Jupyter notebook** (`*.ipynb`) exported from the Live Script.  
You can:

- Inspect the notebook directly on GitHub, or
- Download and run it in your preferred Jupyter environment.

---

## Contributing

Contributions are very welcome.

If you implement functions that could be useful to others (e.g. new readers, plotting routines, MT-specific helpers):

1. Fork this repository on GitHub.
2. Create a feature branch.
3. Add your functions with appropriate documentation and, if possible, a small example.
4. Open a **Pull Request** describing your changes.

Bug reports, feature requests, and discussion are also encouraged via GitHub Issues.

---

## Credits and origin

- **GoFEM** and the **pyGoFEM** Python front-end are developed by **Alexander Grayver**.  
  The original pyGoFEM repository is available at:  
  https://github.com/GoFEM/pyGoFEM/

- **matGoFEM** is an independent MATLAB front-end that is **conceptually inspired by pyGoFEM**, but implemented in MATLAB for users and workflows that are MATLAB-based.

Please acknowledge the original GoFEM / pyGoFEM work (and deal.II) in any publication where you use matGoFEM for visualisation or pre/post-processing.

---

## Referencing

Please acknowledge GoFEM and deal.II by citing one or more of the following publications:

- Grayver, A. V., 2015. Parallel 3D magnetotelluric inversion using adaptive finite-element method. Part I: theory and synthetic study. *Geophysical Journal International*, 202(1), 584–603, doi: 10.1093/gji/ggv165

- Grayver, A. V., & Kolev, T. V., 2015. Large-scale 3D geo-electromagnetic modeling using parallel adaptive high-order finite element method. *Geophysics*, 80(6), 277–291, doi: 10.1190/GEO2015-0013.1

- Arndt, D., Bangerth, W., Blais, B., Clevenger, T. C., Fehling, M., Grayver, A. V., *et al.*, 2020. The deal.II library, version 9.2. *Journal of Numerical Mathematics*, 28(3), 131–146.

If matGoFEM becomes useful for you in a publication, consider also referencing this repository and briefly describing how it was used (e.g. for visualisation or pre/post-processing).

---

## License

> _To be defined._  
> For now, please treat this code as **research/experimental software**. A proper open-source license (e.g. MIT/GPL/BSD) will be added once the structure of matGoFEM stabilises.