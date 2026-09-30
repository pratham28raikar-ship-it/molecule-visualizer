# 🧬 Molecule Visualizer

An interactive 3D web application for chemical and biological molecular visualization built with **Python**, **Streamlit**, and **Py3Dmol** (3Dmol.js).

![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)
![Streamlit](https://img.shields.io/badge/Streamlit-1.30%2B-red.svg)
![Py3Dmol](https://img.shields.io/badge/Py3Dmol-2.0%2B-green.svg)
![License](https://img.shields.io/badge/License-MIT-purple.svg)

---

## 🌟 Overview

**Molecule Visualizer** empowers researchers, students, and educators to upload molecular structure files and explore their 3D spatial conformations in real-time right in their web browser. The application parses atomic records, calculates physicochemical metrics, detects bonds, and renders high-fidelity WebGL graphics.

---

## ✨ Core Features

1. **Multi-Format Molecular File Support**:
   - **`.pdb` (Protein Data Bank)**: Macromolecules (proteins, enzymes, DNA/RNA) and small-molecule ligands.
   - **`.xyz` (Cartesian Coordinates)**: Standard computational chemistry format.

2. **Interactive 3D WebGL Canvas**:
   - 🖱️ **Rotate**: Left-click + drag to orbit in 3D.
   - 🔍 **Zoom**: Mouse wheel or trackpad pinch to zoom in/out.
   - ✋ **Pan**: Right-click + drag (or Shift + Left-click) to translate across the viewport.
   - 🔄 **Reset View**: Re-center and reset camera perspective with one click.
   - ✨ **Auto-Spin**: Continuous smooth rotation around the vertical axis.
   - 📸 **Snapshot Export**: Capture and download high-resolution PNG images directly from the 3D canvas.

3. **Versatile Visualization Styles**:
   - **Ball-and-stick**: Atomic centers with connecting cylindrical bonds.
   - **Stick**: Wireframe/cylinder bond skeleton emphasizing connectivity.
   - **Sphere (Spacefill)**: Van der Waals spheres demonstrating molecular volume and steric packing.
   - **Cartoon**: Secondary structure ribbons ($\alpha$-helices and $\beta$-sheets) for protein backbones.
   - **Wireframe**: Clean geometric line representation.
   - **Molecular Surface**: Optional overlay of Van der Waals (VDW), Solvent Accessible (SAS), or Solvent Excluded (SES) surfaces with adjustable opacity.

4. **Dynamic Color Schemes & Backgrounds**:
   - Standard CPK element palette.
   - Custom carbon accent palettes (Cyan, Green, Magenta, Orange).
   - Secondary structure rainbow spectrum.
   - Background options: Dark Slate, Charcoal, Midnight Blue, Pitch Black, Pure White, Light Slate, or any custom hex color.

5. **Physicochemical Metrics & Structure Info**:
   - Structure file name and format badge.
   - Total atom count and number of distinct elements.
   - Bond counts (from explicit `CONECT` records or estimated from covalent radii).
   - Molecular formula (Hill system notation).
   - Estimated molecular weight (g/mol).
   - Macromolecular metadata (chains and residue counts for PDB).

6. **In-Depth Analysis Tabs**:
   - **📊 Elemental Composition**: Breakdown table with atomic counts, atomic weights, mass percentages, and bar chart.
   - **🔬 Atom Coordinates Explorer**: Interactive table of atom coordinates $(x, y, z)$ with search filter and CSV download.
   - **📝 Raw Molecular File**: Syntax-highlighted source file viewer with download option.
   - **📖 Structure & Guide**: Educational guide for formats and navigation controls.

---

## 📁 Project Structure

```text
├── app.py                  # Main Streamlit application with Upload and Builder modes
├── builder.py              # Molecule Builder engine (state, chemistry validation, 3D optimizer, export)
├── parsers.py              # Molecular file parser (PDB & XYZ) with bond estimation
├── visualizer.py           # Py3Dmol WebGL rendering engine with interactive HUD
├── requirements.txt        # Python package dependencies (including RDKit)
├── README.md               # Documentation and usage guide
└── sample_files/           # Ready-to-use sample molecular files
    ├── sample.pdb          # 1CRN Crambin protein (secondary structure benchmark)
    ├── sample.xyz          # Caffeine (C8H10N4O2)
    ├── aspirin.pdb         # Acetylsalicylic acid with CONECT records
    └── benzene.xyz         # Benzene aromatic ring (C6H6)
```

---

## 🛠️ Interactive Molecule Builder (Avogadro Inspired)

Switch between **📁 Upload Molecule** and **🛠️ Build Molecule** modes using the top-level mode selector.

### Key Capabilities:
- **1. Elemental Construction**: Add common elements (**H, C, N, O, F, P, S, Cl, Br, I**) with standard CPK element colors and recognized sizes.
- **2. Smart Bond Creation**: Form and display **Single ($1$)**, **Double ($2$)**, and **Triple ($3$)** bonds rendered with distinct multiple cylinders in WebGL.
- **3. Editing & Selection**: Select atoms, change elements, tweak Cartesian coordinates, change bond orders, delete atoms and bonds, or clear the workspace.
- **4. Undo & Redo History**: Full multi-step snapshot undo/redo support.
- **5. Real-Time Chemistry Validation**: Automatic valence evaluation (H: 1, C: 4, N: 3, O: 2, Halogens: 1, P: 3/5, S: 2/4/6) with descriptive educational warnings if typical valences are exceeded.
- **6. 3D Structure Optimization**:
  - Leverages **RDKit ETKDGv3 conformer embedding** and **UFF (Universal Force Field)** energy minimization.
  - Option to auto-saturate open valences with explicit Hydrogens.
  - Built-in **pure-Python force-directed VSEPR relaxation** fallback guarantees robust 3D generation even for unusual or radical structures.
- **7. Export Formats**: Export built structures as valid **PDB** (with standard `CONECT` records) or **XYZ** Cartesian files.
- **8. Direct Visualizer Integration**: Click **🧬 Open in Visualizer** to transfer the built structure directly into the full multi-style visualizer!

---

## 🚀 Getting Started

### 1. Clone or Open the Repository
```bash
cd "CC PROJECT"
```

### 2. Set Up a Virtual Environment (Optional but Recommended)
```bash
# Windows
python -m venv .venv
.\.venv\Scripts\activate

# Linux / macOS
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

### 4. Run the Application
```bash
streamlit run app.py
```

The web application will launch automatically at `http://localhost:8501`.

---

## 🧪 Included Sample Files

You can test the visualizer instantly using the preloaded samples in `sample_files/`:
- **`sample.pdb`**: High-resolution Crambin plant protein (PDB ID: 1CRN, 327 atoms, $\alpha$-helices, $\beta$-sheets, and disulfide bridges).
- **`sample.xyz`**: Caffeine stimulant ($C_8H_{10}N_4O_2$, 24 atoms, 25 bonds).
- **`aspirin.pdb`**: Aspirin analgesic ($C_9H_8O_4$, 21 atoms with explicit `CONECT` bonds).
- **`benzene.xyz`**: Planar aromatic hydrocarbon ($C_6H_6$, 12 atoms).

---

## 🎮 Mouse & Navigation Controls

| Action | Mouse Gesture | Description |
| :--- | :--- | :--- |
| **Rotate** | `Left Click + Drag` | Rotates the molecule in 3D space |
| **Zoom** | `Scroll Wheel / Pinch` | Zooms in or out relative to camera |
| **Pan** | `Right Click + Drag` *(or `Shift + Left Click`)* | Translates/moves molecule horizontally and vertically |
| **Hover** | `Move Cursor over Atom` | Shows tooltip with element, serial #, and residue |
| **Reset** | Click **🔄 Reset** HUD or Sidebar | Centers molecule and restores initial view scale |
| **Spin** | Click **✨ Spin** HUD or Sidebar toggle | Toggles continuous auto-rotation |
| **Snapshot** | Click **📸 PNG** HUD | Downloads screenshot PNG of current view |

---

## 🛡️ Error Handling
- Uploading an invalid, corrupted, or non-molecular file presents a clear and helpful error notification instead of crashing the server.
- Selecting **Cartoon** style for an XYZ file gracefully falls back to Ball-and-Stick while displaying an informative tip explaining that XYZ format lacks secondary structure definitions.
