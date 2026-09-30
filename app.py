"""
Molecule Visualizer & Interactive Molecule Builder
Built with Python, Streamlit, Py3Dmol, and RDKit.
Supports 3D molecular inspection (PDB/XYZ) and interactive molecule construction.
"""

from __future__ import annotations
import os
import json
import streamlit as st
import streamlit.components.v1 as components

import parsers
import visualizer
import builder

# Set page configuration
st.set_page_config(
    page_title="Molecule Visualizer & Builder",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for modern scientific UI
st.markdown("""
<style>
    /* Global styling enhancements */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
    
    html, body, [class*="css"] {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
    }
    
    /* Top Mode Navigation Bar */
    .mode-nav-container {
        display: flex;
        align-items: center;
        gap: 0.75rem;
        background: rgba(15, 23, 42, 0.7);
        border: 1px solid rgba(56, 189, 248, 0.25);
        border-radius: 12px;
        padding: 0.4rem 0.5rem;
        margin-bottom: 1rem;
        backdrop-filter: blur(8px);
    }
    
    /* Hero section styling */
    .hero-container {
        padding: 1.1rem 1.4rem;
        background: linear-gradient(135deg, rgba(15, 23, 42, 0.95) 0%, rgba(30, 41, 59, 0.9) 100%);
        border: 1px solid rgba(56, 189, 248, 0.2);
        border-radius: 14px;
        margin-bottom: 1.1rem;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.15);
    }
    .hero-title {
        font-size: 1.85rem;
        font-weight: 700;
        letter-spacing: -0.02em;
        background: linear-gradient(90deg, #38bdf8, #818cf8, #c084fc);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin: 0 0 0.35rem 0;
        display: flex;
        align-items: center;
        gap: 0.6rem;
    }
    .hero-desc {
        color: #94a3b8;
        font-size: 0.92rem;
        margin: 0;
        line-height: 1.45;
    }
    
    /* Metric cards styling */
    .metric-card {
        background: rgba(30, 41, 59, 0.6);
        border: 1px solid rgba(148, 163, 184, 0.15);
        border-radius: 10px;
        padding: 0.85rem 1rem;
        transition: transform 0.15s ease, border-color 0.15s ease;
    }
    .metric-card:hover {
        border-color: rgba(56, 189, 248, 0.4);
        transform: translateY(-2px);
    }
    .metric-label {
        font-size: 0.75rem;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        color: #94a3b8;
        margin-bottom: 0.2rem;
        font-weight: 600;
    }
    .metric-value {
        font-size: 1.35rem;
        font-weight: 700;
        color: #f8fafc;
        font-family: 'JetBrains Mono', monospace;
    }
    .metric-sub {
        font-size: 0.74rem;
        color: #64748b;
        margin-top: 0.2rem;
    }
    
    /* Tag pills */
    .badge-pill {
        display: inline-block;
        padding: 0.18rem 0.6rem;
        border-radius: 9999px;
        font-size: 0.73rem;
        font-weight: 600;
        margin-right: 0.35rem;
    }
    .badge-pdb {
        background: rgba(14, 165, 233, 0.15);
        color: #38bdf8;
        border: 1px solid rgba(14, 165, 233, 0.3);
    }
    .badge-xyz {
        background: rgba(168, 85, 247, 0.15);
        color: #c084fc;
        border: 1px solid rgba(168, 85, 247, 0.3);
    }
    .badge-builder {
        background: rgba(16, 185, 129, 0.15);
        color: #34d399;
        border: 1px solid rgba(16, 185, 129, 0.3);
    }
    
    /* Canvas instruction notice */
    .canvas-notice {
        background: rgba(15, 23, 42, 0.5);
        border: 1px solid rgba(148, 163, 184, 0.12);
        border-radius: 8px;
        padding: 0.45rem 0.9rem;
        margin-bottom: 0.5rem;
        font-size: 0.8rem;
        color: #94a3b8;
        display: flex;
        align-items: center;
        justify-content: space-between;
    }

    /* Avogadro tool container styling */
    .tool-box {
        background: rgba(30, 41, 59, 0.45);
        border: 1px solid rgba(148, 163, 184, 0.15);
        border-radius: 10px;
        padding: 0.9rem;
        margin-bottom: 0.8rem;
    }
    
    .elem-chip {
        display: inline-flex;
        align-items: center;
        gap: 5px;
        padding: 4px 10px;
        border-radius: 6px;
        background: rgba(15, 23, 42, 0.6);
        border: 1px solid rgba(148, 163, 184, 0.2);
        font-size: 0.85rem;
        font-weight: 600;
        margin-right: 6px;
        margin-bottom: 6px;
    }

    .elem-dot {
        width: 10px;
        height: 10px;
        border-radius: 50%;
        display: inline-block;
    }

    /* Code block styling */
    pre, code {
        font-family: 'JetBrains Mono', monospace !important;
    }
</style>
""", unsafe_allow_html=True)

# ----------------- SESSION STATE INITIALIZATION -----------------
if "app_mode" not in st.session_state:
    st.session_state.app_mode = "📁 Upload Molecule"
if "reset_count" not in st.session_state:
    st.session_state.reset_count = 0
if "builder_reset_count" not in st.session_state:
    st.session_state.builder_reset_count = 0
if "builder_state" not in st.session_state:
    st.session_state.builder_state = builder.create_initial_state()
if "builder_opt_message" not in st.session_state:
    st.session_state.builder_opt_message = None
if "built_mol_content" not in st.session_state:
    st.session_state.built_mol_content = None
if "built_mol_filename" not in st.session_state:
    st.session_state.built_mol_filename = None
if "select_built_source" not in st.session_state:
    st.session_state.select_built_source = False

# Sample Files definition
SAMPLE_DIR = "sample_files"
SAMPLES = {
    "Crambin Protein (PDB)": os.path.join(SAMPLE_DIR, "sample.pdb"),
    "Caffeine Molecule (XYZ)": os.path.join(SAMPLE_DIR, "sample.xyz"),
    "Aspirin Molecule (PDB)": os.path.join(SAMPLE_DIR, "aspirin.pdb"),
    "Benzene Ring (XYZ)": os.path.join(SAMPLE_DIR, "benzene.xyz")
}

# ----------------- TOP MODE SELECTOR -----------------
# Prominent top tabs: [ Upload Molecule ] [ Build Molecule ]
st.write("")
col_nav_left, col_nav_right = st.columns([1, 1], gap="small")

with col_nav_left:
    is_upload_active = (st.session_state.app_mode == "📁 Upload Molecule")
    if st.button(
        "📁 Upload Molecule",
        key="nav_upload_btn",
        use_container_width=True,
        type="primary" if is_upload_active else "secondary",
        help="View preloaded samples or upload custom PDB / XYZ files"
    ):
        st.session_state.app_mode = "📁 Upload Molecule"
        st.rerun()

with col_nav_right:
    is_builder_active = (st.session_state.app_mode == "🛠️ Build Molecule")
    if st.button(
        "🛠️ Build Molecule",
        key="nav_builder_btn",
        use_container_width=True,
        type="primary" if is_builder_active else "secondary",
        help="Interactive Avogadro-style 3D molecule builder workspace"
    ):
        st.session_state.app_mode = "🛠️ Build Molecule"
        st.rerun()

st.write("")


# =========================================================================
# MODE 1: UPLOAD & EXPLORE MOLECULE (ORIGINAL FUNCTIONALITY UNCHANGED)
# =========================================================================
if st.session_state.app_mode == "📁 Upload Molecule":

    # ----------------- SIDEBAR: VISUALIZER -----------------
    with st.sidebar:
        st.markdown("### 🧬 Molecule Visualizer")
        st.caption("3D Structural Chemistry & Biology Explorer")
        st.divider()

        # Input Mode Selection
        source_options = ["🧪 Sample Molecules", "📁 Upload Custom File"]
        if st.session_state.built_mol_content is not None:
            source_options.append("🛠️ Built Molecule")

        default_idx = 0
        if st.session_state.select_built_source and "🛠️ Built Molecule" in source_options:
            default_idx = source_options.index("🛠️ Built Molecule")
            st.session_state.select_built_source = False

        input_source = st.radio(
            "Molecule Source",
            source_options,
            index=default_idx,
            help="Choose a preloaded scientific molecule, upload your own .pdb or .xyz file, or load the built molecule"
        )

        filename = ""
        file_content = ""
        file_bytes = None

        if input_source == "🧪 Sample Molecules":
            selected_sample = st.selectbox(
                "Select Sample",
                list(SAMPLES.keys()),
                index=0,
                help="Preloaded molecular structures covering both proteins (PDB) and small molecules (XYZ/PDB)"
            )
            sample_path = SAMPLES[selected_sample]
            if os.path.exists(sample_path):
                with open(sample_path, "r", encoding="utf-8") as f:
                    file_content = f.read()
                filename = os.path.basename(sample_path)
            else:
                st.error(f"Sample file '{sample_path}' not found.")

        elif input_source == "🛠️ Built Molecule":
            file_content = st.session_state.built_mol_content or ""
            filename = st.session_state.built_mol_filename or "builder_molecule.pdb"
            st.success(f"Loaded structure from Molecule Builder: `{filename}`")
            if st.button("✏️ Return to Builder Workspace", use_container_width=True):
                st.session_state.app_mode = "🛠️ Build Molecule"
                st.rerun()

        else:
            uploaded_file = st.file_uploader(
                "Upload Molecular File",
                type=["pdb", "xyz"],
                help="Upload a valid .pdb (Protein Data Bank) or .xyz (Cartesian Coordinates) file"
            )
            if uploaded_file is not None:
                filename = uploaded_file.name
                try:
                    file_bytes = uploaded_file.getvalue()
                    file_content = file_bytes.decode("utf-8")
                except UnicodeDecodeError:
                    try:
                        file_content = file_bytes.decode("latin-1")
                    except Exception as e:
                        st.error(f"Error reading file encoding: {e}")
            else:
                st.info("👆 Please upload a .pdb or .xyz file above to begin 3D visualization.")

        st.divider()
        st.markdown("#### 🎨 Visualization Style")

        viz_style = st.selectbox(
            "Rendering Style",
            [
                "Ball-and-stick",
                "Stick",
                "Sphere (Spacefill)",
                "Cartoon (Proteins/Polymers)",
                "Wireframe"
            ],
            index=0,
            help="Choose molecular representation: Ball-and-stick, Stick, Sphere/Spacefill, Cartoon, or Wireframe."
        )

        color_scheme = st.selectbox(
            "Color Scheme",
            list(visualizer.COLOR_SCHEMES.keys()),
            index=0,
            help="Color atoms by standard CPK element table, carbon accent, or secondary structure spectrum"
        )

        st.markdown("#### 🖼️ Background Color")
        bg_selection = st.selectbox(
            "Preset Backgrounds",
            list(visualizer.BG_PRESETS.keys()) + ["Custom Color..."],
            index=0
        )

        if bg_selection == "Custom Color...":
            bg_color = st.color_picker("Pick Custom Background", "#0f172a")
        else:
            bg_color = visualizer.BG_PRESETS[bg_selection]

        with st.expander("⚙️ Advanced Display Options", expanded=False):
            show_spin = st.toggle("Auto-Spin Animation", value=False, help="Continuous rotation around the Y-axis")
            show_hover = st.toggle("Atom Hover Tooltips", value=True, help="Display element and serial when hovering")
            show_labels = st.toggle("Atom 3D Text Labels", value=False, help="Render persistent atom labels on the structure")
            
            st.markdown("**Molecular Surface Overlay**")
            surface_type = st.selectbox(
                "Surface Type",
                ["None", "VDW (Van der Waals)", "SAS (Solvent Accessible)", "SES (Solvent Excluded)"],
                index=0
            )
            surface_opacity = 0.5
            if surface_type != "None":
                surface_opacity = st.slider("Surface Opacity", 0.1, 1.0, 0.45, 0.05)

            viewer_height = st.slider("Viewer Height (px)", min_value=400, max_value=850, value=580, step=20)

        st.divider()
        # Reset View Button
        if st.sidebar.button("🔄 Reset 3D View", use_container_width=True, help="Reset camera angle, zoom, and orientation"):
            st.session_state.reset_count += 1
            st.rerun()

    # ----------------- MAIN VIEW: VISUALIZER -----------------
    st.markdown("""
    <div class="hero-container">
        <div class="hero-title">
            <span>🧬</span> Molecule Visualizer
        </div>
        <p class="hero-desc">
            Explore 3D chemical and biomolecular architectures interactively. Supports atomic structures in 
            <b>.pdb</b> (Protein Data Bank) and <b>.xyz</b> (Cartesian coordinate) formats with WebGL rendering.
        </p>
    </div>
    """, unsafe_allow_html=True)

    mol_data = None
    parse_error = None

    if file_content and filename:
        try:
            mol_data = parsers.parse_molecular_file(filename, file_content)
        except Exception as e:
            parse_error = str(e)

    if parse_error:
        st.error(f"❌ Failed to parse '{filename}': {parse_error}")
        st.info("Please verify that the uploaded file follows standard PDB or XYZ specifications.")

    elif mol_data:
        # Structure Metrics Bar
        badge_class = "badge-pdb" if mol_data["file_type"] == "PDB" else "badge-xyz"
        if input_source == "🛠️ Built Molecule":
            badge_class = "badge-builder"

        col1, col2, col3, col4, col5 = st.columns(5)

        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">Structure File</div>
                <div class="metric-value" style="font-size: 1.1rem; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;" title="{mol_data['filename']}">
                    {mol_data['filename']}
                </div>
                <div class="metric-sub"><span class="badge-pill {badge_class}">{mol_data['file_type']} FORMAT</span></div>
            </div>
            """, unsafe_allow_html=True)

        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">Total Atoms</div>
                <div class="metric-value">{mol_data['num_atoms']:,}</div>
                <div class="metric-sub">{len(mol_data['element_counts'])} distinct elements</div>
            </div>
            """, unsafe_allow_html=True)

        with col3:
            bonds_val = mol_data['num_bonds']
            bonds_disp = f"{bonds_val:,}" if isinstance(bonds_val, int) and bonds_val >= 0 else "Dynamic"
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">Bonds</div>
                <div class="metric-value">{bonds_disp}</div>
                <div class="metric-sub">{mol_data['bonds_source']}</div>
            </div>
            """, unsafe_allow_html=True)

        with col4:
            mw = mol_data['molecular_weight']
            mw_str = f"{mw:,.2f}" if mw > 0 else "N/A"
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">Mol. Weight</div>
                <div class="metric-value">{mw_str}</div>
                <div class="metric-sub">g / mol (estimated)</div>
            </div>
            """, unsafe_allow_html=True)

        with col5:
            if mol_data["file_type"] == "PDB" and mol_data.get("num_residues", 0) > 1:
                res_count = mol_data.get('num_residues', 'N/A')
                chain_list = ", ".join(mol_data.get('chains', ['A']))
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-label">Residues & Chains</div>
                    <div class="metric-value">{res_count}</div>
                    <div class="metric-sub">Chain(s): {chain_list}</div>
                </div>
                """, unsafe_allow_html=True)
            else:
                formula_disp = mol_data['formula'] if mol_data['formula'] else "N/A"
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-label">Chemical Formula</div>
                    <div class="metric-value" style="font-size: 1.15rem;">{formula_disp}</div>
                    <div class="metric-sub">Hill system notation</div>
                </div>
                """, unsafe_allow_html=True)

        st.write("")

        # 3D Viewer Container
        st.markdown("""
        <div class="canvas-notice">
            <div><b>Interactive Controls:</b> 🖱️ <b>Left-Click + Drag:</b> Rotate &nbsp;|&nbsp; 🔍 <b>Scroll / Pinch:</b> Zoom &nbsp;|&nbsp; ✋ <b>Right-Click / Shift + Drag:</b> Pan</div>
            <div><i>Tip: Use HUD buttons inside the viewer to Reset, Spin, or capture PNG</i></div>
        </div>
        """, unsafe_allow_html=True)

        if "cartoon" in viz_style.lower() and mol_data["file_type"] == "XYZ":
            st.warning("ℹ️ **Cartoon representation** requires secondary structure records from macromolecular PDB files. Small molecule XYZ files are displayed in Ball-and-Stick mode instead.")

        try:
            viewer_html = visualizer.create_molecule_viewer(
                mol_data=mol_data,
                style=viz_style,
                bg_color=bg_color,
                color_scheme=color_scheme,
                show_spin=show_spin,
                show_hover=show_hover,
                show_labels=show_labels,
                surface_type=surface_type,
                surface_opacity=surface_opacity,
                viewer_height=viewer_height
            )
            render_key = f"mol_viewer_{mol_data['filename']}_{st.session_state.reset_count}_{viz_style}_{bg_color}"
            components.html(viewer_html, height=viewer_height + 8, scrolling=False)
        except Exception as e:
            st.error(f"Error initializing 3D viewer: {e}")

        st.write("")

        # Structure Details Tabs
        tab_comp, tab_coords, tab_raw, tab_about = st.tabs([
            "📊 Elemental Composition",
            "🔬 Atom Coordinates Explorer",
            "📝 Raw Molecular File",
            "📖 Structure & Guide"
        ])

        with tab_comp:
            st.markdown("#### Elemental Breakdown")
            elem_counts = mol_data["element_counts"]
            
            if elem_counts:
                total_atoms = mol_data["num_atoms"]
                total_mw = mol_data["molecular_weight"] or 1.0
                
                comp_rows = []
                for elem, count in sorted(elem_counts.items(), key=lambda x: -x[1]):
                    weight_single = parsers.ATOMIC_WEIGHTS.get(elem.upper(), 12.0)
                    mass_subtotal = weight_single * count
                    atom_pct = (count / total_atoms) * 100.0
                    mass_pct = (mass_subtotal / total_mw) * 100.0 if total_mw > 0 else 0.0
                    cpk_color = parsers.CPK_COLORS.get(elem.upper(), "#909090")
                    
                    comp_rows.append({
                        "Element": elem,
                        "Count": count,
                        "Atom %": f"{atom_pct:.1f}%",
                        "Atomic Mass (u)": f"{weight_single:.3f}",
                        "Total Mass (g/mol)": f"{mass_subtotal:.2f}",
                        "Mass %": f"{mass_pct:.1f}%",
                        "Color": cpk_color
                    })

                c_left, c_right = st.columns([1, 1])
                with c_left:
                    st.markdown("**Element Distribution Table**")
                    table_display = [{k: v for k, v in row.items() if k != "Color"} for row in comp_rows]
                    st.table(table_display)

                with c_right:
                    st.markdown("**Atom Count Chart**")
                    chart_data = {row["Element"]: row["Count"] for row in comp_rows}
                    st.bar_chart(chart_data)

        with tab_coords:
            st.markdown("#### Atom Coordinates Table")
            atoms_list = mol_data["atoms"]
            
            search_query = st.text_input("Filter by element or atom name (e.g., 'C', 'N', 'CA')", "")
            filtered_atoms = atoms_list
            if search_query:
                q = search_query.strip().upper()
                filtered_atoms = [a for a in atoms_list if q in a["element"].upper() or q in a["name"].upper()]

            st.caption(f"Showing {len(filtered_atoms)} of {len(atoms_list)} atoms")
            
            display_cols = ["serial", "element", "name", "x", "y", "z"]
            if mol_data["file_type"] == "PDB":
                display_cols.extend(["resName", "chain", "resSeq"])
                
            display_data = [{k: a.get(k, "") for k in display_cols} for a in filtered_atoms[:500]]
            st.dataframe(display_data, use_container_width=True, height=360)
            
            if len(filtered_atoms) > 500:
                st.info(f"Showing first 500 atoms out of {len(filtered_atoms)}.")

            csv_header = ",".join(display_cols) + "\n"
            csv_lines = [",".join(str(a.get(k, "")) for k in display_cols) for a in atoms_list]
            csv_content = csv_header + "\n".join(csv_lines)
            
            st.download_button(
                label="📥 Download Atom Coordinates (CSV)",
                data=csv_content,
                file_name=f"{mol_data['filename']}_coordinates.csv",
                mime="text/csv"
            )

        with tab_raw:
            st.markdown(f"#### Source File: `{mol_data['filename']}`")
            st.download_button(
                label=f"📥 Download {mol_data['filename']}",
                data=mol_data["raw_content"],
                file_name=mol_data["filename"],
                mime="text/plain"
            )
            raw_text = mol_data["raw_content"]
            lines = raw_text.splitlines()
            preview_text = "\n".join(lines[:250])
            if len(lines) > 250:
                preview_text += f"\n\n... [{len(lines) - 250} more lines omitted for preview]"
            st.code(preview_text, language="text")

        with tab_about:
            st.markdown("""
            ### About Molecule Visualizer
            
            **Molecule Visualizer** is a high-performance web tool designed for researchers, educators, and students 
            to inspect 3D spatial conformations of chemicals, drugs, and biomacromolecules.
            
            #### Supported Formats:
            - **`.pdb` (Protein Data Bank)**: Standard coordinate file for biological molecules including proteins, nucleic acids, and small molecule complexes.
            - **`.xyz` (Cartesian Coordinates)**: Classic computational chemistry format containing atom counts and 3D $(x, y, z)$ coordinates.
            
            #### Rendering Capabilities:
            - **Ball-and-Stick**: Explicit representation of atomic centers and chemical bonds with multiple bond cylinders.
            - **Stick**: Cylindrical bond skeleton emphasizing connectivity.
            - **Sphere (Spacefill)**: Van der Waals spheres illustrating molecular volume and steric hindrance.
            - **Cartoon**: Ribbon representation highlighting secondary structure ($\alpha$-helices and $\beta$-sheets).
            - **Wireframe**: Minimalist line drawing for rapid geometric inspection.
            """)

    else:
        st.info("👈 Select a sample molecule or upload a `.pdb` or `.xyz` file in the sidebar to begin.")


# =========================================================================
# MODE 2: MOLECULE BUILDER (AVOGADRO-INSPIRED WORKSPACE)
# =========================================================================
else:
    b_state = st.session_state.builder_state

    # ----------------- SIDEBAR: BUILDER SETTINGS & TEMPLATES -----------------
    with st.sidebar:
        st.markdown("### 🛠️ Molecule Builder")
        st.caption("Avogadro-Inspired 3D Construction Workspace")
        st.divider()

        st.markdown("#### ⚡ Quick Starting Templates")
        selected_tmpl = st.selectbox(
            "Load Molecular Template",
            ["-- Select a Template --"] + list(builder.TEMPLATES.keys()),
            index=0,
            help="Pre-loaded starter molecules (Water, Methane, Ethylene, Acetylene, Ethanol, Benzene)"
        )
        if selected_tmpl != "-- Select a Template --":
            if st.button(f"📥 Load {selected_tmpl.split()[0]}", use_container_width=True):
                builder.load_template_molecule(b_state, selected_tmpl)
                st.session_state.builder_opt_message = f"Loaded template: {selected_tmpl}"
                st.rerun()

        st.divider()
        st.markdown("#### 🎨 3D Canvas Settings")

        builder_bg = st.selectbox(
            "Background Color",
            list(visualizer.BG_PRESETS.keys()) + ["Custom Color..."],
            index=0,
            key="builder_bg_sel"
        )
        if builder_bg == "Custom Color...":
            b_bg_color = st.color_picker("Pick Custom Background", "#0f172a", key="builder_bg_cp")
        else:
            b_bg_color = visualizer.BG_PRESETS[builder_bg]

        b_viewer_height = st.slider("Canvas Height (px)", min_value=420, max_value=850, value=580, step=20, key="builder_h_slider")

        st.divider()
        if st.button("🔄 Reset Builder View", use_container_width=True, help="Reset camera zoom and rotation in builder"):
            st.session_state.builder_reset_count += 1
            st.rerun()

        st.divider()
        with st.expander("📖 Builder Workflow Help", expanded=False):
            st.markdown("""
            **Avogadro Workflow:**
            1. **Select an Element**: Choose H, C, N, O, F, P, S, Cl, Br, or I.
            2. **Add Atoms**: Attach directly to an existing atom with Single, Double, or Triple bonds, or place isolated atoms.
            3. **Form Bonds**: Connect any two atoms and set the bond order.
            4. **Valence Check**: Watch real-time chemistry validation. If an atom exceeds its typical valence, a warning is raised.
            5. **Optimize 3D**: Click **Optimize 3D Structure** to run RDKit UFF / VSEPR geometry relaxation!
            6. **Export / Visualize**: Export as PDB or XYZ, or open directly in the Visualizer!
            """)

    # ----------------- MAIN VIEW: BUILDER -----------------
    st.markdown("""
    <div class="hero-container">
        <div class="hero-title">
            <span>🛠️</span> Molecule Builder
        </div>
        <p class="hero-desc">
            Construct and manipulate custom 3D molecules interactively. Assemble atoms, form single, double, 
            and triple bonds, inspect chemistry valence rules, and optimize realistic 3D geometries with RDKit.
        </p>
    </div>
    """, unsafe_allow_html=True)

    # ----------------- TOP ACTION TOOLBAR -----------------
    num_history = len(b_state["history"])
    num_redo = len(b_state["redo_history"])
    mol_info = builder.compute_molecular_info(b_state)

    tool_c1, tool_c2, tool_c3, tool_c4, tool_c5, tool_c6 = st.columns([1, 1, 1, 1.3, 1.2, 1.5], gap="small")

    with tool_c1:
        if st.button(f"↩️ Undo ({num_history})", disabled=(num_history == 0), use_container_width=True, help="Revert the last modification"):
            if builder.undo(b_state):
                st.session_state.builder_opt_message = "Action undone."
                st.rerun()

    with tool_c2:
        if st.button(f"↪️ Redo ({num_redo})", disabled=(num_redo == 0), use_container_width=True, help="Redo the previously undone action"):
            if builder.redo(b_state):
                st.session_state.builder_opt_message = "Action redone."
                st.rerun()

    with tool_c3:
        if st.button("🗑️ Clear", disabled=(len(b_state["atoms"]) == 0), use_container_width=True, help="Clear all atoms and bonds in the molecule"):
            builder.clear_molecule(b_state)
            st.session_state.builder_opt_message = "Molecule cleared."
            st.rerun()

    with tool_c4:
        pdb_data_export = builder.export_to_pdb(b_state)
        st.download_button(
            "📥 Export PDB",
            data=pdb_data_export,
            file_name=f"{mol_info['formula'] or 'molecule'}.pdb",
            mime="chemical/x-pdb",
            disabled=(len(b_state["atoms"]) == 0),
            use_container_width=True,
            help="Download as standard Protein Data Bank file with CONECT records"
        )

    with tool_c5:
        xyz_data_export = builder.export_to_xyz(b_state)
        st.download_button(
            "📥 Export XYZ",
            data=xyz_data_export,
            file_name=f"{mol_info['formula'] or 'molecule'}.xyz",
            mime="chemical/x-xyz",
            disabled=(len(b_state["atoms"]) == 0),
            use_container_width=True,
            help="Download as Cartesian coordinate format (.xyz)"
        )

    with tool_c6:
        if st.button("🧬 Open in Visualizer", type="primary", disabled=(len(b_state["atoms"]) == 0), use_container_width=True, help="Load this constructed molecule into the main 3D Visualizer"):
            pdb_str = builder.export_to_pdb(b_state, title="Constructed Molecule")
            form_str = mol_info['formula'] if mol_info['formula'] != "Empty" else "molecule"
            st.session_state.built_mol_content = pdb_str
            st.session_state.built_mol_filename = f"{form_str}_built.pdb"
            st.session_state.select_built_source = True
            st.session_state.app_mode = "📁 Upload Molecule"
            st.rerun()

    # Feedback / Action Notice
    if st.session_state.builder_opt_message:
        st.info(st.session_state.builder_opt_message)
        st.session_state.builder_opt_message = None

    # Valence Warnings & Suggestions
    val_feedback = builder.validate_valences(b_state)
    val_warnings = [v for v in val_feedback if v["status"] == "warning"]
    val_infos = [v for v in val_feedback if v["status"] == "info"]

    if val_warnings:
        for w in val_warnings:
            st.warning(w["message"])

    # ----------------- 2-COLUMN WORKSPACE -----------------
    canvas_col, tools_col = st.columns([1.2, 1.0], gap="medium")

    # ----- LEFT COLUMN: 3D INTERACTIVE CANVAS & METRICS -----
    with canvas_col:
        # Navigation guidance banner
        st.markdown("""
        <div class="canvas-notice">
            <div>🖱️ <b>Rotate:</b> Left-Click &nbsp;|&nbsp; 🔍 <b>Zoom:</b> Scroll &nbsp;|&nbsp; ✋ <b>Pan:</b> Right-Click &nbsp;|&nbsp; ⭐ <b>Yellow Halo:</b> Selected Atom</div>
            <div><i>Multiple sticks denote Double & Triple bonds</i></div>
        </div>
        """, unsafe_allow_html=True)

        # 3D Py3Dmol Viewer
        if len(b_state["atoms"]) == 0:
            st.markdown(f"""
            <div style="height: {b_viewer_height}px; background: {b_bg_color}; border: 2px dashed rgba(56, 189, 248, 0.25); border-radius: 12px; display: flex; flex-direction: column; align-items: center; justify-content: center; color: #94a3b8; text-align: center; padding: 2rem;">
                <div style="font-size: 3rem; margin-bottom: 0.8rem;">🧪</div>
                <div style="font-size: 1.25rem; font-weight: 600; color: #f8fafc; margin-bottom: 0.5rem;">3D Workspace is Empty</div>
                <div style="font-size: 0.9rem; max-width: 420px; line-height: 1.5;">
                    Use the <b>➕ Add Atom</b> panel on the right to place your first atom, or select a starter template from the sidebar.
                </div>
            </div>
            """, unsafe_allow_html=True)
        else:
            try:
                builder_html = builder.create_builder_viewer(
                    state=b_state,
                    selected_atom_id=b_state.get("selected_atom_id"),
                    style="Ball-and-stick",
                    bg_color=b_bg_color,
                    viewer_height=b_viewer_height
                )
                components.html(builder_html, height=b_viewer_height + 8, scrolling=False)
            except Exception as e:
                st.error(f"Error rendering builder 3D view: {e}")

        st.write("")

        # Molecular Metrics Bar
        m_col1, m_col2, m_col3, m_col4 = st.columns(4)

        with m_col1:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">Formula</div>
                <div class="metric-value">{mol_info['formula']}</div>
                <div class="metric-sub">Hill system notation</div>
            </div>
            """, unsafe_allow_html=True)

        with m_col2:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">Total Atoms</div>
                <div class="metric-value">{mol_info['num_atoms']}</div>
                <div class="metric-sub">{len(mol_info['elements_present'])} distinct element(s)</div>
            </div>
            """, unsafe_allow_html=True)

        with m_col3:
            b_breakdown = f"{mol_info['single_bonds']}S / {mol_info['double_bonds']}D / {mol_info['triple_bonds']}T"
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">Total Bonds</div>
                <div class="metric-value">{mol_info['num_bonds']}</div>
                <div class="metric-sub">{b_breakdown}</div>
            </div>
            """, unsafe_allow_html=True)

        with m_col4:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">Mol. Weight</div>
                <div class="metric-value">{mol_info['molecular_weight']:.2f}</div>
                <div class="metric-sub">g / mol (estimated)</div>
            </div>
            """, unsafe_allow_html=True)

    # ----- RIGHT COLUMN: AVOGADRO-STYLE INTERACTIVE TOOLS -----
    with tools_col:
        st.markdown("#### 🛠️ Construction Tools")

        tab_add, tab_bond, tab_edit, tab_opt = st.tabs([
            "➕ Add Atom",
            "🔗 Create Bond",
            "✏️ Edit & Delete",
            "✨ Optimize 3D"
        ])

        # ---------- TAB 1: ADD ATOM TOOL ----------
        with tab_add:
            st.markdown("**1. Select Element:**")
            
            # Element picker chips with colors and typical valences
            selected_element = st.selectbox(
                "Element Symbol",
                builder.COMMON_ELEMENTS,
                index=1, # Default C
                format_func=lambda elem: f"{elem} — {parsers.clean_element_symbol(elem)} (Typical Valence: {', '.join(str(v) for v in builder.ELEMENT_VALENCES.get(elem.upper(), [1]))})",
                help="Choose from common chemical elements: H, C, N, O, F, P, S, Cl, Br, I"
            )

            # Show CPK color badge preview
            cpk_hex = parsers.CPK_COLORS.get(selected_element.upper(), "#909090")
            st.markdown(f"""
            <div style="display: flex; align-items: center; gap: 8px; margin-bottom: 0.8rem;">
                <span class="elem-dot" style="background-color: {cpk_hex}; width: 14px; height: 14px;"></span>
                <span style="font-size: 0.88rem; color: #f8fafc; font-weight: 600;">{selected_element} Selected</span>
                <span style="font-size: 0.78rem; color: #94a3b8;">(Atomic Mass: {parsers.ATOMIC_WEIGHTS.get(selected_element.upper(), 12.0):.3f} u)</span>
            </div>
            """, unsafe_allow_html=True)

            # Placement Mode
            if len(b_state["atoms"]) == 0:
                st.info("Placing the first atom into the origin `(0.0, 0.0, 0.0)`.")
                if st.button("➕ Place First Atom", type="primary", use_container_width=True):
                    builder.add_atom(b_state, element=selected_element, x=0.0, y=0.0, z=0.0)
                    st.rerun()
            else:
                placement_mode = st.radio(
                    "Placement Method",
                    ["🔗 Attach to an Atom (Avogadro Style)", "📍 Place as Isolated Atom"],
                    index=0,
                    help="Attach bonded to an existing atom or position as an independent coordinate"
                )

                if "Attach" in placement_mode:
                    atom_options = {a["id"]: f"Atom #{a['id']} ({a['element']}) — Bonds: {builder.get_atom_bond_order_sum(b_state, a['id'])}" for a in b_state["atoms"]}
                    
                    # Pre-select currently selected atom if valid
                    curr_sel = b_state.get("selected_atom_id")
                    default_parent_idx = 0
                    if curr_sel in atom_options:
                        default_parent_idx = list(atom_options.keys()).index(curr_sel)

                    parent_id = st.selectbox(
                        "Attach To Parent Atom:",
                        options=list(atom_options.keys()),
                        index=default_parent_idx,
                        format_func=lambda x: atom_options[x]
                    )

                    bond_order_choice = st.radio(
                        "Connecting Bond Order:",
                        ["Single (1)", "Double (2)", "Triple (3)"],
                        index=0,
                        horizontal=True
                    )
                    order_val = 1 if "Single" in bond_order_choice else (2 if "Double" in bond_order_choice else 3)

                    if st.button("➕ Add & Attach Atom", type="primary", use_container_width=True):
                        new_id = builder.add_atom(
                            b_state,
                            element=selected_element,
                            attach_to_id=parent_id,
                            bond_order=order_val
                        )
                        st.session_state.builder_opt_message = f"Added {selected_element} #{new_id} bonded to Atom #{parent_id} (Order: {order_val})."
                        st.rerun()

                else:
                    st.markdown("**Isolated Atom Coordinates:**")
                    c_x, c_y, c_z = st.columns(3)
                    max_x = max((a["x"] for a in b_state["atoms"]), default=0.0)
                    with c_x:
                        iso_x = st.number_input("X (Å)", value=round(max_x + 1.8, 2), step=0.2)
                    with c_y:
                        iso_y = st.number_input("Y (Å)", value=0.0, step=0.2)
                    with c_z:
                        iso_z = st.number_input("Z (Å)", value=0.0, step=0.2)

                    if st.button("➕ Place Isolated Atom", type="primary", use_container_width=True):
                        new_id = builder.add_atom(b_state, element=selected_element, x=iso_x, y=iso_y, z=iso_z)
                        st.session_state.builder_opt_message = f"Placed isolated {selected_element} #{new_id} at ({iso_x}, {iso_y}, {iso_z})."
                        st.rerun()

        # ---------- TAB 2: CREATE BOND TOOL ----------
        with tab_bond:
            if len(b_state["atoms"]) < 2:
                st.info("Please add at least 2 atoms to create bonds between them.")
            else:
                st.markdown("**Select Two Atoms to Bond:**")
                atom_labels = {a["id"]: f"Atom #{a['id']} ({a['element']}) at ({a['x']:.1f}, {a['y']:.1f}, {a['z']:.1f})" for a in b_state["atoms"]}
                atom_ids = list(atom_labels.keys())

                b_col1, b_col2 = st.columns(2)
                with b_col1:
                    sel_a1 = st.selectbox("Atom 1", atom_ids, index=0, format_func=lambda x: atom_labels[x], key="bond_a1")
                with b_col2:
                    default_a2_idx = 1 if len(atom_ids) > 1 else 0
                    sel_a2 = st.selectbox("Atom 2", atom_ids, index=default_a2_idx, format_func=lambda x: atom_labels[x], key="bond_a2")

                # Check if bond already exists
                existing_bond = builder.get_bond_between(b_state, sel_a1, sel_a2)
                if existing_bond:
                    st.info(f"ℹ️ Atoms #{sel_a1} and #{sel_a2} currently have a **Bond (Order {existing_bond.get('order', 1)})**.")

                new_bond_order = st.radio(
                    "Bond Order:",
                    ["Single (1)", "Double (2)", "Triple (3)"],
                    index=(existing_bond.get("order", 1) - 1) if existing_bond else 0,
                    horizontal=True,
                    key="bond_order_radio"
                )
                b_order_val = 1 if "Single" in new_bond_order else (2 if "Double" in new_bond_order else 3)

                btn_bond_label = "🔄 Update Bond Order" if existing_bond else "🔗 Create Bond"
                if st.button(btn_bond_label, type="primary", use_container_width=True):
                    if sel_a1 == sel_a2:
                        st.error("Cannot create a bond between an atom and itself.")
                    else:
                        builder.create_or_update_bond(b_state, sel_a1, sel_a2, order=b_order_val)
                        st.session_state.builder_opt_message = f"Bond between Atom #{sel_a1} and Atom #{sel_a2} set to Order {b_order_val}."
                        st.rerun()

                if existing_bond:
                    if st.button("❌ Remove Bond", use_container_width=True):
                        builder.delete_bond(b_state, sel_a1, sel_a2)
                        st.session_state.builder_opt_message = f"Removed bond between Atom #{sel_a1} and #{sel_a2}."
                        st.rerun()

        # ---------- TAB 3: EDIT & DELETE TOOL ----------
        with tab_edit:
            if not b_state["atoms"]:
                st.info("No atoms currently exist to edit.")
            else:
                st.markdown("**Inspect / Edit Active Atom:**")
                all_atom_ids = [a["id"] for a in b_state["atoms"]]
                curr_selected = b_state.get("selected_atom_id")
                curr_idx = all_atom_ids.index(curr_selected) if curr_selected in all_atom_ids else 0

                picked_atom_id = st.selectbox(
                    "Select Atom to Edit",
                    all_atom_ids,
                    index=curr_idx,
                    format_func=lambda x: f"Atom #{x} ({builder.get_atom_by_id(b_state, x)['element']})",
                    key="edit_atom_select"
                )
                b_state["selected_atom_id"] = picked_atom_id
                target_atom = builder.get_atom_by_id(b_state, picked_atom_id)

                if target_atom:
                    curr_order_sum = builder.get_atom_bond_order_sum(b_state, picked_atom_id)
                    nbrs = builder.get_atom_neighbors(b_state, picked_atom_id)
                    nbrs_str = ", ".join(f"#{nid} (order {ord_val})" for nid, ord_val in nbrs) if nbrs else "None"

                    st.markdown(f"""
                    <div class="tool-box">
                        <b>Atom #{picked_atom_id} Details:</b><br>
                        • Element: <b>{target_atom['element']}</b><br>
                        • Total Bond Order: <b>{curr_order_sum}</b><br>
                        • Connected Neighbors: <b>{nbrs_str}</b><br>
                        • Coordinates: <code>({target_atom['x']:.3f}, {target_atom['y']:.3f}, {target_atom['z']:.3f})</code>
                    </div>
                    """, unsafe_allow_html=True)

                    # Change element
                    c_elem1, c_elem2 = st.columns([1.5, 1])
                    with c_elem1:
                        new_elem = st.selectbox(
                            "Change Element To:",
                            builder.COMMON_ELEMENTS,
                            index=builder.COMMON_ELEMENTS.index(target_atom["element"]) if target_atom["element"] in builder.COMMON_ELEMENTS else 0,
                            key="change_elem_sel"
                        )
                    with c_elem2:
                        st.write("")
                        if st.button("🔄 Apply Element", use_container_width=True):
                            builder.update_atom_element(b_state, picked_atom_id, new_elem)
                            st.session_state.builder_opt_message = f"Changed Atom #{picked_atom_id} element to {new_elem}."
                            st.rerun()

                    # Delete atom button
                    if st.button(f"🗑️ Delete Atom #{picked_atom_id}", type="secondary", use_container_width=True):
                        builder.delete_atom(b_state, picked_atom_id)
                        st.session_state.builder_opt_message = f"Deleted Atom #{picked_atom_id} and its bonds."
                        st.rerun()

                    # Manual coordinate tweak
                    with st.expander("Adjust Coordinates (X, Y, Z)", expanded=False):
                        cx, cy, cz = st.columns(3)
                        with cx:
                            new_x = st.number_input("X", value=float(target_atom["x"]), step=0.1, key="adj_x")
                        with cy:
                            new_y = st.number_input("Y", value=float(target_atom["y"]), step=0.1, key="adj_y")
                        with cz:
                            new_z = st.number_input("Z", value=float(target_atom["z"]), step=0.1, key="adj_z")
                        if st.button("Save Coordinates", use_container_width=True):
                            builder.update_atom_coordinates(b_state, picked_atom_id, new_x, new_y, new_z)
                            st.session_state.builder_opt_message = f"Updated coordinates for Atom #{picked_atom_id}."
                            st.rerun()

                st.divider()
                st.markdown("**Existing Bonds Management:**")
                if not b_state["bonds"]:
                    st.caption("No bonds in the molecule.")
                else:
                    bond_options = {
                        b["id"]: f"Bond #{b['id']}: Atom #{b['atom1']} ({builder.get_atom_by_id(b_state, b['atom1'])['element']}) — Atom #{b['atom2']} ({builder.get_atom_by_id(b_state, b['atom2'])['element']}) [Order {b.get('order', 1)}]"
                        for b in b_state["bonds"]
                    }
                    sel_bond_id = st.selectbox("Select Bond", list(bond_options.keys()), format_func=lambda x: bond_options[x])
                    target_b = next((b for b in b_state["bonds"] if b["id"] == sel_bond_id), None)
                    if target_b:
                        bo_col1, bo_col2 = st.columns([1.5, 1])
                        with bo_col1:
                            new_b_order = st.selectbox("Change Order", [1, 2, 3], index=target_b.get("order", 1) - 1, format_func=lambda x: f"{['Single', 'Double', 'Triple'][x-1]} ({x})")
                        with bo_col2:
                            st.write("")
                            if st.button("Update Order", use_container_width=True):
                                builder.update_bond_order(b_state, target_b["atom1"], target_b["atom2"], new_b_order)
                                st.session_state.builder_opt_message = f"Updated Bond #{sel_bond_id} to Order {new_b_order}."
                                st.rerun()

                        if st.button(f"🗑️ Delete Bond #{sel_bond_id}", use_container_width=True):
                            builder.delete_bond_by_id(b_state, sel_bond_id)
                            st.session_state.builder_opt_message = f"Deleted Bond #{sel_bond_id}."
                            st.rerun()

        # ---------- TAB 4: 3D STRUCTURE OPTIMIZATION ----------
        with tab_opt:
            st.markdown("**3D Structure Optimization & Minimization:**")
            st.caption("Generates realistic 3D spatial conformations using RDKit's ETKDG conformer generator and UFF energy minimization. If unusual valences or radicals are present, automatically engages the pure-Python force-directed VSEPR relaxation algorithm.")

            opt_add_hs = st.checkbox(
                "Automatically Add Hydrogens to Complete Valences",
                value=False,
                help="Saturates open valences with Hydrogens before 3D energy minimization"
            )

            if st.button("✨ Optimize 3D Structure", type="primary", use_container_width=True):
                success, msg = builder.optimize_3d_geometry(b_state, add_hydrogens=opt_add_hs)
                st.session_state.builder_opt_message = msg
                st.rerun()

            st.write("")
            st.markdown("**Auto-Saturation:**")
            if st.button("➕ Add Hydrogens Only", use_container_width=True, help="Attach explicit Hydrogens to all atoms with open valences"):
                added_h, h_msg = builder.auto_add_hydrogens(b_state)
                st.session_state.builder_opt_message = h_msg
                st.rerun()

    st.write("")

    # ----------------- BOTTOM EXPLORER TABS (ATOM & BOND TABLES) -----------------
    b_tab_atoms, b_tab_bonds, b_tab_comp, b_tab_help = st.tabs([
        "📋 Atom Coordinates Table",
        "🔗 Bond Connectivity Table",
        "📊 Chemical Composition",
        "📖 Avogadro Guide & Chemistry Rules"
    ])

    with b_tab_atoms:
        st.markdown("#### Atom Coordinates")
        if not b_state["atoms"]:
            st.info("No atoms in the workspace.")
        else:
            atom_table_data = []
            for a in b_state["atoms"]:
                v_sum = builder.get_atom_bond_order_sum(b_state, a["id"])
                atom_table_data.append({
                    "ID": a["id"],
                    "Element": a["element"],
                    "X (Å)": round(a["x"], 3),
                    "Y (Å)": round(a["y"], 3),
                    "Z (Å)": round(a["z"], 3),
                    "Valence (Bonds)": v_sum
                })
            st.dataframe(atom_table_data, use_container_width=True, height=280)

    with b_tab_bonds:
        st.markdown("#### Bond Connections")
        if not b_state["bonds"]:
            st.info("No chemical bonds formed yet.")
        else:
            bond_table_data = []
            for b in b_state["bonds"]:
                a1 = builder.get_atom_by_id(b_state, b["atom1"])
                a2 = builder.get_atom_by_id(b_state, b["atom2"])
                order_name = {1: "Single", 2: "Double", 3: "Triple"}.get(b.get("order", 1), "Single")
                bond_table_data.append({
                    "Bond ID": b["id"],
                    "Atom 1": f"#{b['atom1']} ({a1['element'] if a1 else '?'})",
                    "Atom 2": f"#{b['atom2']} ({a2['element'] if a2 else '?'})",
                    "Order": b.get("order", 1),
                    "Type": order_name
                })
            st.dataframe(bond_table_data, use_container_width=True, height=280)

    with b_tab_comp:
        st.markdown("#### Chemical Breakdown")
        if not b_state["atoms"]:
            st.info("Add atoms to view composition.")
        else:
            elem_counts = mol_info["element_counts"]
            total_atoms = mol_info["num_atoms"]
            total_mw = mol_info["molecular_weight"] or 1.0

            comp_rows = []
            for elem, count in sorted(elem_counts.items(), key=lambda x: -x[1]):
                weight_single = parsers.ATOMIC_WEIGHTS.get(elem.upper(), 12.0)
                mass_subtotal = weight_single * count
                atom_pct = (count / total_atoms) * 100.0
                mass_pct = (mass_subtotal / total_mw) * 100.0 if total_mw > 0 else 0.0

                comp_rows.append({
                    "Element": elem,
                    "Count": count,
                    "Atom %": f"{atom_pct:.1f}%",
                    "Atomic Mass (u)": f"{weight_single:.3f}",
                    "Total Mass (g/mol)": f"{mass_subtotal:.2f}",
                    "Mass %": f"{mass_pct:.1f}%"
                })

            c_t1, c_t2 = st.columns([1, 1])
            with c_t1:
                st.table(comp_rows)
            with c_t2:
                st.bar_chart({r["Element"]: r["Count"] for r in comp_rows})

    with b_tab_help:
        st.markdown("""
        ### About the Molecule Builder
        
        The **Molecule Builder** provides an interactive educational workspace inspired by Avogadro.
        
        #### Key Features:
        - **Elemental Palette**: Supports common elements: **H, C, N, O, F, P, S, Cl, Br, I**.
        - **Bond Multiplicity**: Single ($1$), Double ($2$), and Triple ($3$) bonds rendered with distinct multiple cylinders in WebGL.
        - **Valence Rules & Chemistry Validation**:
          - **H**: 1 bond
          - **C**: 4 bonds
          - **N**: 3 bonds (amines) or 4 (ammonium)
          - **O**: 2 bonds
          - **F, Cl, Br, I**: 1 bond (halogens)
          - **P**: 3 or 5 bonds (phosphines, phosphates)
          - **S**: 2, 4, or 6 bonds (thiols, sulfoxides, sulfones)
        - **3D Optimization Engine**:
          - Uses **RDKit ETKDGv3 conformer embedding** and **UFF (Universal Force Field)** energy minimization.
          - Built-in **pure-Python VSEPR force-directed relaxation** fallback prevents crashes on radical or hypervalent structures.
        - **Seamless Export & Inspection**:
          - Export valid **PDB** (with standard CONECT records) and **XYZ** files.
          - Click **🧬 Open in Visualizer** to transfer the built structure into the full Molecule Visualizer!
        """)
