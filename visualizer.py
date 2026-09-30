"""
Py3Dmol-based 3D Molecular Visualizer module.
Builds and customizes interactive WebGL views for PDB and XYZ structures.
"""

from __future__ import annotations
import py3Dmol
from typing import Dict, Any, Optional

# Pre-defined background palette options
BG_PRESETS: Dict[str, str] = {
    "Dark Slate (#0f172a)": "#0f172a",
    "Charcoal (#18181b)": "#18181b",
    "Deep Navy (#0a0f1d)": "#0a0f1d",
    "Pitch Black (#000000)": "#000000",
    "Pure White (#ffffff)": "#ffffff",
    "Light Slate (#f8fafc)": "#f8fafc",
    "Soft Gray (#e2e8f0)": "#e2e8f0",
}

# Color scheme definitions for py3Dmol
COLOR_SCHEMES: Dict[str, Optional[str]] = {
    "CPK Standard (Element)": None,
    "Cyan Carbon": "cyanCarbon",
    "Green Carbon": "greenCarbon",
    "Magenta Carbon": "magentaCarbon",
    "Orange Carbon": "orangeCarbon",
    "Secondary Structure (Spectrum)": "spectrum",
}

def create_molecule_viewer(
    mol_data: Dict[str, Any],
    style: str = "Ball-and-stick",
    bg_color: str = "#0f172a",
    color_scheme: str = "CPK Standard (Element)",
    show_spin: bool = False,
    show_hover: bool = True,
    show_labels: bool = False,
    surface_type: str = "None",
    surface_opacity: float = 0.5,
    viewer_height: int = 580
) -> str:
    """
    Generate an interactive 3D WebGL viewer using Py3Dmol and return self-contained HTML.
    
    Parameters
    ----------
    mol_data : dict
        Parsed molecular data containing 'raw_content', 'file_type', 'atoms', etc.
    style : str
        Rendering style: 'Ball-and-stick', 'Stick', 'Sphere', 'Cartoon', 'Wireframe'.
    bg_color : str
        Hex color for the canvas background.
    color_scheme : str
        Color scheme key from COLOR_SCHEMES.
    show_spin : bool
        Whether to enable continuous auto-rotation.
    show_hover : bool
        Whether to display atom label tooltips on mouse hover.
    show_labels : bool
        Whether to display persistent text labels for atoms.
    surface_type : str
        'None', 'VDW (Van der Waals)', 'SAS (Solvent Accessible)', 'SES (Solvent Excluded)'.
    surface_opacity : float
        Opacity of molecular surface between 0.1 and 1.0.
    viewer_height : int
        Height of the viewer in pixels.
        
    Returns
    -------
    str
        Full HTML markup to be embedded in Streamlit iframe.
    """
    raw_content = mol_data.get("raw_content", "")
    file_type = mol_data.get("file_type", "PDB").lower()
    
    # Initialize py3Dmol view with full container width
    view = py3Dmol.view(width="100%", height=f"{viewer_height}px")
    
    # Add model into the viewer
    view.addModel(raw_content, file_type)
    
    # Determine colorscheme parameter for py3Dmol
    scheme_param = COLOR_SCHEMES.get(color_scheme, None)
    scheme_dict = {"colorscheme": scheme_param} if scheme_param else {}
    
    # Apply rendering styles
    style_lower = style.lower()
    
    if "ball" in style_lower:
        # Ball-and-stick: combination of stick and scaled sphere
        stick_spec = {"radius": 0.14, "multipleBonds": True}
        sphere_spec = {"scale": 0.28}
        if scheme_param:
            stick_spec["colorscheme"] = scheme_param
            sphere_spec["colorscheme"] = scheme_param
        view.setStyle({"stick": stick_spec, "sphere": sphere_spec})
        
    elif "stick" in style_lower:
        # Stick representation
        stick_spec = {"radius": 0.22, "multipleBonds": True}
        if scheme_param:
            stick_spec["colorscheme"] = scheme_param
        view.setStyle({"stick": stick_spec})
        
    elif "sphere" in style_lower or "spacefill" in style_lower:
        # Full Van der Waals spacefill sphere
        sphere_spec = {"scale": 0.82}
        if scheme_param:
            sphere_spec["colorscheme"] = scheme_param
        view.setStyle({"sphere": sphere_spec})
        
    elif "cartoon" in style_lower:
        # Cartoon representation for macromolecules
        if file_type == "pdb":
            cartoon_spec = {"color": "spectrum"} if (scheme_param == "spectrum" or not scheme_param) else {"colorscheme": scheme_param}
            view.setStyle({"cartoon": cartoon_spec})
            # Also render ligands and non-protein heteroatoms as ball-and-stick so they remain visible
            view.addStyle({"hetflag": True}, {"stick": {"radius": 0.16}, "sphere": {"scale": 0.28}})
        else:
            # Fallback for XYZ (which has no secondary structure records)
            stick_spec = {"radius": 0.14}
            sphere_spec = {"scale": 0.28}
            if scheme_param:
                stick_spec["colorscheme"] = scheme_param
                sphere_spec["colorscheme"] = scheme_param
            view.setStyle({"stick": stick_spec, "sphere": sphere_spec})
            
    elif "wireframe" in style_lower or "line" in style_lower:
        # Line / wireframe representation
        view.setStyle({"line": {"linewidth": 2}})
        
    else:
        # Default fallback
        view.setStyle({"stick": {"radius": 0.15}, "sphere": {"scale": 0.25}})
        
    # Background color
    view.setBackgroundColor(bg_color)
    
    # Persistent atom labels if requested
    if show_labels:
        # Show atom symbols
        view.addPropertyLabels("elem", {}, {
            "fontColor": "#ffffff" if bg_color.lower() in ("#0f172a", "#000000", "#18181b", "#0a0f1d") else "#1e293b",
            "fontSize": 10,
            "showBackground": True,
            "backgroundColor": "#1e293b",
            "backgroundOpacity": 0.6,
            "alignment": "center"
        })
        
    # Hover tooltips
    if show_hover:
        hover_js = """function(atom, viewer) {
            if (!atom.label) {
                var txt = (atom.elem || '') + (atom.serial ? (' #' + atom.serial) : '');
                if (atom.resn) txt += ' (' + atom.resn + (atom.resi ? atom.resi : '') + ')';
                atom.label = viewer.addLabel(txt, {
                    position: atom,
                    backgroundColor: '#0f172a',
                    backgroundOpacity: 0.88,
                    fontColor: '#38bdf8',
                    fontSize: 12,
                    borderThickness: 1,
                    borderColor: '#38bdf8'
                });
            }
        }"""
        unhover_js = """function(atom, viewer) {
            if (atom.label) {
                viewer.removeLabel(atom.label);
                delete atom.label;
            }
        }"""
        view.setHoverable({}, True, hover_js, unhover_js)
        
    # Molecular surface overlay
    if surface_type != "None":
        surf_map = {
            "VDW (Van der Waals)": 1,
            "SAS (Solvent Accessible)": 2,
            "SES (Solvent Excluded)": 3
        }
        stype_val = surf_map.get(surface_type, 1)
        surf_color = "#38bdf8" if "dark" in bg_color or bg_color in ("#0f172a", "#000000", "#18181b", "#0a0f1d") else "#64748b"
        view.addSurface(stype_val, {"opacity": float(surface_opacity), "color": surf_color})
        
    # Center and scale to fill viewport
    view.zoomTo()
    
    # Auto-rotation if enabled
    if show_spin:
        view.spin(True, 1)
        
    # Generate standard py3Dmol HTML
    base_html = view._make_html()
    uniqueid = view.uniqueid
    
    # Enhance base HTML with an integrated, floating HUD control bar
    # and custom styles for a modern visual interface
    is_dark = bg_color.lower() in ("#0f172a", "#18181b", "#0a0f1d", "#000000", "#111827", "#1e1e2f")
    hud_bg = "rgba(15, 23, 42, 0.75)" if is_dark else "rgba(255, 255, 255, 0.85)"
    hud_border = "rgba(56, 189, 248, 0.3)" if is_dark else "rgba(203, 213, 225, 0.8)"
    btn_color = "#f1f5f9" if is_dark else "#1e293b"
    btn_hover = "rgba(56, 189, 248, 0.25)" if is_dark else "rgba(56, 189, 248, 0.15)"
    
    hud_html = f"""
    <style>
        body, html {{
            margin: 0;
            padding: 0;
            overflow: hidden;
            width: 100%;
            height: 100%;
            background-color: {bg_color};
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        }}
        .hud-toolbar {{
            position: absolute;
            top: 14px;
            right: 14px;
            z-index: 1000;
            display: flex;
            align-items: center;
            gap: 6px;
            background: {hud_bg};
            backdrop-filter: blur(8px);
            -webkit-backdrop-filter: blur(8px);
            padding: 5px 8px;
            border-radius: 10px;
            border: 1px solid {hud_border};
            box-shadow: 0 4px 16px rgba(0, 0, 0, 0.25);
        }}
        .hud-btn {{
            background: transparent;
            border: none;
            color: {btn_color};
            cursor: pointer;
            padding: 6px 10px;
            border-radius: 6px;
            font-size: 13px;
            font-weight: 500;
            display: flex;
            align-items: center;
            gap: 5px;
            transition: all 0.15s ease;
        }}
        .hud-btn:hover {{
            background: {btn_hover};
            color: #38bdf8;
            transform: translateY(-1px);
        }}
        .hud-btn:active {{
            transform: translateY(0);
        }}
        .hud-badge {{
            position: absolute;
            bottom: 14px;
            left: 14px;
            z-index: 1000;
            background: {hud_bg};
            backdrop-filter: blur(8px);
            border: 1px solid {hud_border};
            border-radius: 8px;
            padding: 4px 10px;
            font-size: 11px;
            color: {btn_color};
            opacity: 0.85;
            display: flex;
            gap: 12px;
            pointer-events: none;
        }}
    </style>
    
    <div class="hud-toolbar">
        <button class="hud-btn" id="btn-reset" title="Reset View Angle & Zoom" onclick="reset3DView()">
            <span>🔄</span> Reset
        </button>
        <button class="hud-btn" id="btn-spin" title="Toggle Auto-Spin Animation" onclick="toggle3DSpin()">
            <span>✨</span> Spin
        </button>
        <button class="hud-btn" id="btn-zoomin" title="Zoom In" onclick="zoomIn3D()">
            <span>➕</span>
        </button>
        <button class="hud-btn" id="btn-zoomout" title="Zoom Out" onclick="zoomOut3D()">
            <span>➖</span>
        </button>
        <button class="hud-btn" id="btn-screenshot" title="Download Screenshot PNG" onclick="export3DPNG()">
            <span>📸</span> PNG
        </button>
    </div>
    
    <div class="hud-badge">
        <span>🖱️ <b>Left Click:</b> Rotate</span>
        <span>✋ <b>Right Click / Shift:</b> Pan</span>
        <span>🔍 <b>Scroll:</b> Zoom</span>
    </div>
    
    <script>
    var spinningState = {"true" if show_spin else "false"};
    
    function reset3DView() {{
        if (typeof viewer_{uniqueid} !== 'undefined') {{
            viewer_{uniqueid}.zoomTo(600);
            viewer_{uniqueid}.render();
        }}
    }}
    
    function toggle3DSpin() {{
        if (typeof viewer_{uniqueid} !== 'undefined') {{
            spinningState = !spinningState;
            viewer_{uniqueid}.spin(spinningState, 1);
        }}
    }}
    
    function zoomIn3D() {{
        if (typeof viewer_{uniqueid} !== 'undefined') {{
            viewer_{uniqueid}.zoom(1.25);
            viewer_{uniqueid}.render();
        }}
    }}
    
    function zoomOut3D() {{
        if (typeof viewer_{uniqueid} !== 'undefined') {{
            viewer_{uniqueid}.zoom(0.8);
            viewer_{uniqueid}.render();
        }}
    }}
    
    function export3DPNG() {{
        if (typeof viewer_{uniqueid} !== 'undefined') {{
            var uri = viewer_{uniqueid}.pngURI();
            var link = document.createElement('a');
            link.download = 'molecule_{uniqueid}.png';
            link.href = uri;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}
    }}
    </script>
    """
    
    # Combine py3Dmol base HTML and the HUD overlay
    full_html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
</head>
<body>
    {base_html}
    {hud_html}
</body>
</html>
"""
    return full_html
