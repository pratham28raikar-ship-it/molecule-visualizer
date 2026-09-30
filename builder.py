"""
Molecule Builder Module
Provides interactive 3D molecule building, editing, chemistry valence validation,
structure optimization (RDKit with pure-Python VSEPR force-directed fallback),
and export capabilities (PDB, XYZ, V2000 Molfile).
"""

from __future__ import annotations
import math
import random
import copy
from typing import Dict, List, Tuple, Any, Optional, Set

import parsers

# Common elements for builder with standard properties
COMMON_ELEMENTS: List[str] = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"]

# Standard valence rules (typical number of covalent bonds)
ELEMENT_VALENCES: Dict[str, List[int]] = {
    "H": [1],
    "C": [4],
    "N": [3, 4],     # 3 typical (amines), 4 in ammonium
    "O": [2],        # 2 typical, 3 in hydronium/oxonium
    "F": [1],
    "CL": [1, 3, 5, 7],
    "BR": [1, 3, 5],
    "I": [1, 3, 5, 7],
    "P": [3, 5],     # Phosphine, Phosphates
    "S": [2, 4, 6],  # Thiol/sulfide, Sulfoxide, Sulfone/sulfate
}

# Default covalent bond lengths (Angstroms) for standard single bonds
DEFAULT_BOND_LENGTHS: Dict[str, float] = {
    "H": 0.37, "C": 0.77, "N": 0.74, "O": 0.73, "F": 0.71,
    "P": 1.10, "S": 1.03, "CL": 0.99, "BR": 1.14, "I": 1.33
}

# Van der Waals radii (Angstroms)
VDW_RADII: Dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98
}

# Pre-defined molecular templates
TEMPLATES: Dict[str, Dict[str, Any]] = {
    "Water (H2O)": {
        "atoms": [
            {"id": 1, "element": "O", "x": 0.0, "y": 0.0, "z": 0.0},
            {"id": 2, "element": "H", "x": 0.757, "y": 0.586, "z": 0.0},
            {"id": 3, "element": "H", "x": -0.757, "y": 0.586, "z": 0.0},
        ],
        "bonds": [
            {"id": 1, "atom1": 1, "atom2": 2, "order": 1},
            {"id": 2, "atom1": 1, "atom2": 3, "order": 1},
        ]
    },
    "Methane (CH4)": {
        "atoms": [
            {"id": 1, "element": "C", "x": 0.0, "y": 0.0, "z": 0.0},
            {"id": 2, "element": "H", "x": 0.629, "y": 0.629, "z": 0.629},
            {"id": 3, "element": "H", "x": -0.629, "y": -0.629, "z": 0.629},
            {"id": 4, "element": "H", "x": -0.629, "y": 0.629, "z": -0.629},
            {"id": 5, "element": "H", "x": 0.629, "y": -0.629, "z": -0.629},
        ],
        "bonds": [
            {"id": 1, "atom1": 1, "atom2": 2, "order": 1},
            {"id": 2, "atom1": 1, "atom2": 3, "order": 1},
            {"id": 3, "atom1": 1, "atom2": 4, "order": 1},
            {"id": 4, "atom1": 1, "atom2": 5, "order": 1},
        ]
    },
    "Ethylene (C2H4)": {
        "atoms": [
            {"id": 1, "element": "C", "x": -0.665, "y": 0.0, "z": 0.0},
            {"id": 2, "element": "C", "x": 0.665, "y": 0.0, "z": 0.0},
            {"id": 3, "element": "H", "x": -1.23, "y": 0.928, "z": 0.0},
            {"id": 4, "element": "H", "x": -1.23, "y": -0.928, "z": 0.0},
            {"id": 5, "element": "H", "x": 1.23, "y": 0.928, "z": 0.0},
            {"id": 6, "element": "H", "x": 1.23, "y": -0.928, "z": 0.0},
        ],
        "bonds": [
            {"id": 1, "atom1": 1, "atom2": 2, "order": 2},
            {"id": 2, "atom1": 1, "atom2": 3, "order": 1},
            {"id": 3, "atom1": 1, "atom2": 4, "order": 1},
            {"id": 4, "atom1": 2, "atom2": 5, "order": 1},
            {"id": 5, "atom1": 2, "atom2": 6, "order": 1},
        ]
    },
    "Acetylene (C2H2)": {
        "atoms": [
            {"id": 1, "element": "C", "x": -0.60, "y": 0.0, "z": 0.0},
            {"id": 2, "element": "C", "x": 0.60, "y": 0.0, "z": 0.0},
            {"id": 3, "element": "H", "x": -1.66, "y": 0.0, "z": 0.0},
            {"id": 4, "element": "H", "x": 1.66, "y": 0.0, "z": 0.0},
        ],
        "bonds": [
            {"id": 1, "atom1": 1, "atom2": 2, "order": 3},
            {"id": 2, "atom1": 1, "atom2": 3, "order": 1},
            {"id": 3, "atom1": 2, "atom2": 4, "order": 1},
        ]
    },
    "Ethanol (C2H6O)": {
        "atoms": [
            {"id": 1, "element": "C", "x": -1.18, "y": -0.40, "z": 0.0},
            {"id": 2, "element": "C", "x": 0.0, "y": 0.54, "z": 0.0},
            {"id": 3, "element": "O", "x": 1.21, "y": -0.22, "z": 0.0},
            {"id": 4, "element": "H", "x": 1.95, "y": 0.39, "z": 0.0},
            {"id": 5, "element": "H", "x": -1.18, "y": -1.04, "z": 0.89},
            {"id": 6, "element": "H", "x": -1.18, "y": -1.04, "z": -0.89},
            {"id": 7, "element": "H", "x": -2.09, "y": 0.20, "z": 0.0},
            {"id": 8, "element": "H", "x": 0.03, "y": 1.18, "z": 0.89},
            {"id": 9, "element": "H", "x": 0.03, "y": 1.18, "z": -0.89},
        ],
        "bonds": [
            {"id": 1, "atom1": 1, "atom2": 2, "order": 1},
            {"id": 2, "atom1": 2, "atom2": 3, "order": 1},
            {"id": 3, "atom1": 3, "atom2": 4, "order": 1},
            {"id": 4, "atom1": 1, "atom2": 5, "order": 1},
            {"id": 5, "atom1": 1, "atom2": 6, "order": 1},
            {"id": 6, "atom1": 1, "atom2": 7, "order": 1},
            {"id": 7, "atom1": 2, "atom2": 8, "order": 1},
            {"id": 8, "atom1": 2, "atom2": 9, "order": 1},
        ]
    },
    "Benzene (C6H6)": {
        "atoms": [
            {"id": 1, "element": "C", "x": 1.397, "y": 0.0, "z": 0.0},
            {"id": 2, "element": "C", "x": 0.698, "y": 1.21, "z": 0.0},
            {"id": 3, "element": "C", "x": -0.698, "y": 1.21, "z": 0.0},
            {"id": 4, "element": "C", "x": -1.397, "y": 0.0, "z": 0.0},
            {"id": 5, "element": "C", "x": -0.698, "y": -1.21, "z": 0.0},
            {"id": 6, "element": "C", "x": 0.698, "y": -1.21, "z": 0.0},
            {"id": 7, "element": "H", "x": 2.48, "y": 0.0, "z": 0.0},
            {"id": 8, "element": "H", "x": 1.24, "y": 2.15, "z": 0.0},
            {"id": 9, "element": "H", "x": -1.24, "y": 2.15, "z": 0.0},
            {"id": 10, "element": "H", "x": -2.48, "y": 0.0, "z": 0.0},
            {"id": 11, "element": "H", "x": -1.24, "y": -2.15, "z": 0.0},
            {"id": 12, "element": "H", "x": 1.24, "y": -2.15, "z": 0.0},
        ],
        "bonds": [
            {"id": 1, "atom1": 1, "atom2": 2, "order": 2},
            {"id": 2, "atom1": 2, "atom2": 3, "order": 1},
            {"id": 3, "atom1": 3, "atom2": 4, "order": 2},
            {"id": 4, "atom1": 4, "atom2": 5, "order": 1},
            {"id": 5, "atom1": 5, "atom2": 6, "order": 2},
            {"id": 6, "atom1": 6, "atom2": 1, "order": 1},
            {"id": 7, "atom1": 1, "atom2": 7, "order": 1},
            {"id": 8, "atom1": 2, "atom2": 8, "order": 1},
            {"id": 9, "atom1": 3, "atom2": 9, "order": 1},
            {"id": 10, "atom1": 4, "atom2": 10, "order": 1},
            {"id": 11, "atom1": 5, "atom2": 11, "order": 1},
            {"id": 12, "atom1": 6, "atom2": 12, "order": 1},
        ]
    }
}


def create_initial_state() -> Dict[str, Any]:
    """Create a fresh builder state dictionary."""
    return {
        "atoms": [],
        "bonds": [],
        "next_atom_id": 1,
        "next_bond_id": 1,
        "selected_atom_id": None,
        "history": [],
        "redo_history": []
    }


def push_history(state: Dict[str, Any]) -> None:
    """Save current state snapshot to history for Undo (keeps up to 35 steps)."""
    snapshot = {
        "atoms": copy.deepcopy(state["atoms"]),
        "bonds": copy.deepcopy(state["bonds"]),
        "next_atom_id": state["next_atom_id"],
        "next_bond_id": state["next_bond_id"],
        "selected_atom_id": state["selected_atom_id"],
    }
    state["history"].append(snapshot)
    if len(state["history"]) > 35:
        state["history"].pop(0)
    # Clear redo stack on new action
    state["redo_history"].clear()


def undo(state: Dict[str, Any]) -> bool:
    """Revert state to previous snapshot in history."""
    if not state["history"]:
        return False
    current_snapshot = {
        "atoms": copy.deepcopy(state["atoms"]),
        "bonds": copy.deepcopy(state["bonds"]),
        "next_atom_id": state["next_atom_id"],
        "next_bond_id": state["next_bond_id"],
        "selected_atom_id": state["selected_atom_id"],
    }
    state["redo_history"].append(current_snapshot)
    previous = state["history"].pop()
    state["atoms"] = previous["atoms"]
    state["bonds"] = previous["bonds"]
    state["next_atom_id"] = previous["next_atom_id"]
    state["next_bond_id"] = previous["next_bond_id"]
    state["selected_atom_id"] = previous["selected_atom_id"]
    return True


def redo(state: Dict[str, Any]) -> bool:
    """Redo previously undone action."""
    if not state["redo_history"]:
        return False
    current_snapshot = {
        "atoms": copy.deepcopy(state["atoms"]),
        "bonds": copy.deepcopy(state["bonds"]),
        "next_atom_id": state["next_atom_id"],
        "next_bond_id": state["next_bond_id"],
        "selected_atom_id": state["selected_atom_id"],
    }
    state["history"].append(current_snapshot)
    future = state["redo_history"].pop()
    state["atoms"] = future["atoms"]
    state["bonds"] = future["bonds"]
    state["next_atom_id"] = future["next_atom_id"]
    state["next_bond_id"] = future["next_bond_id"]
    state["selected_atom_id"] = future["selected_atom_id"]
    return True


def clear_molecule(state: Dict[str, Any]) -> None:
    """Remove all atoms and bonds, saving snapshot to history."""
    if not state["atoms"] and not state["bonds"]:
        return
    push_history(state)
    state["atoms"].clear()
    state["bonds"].clear()
    state["next_atom_id"] = 1
    state["next_bond_id"] = 1
    state["selected_atom_id"] = None


def load_template_molecule(state: Dict[str, Any], template_name: str) -> bool:
    """Load a molecular template into state."""
    if template_name not in TEMPLATES:
        return False
    push_history(state)
    tmpl = TEMPLATES[template_name]
    state["atoms"] = copy.deepcopy(tmpl["atoms"])
    state["bonds"] = copy.deepcopy(tmpl["bonds"])
    state["next_atom_id"] = max((a["id"] for a in state["atoms"]), default=0) + 1
    state["next_bond_id"] = max((b["id"] for b in state["bonds"]), default=0) + 1
    state["selected_atom_id"] = state["atoms"][0]["id"] if state["atoms"] else None
    return True


def get_atom_by_id(state: Dict[str, Any], atom_id: int) -> Optional[Dict[str, Any]]:
    """Retrieve atom dict by ID."""
    for atom in state["atoms"]:
        if atom["id"] == atom_id:
            return atom
    return None


def get_bond_between(state: Dict[str, Any], atom1_id: int, atom2_id: int) -> Optional[Dict[str, Any]]:
    """Find existing bond between two atoms."""
    for b in state["bonds"]:
        if (b["atom1"] == atom1_id and b["atom2"] == atom2_id) or \
           (b["atom1"] == atom2_id and b["atom2"] == atom1_id):
            return b
    return None


def calculate_attached_position(state: Dict[str, Any], parent_atom_id: int, new_element: str, bond_order: int = 1) -> Tuple[float, float, float]:
    """
    Calculate intelligent 3D position for a new atom attached to an existing atom.
    Determines an open direction based on current bonds and typical bond lengths.
    """
    parent = get_atom_by_id(state, parent_atom_id)
    if not parent:
        return (0.0, 0.0, 0.0)

    # Calculate typical bond length based on covalent radii
    p_elem = parent["element"].upper()
    n_elem = new_element.upper()
    r1 = DEFAULT_BOND_LENGTHS.get(p_elem, 0.77)
    r2 = DEFAULT_BOND_LENGTHS.get(n_elem, 0.77)
    bond_dist = r1 + r2
    if bond_order == 2:
        bond_dist *= 0.88
    elif bond_order == 3:
        bond_dist *= 0.78

    # Find vectors of already attached neighbors
    neighbor_vectors = []
    for b in state["bonds"]:
        nbr_id = None
        if b["atom1"] == parent_atom_id:
            nbr_id = b["atom2"]
        elif b["atom2"] == parent_atom_id:
            nbr_id = b["atom1"]
        if nbr_id is not None:
            nbr = get_atom_by_id(state, nbr_id)
            if nbr:
                dx = nbr["x"] - parent["x"]
                dy = nbr["y"] - parent["y"]
                dz = nbr["z"] - parent["z"]
                norm = math.sqrt(dx*dx + dy*dy + dz*dz)
                if norm > 1e-4:
                    neighbor_vectors.append((dx / norm, dy / norm, dz / norm))

    num_nbrs = len(neighbor_vectors)
    if num_nbrs == 0:
        # First attached atom: along X-axis
        direction = (1.0, 0.0, 0.0)
    elif num_nbrs == 1:
        # Linear opposite or tetrahedral tilt:
        # Opposite direction with slight tilt for tetrahedral/trigonal
        vx, vy, vz = neighbor_vectors[0]
        # Choose a perpendicular vector
        if abs(vx) < 0.9:
            px, py, pz = 1.0, 0.0, 0.0
        else:
            px, py, pz = 0.0, 1.0, 0.0
        # Gram-Schmidt perpendicular
        dot = vx*px + vy*py + vz*pz
        cx, cy, cz = px - dot*vx, py - dot*vy, pz - dot*vz
        c_norm = math.sqrt(cx*cx + cy*cy + cz*cz)
        if c_norm > 1e-4:
            cx, cy, cz = cx / c_norm, cy / c_norm, cz / c_norm
        # Standard tetrahedral angle: cos(109.5 deg) ~ -0.334
        # Direction = -0.334 * v + 0.942 * perp
        dx = -0.334 * vx + 0.942 * cx
        dy = -0.334 * vy + 0.942 * cy
        dz = -0.334 * vz + 0.942 * cz
        d_norm = math.sqrt(dx*dx + dy*dy + dz*dz)
        direction = (dx / d_norm, dy / d_norm, dz / d_norm) if d_norm > 1e-4 else (-vx, -vy, -vz)
    elif num_nbrs == 2:
        # Between the two neighbors, pointed outwards
        v1, v2 = neighbor_vectors[0], neighbor_vectors[1]
        sum_x, sum_y, sum_z = v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]
        s_norm = math.sqrt(sum_x*sum_x + sum_y*sum_y + sum_z*sum_z)
        if s_norm > 1e-4:
            # Negative sum vector points away from both neighbors
            # Add a slight out-of-plane z-component to form a tripod/tetrahedron
            cross_x = v1[1]*v2[2] - v1[2]*v2[1]
            cross_y = v1[2]*v2[0] - v1[0]*v2[2]
            cross_z = v1[0]*v2[1] - v1[1]*v2[0]
            cr_norm = math.sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z)
            if cr_norm > 1e-4:
                nx, ny, nz = -sum_x / s_norm, -sum_y / s_norm, -sum_z / s_norm
                px, py, pz = cross_x / cr_norm, cross_y / cr_norm, cross_z / cr_norm
                dx = 0.8 * nx + 0.6 * px
                dy = 0.8 * ny + 0.6 * py
                dz = 0.8 * nz + 0.6 * pz
                d_norm = math.sqrt(dx*dx + dy*dy + dz*dz)
                direction = (dx / d_norm, dy / d_norm, dz / d_norm)
            else:
                direction = (-sum_x / s_norm, -sum_y / s_norm, -sum_z / s_norm)
        else:
            # Opposite colinear
            direction = (0.0, 1.0, 0.0)
    elif num_nbrs == 3:
        # Tripod apex: negative sum of all 3 vectors
        sx = sum(v[0] for v in neighbor_vectors)
        sy = sum(v[1] for v in neighbor_vectors)
        sz = sum(v[2] for v in neighbor_vectors)
        s_norm = math.sqrt(sx*sx + sy*sy + sz*sz)
        if s_norm > 1e-4:
            direction = (-sx / s_norm, -sy / s_norm, -sz / s_norm)
        else:
            direction = (0.0, 0.0, 1.0)
    else:
        # 4 or more: pick random open direction
        theta = random.uniform(0, 2 * math.pi)
        phi = random.uniform(-math.pi/2, math.pi/2)
        direction = (math.cos(phi) * math.cos(theta), math.cos(phi) * math.sin(theta), math.sin(phi))

    return (
        round(parent["x"] + direction[0] * bond_dist, 3),
        round(parent["y"] + direction[1] * bond_dist, 3),
        round(parent["z"] + direction[2] * bond_dist, 3),
    )


def add_atom(
    state: Dict[str, Any],
    element: str,
    x: Optional[float] = None,
    y: Optional[float] = None,
    z: Optional[float] = None,
    attach_to_id: Optional[int] = None,
    bond_order: int = 1
) -> int:
    """
    Add a new atom to the workspace.
    If attach_to_id is provided, automatically positions and bonds the new atom.
    """
    push_history(state)
    atom_id = state["next_atom_id"]
    state["next_atom_id"] += 1
    norm_elem = parsers.clean_element_symbol(element)

    if attach_to_id is not None and get_atom_by_id(state, attach_to_id):
        new_x, new_y, new_z = calculate_attached_position(state, attach_to_id, norm_elem, bond_order)
    elif x is not None and y is not None and z is not None:
        new_x, new_y, new_z = float(x), float(y), float(z)
    else:
        # Automatic isolated position
        if not state["atoms"]:
            new_x, new_y, new_z = 0.0, 0.0, 0.0
        else:
            # Place at comfortable offset from existing atoms
            max_x = max(a["x"] for a in state["atoms"])
            new_x, new_y, new_z = round(max_x + 1.8, 3), 0.0, 0.0

    atom_record = {
        "id": atom_id,
        "element": norm_elem,
        "x": new_x,
        "y": new_y,
        "z": new_z
    }
    state["atoms"].append(atom_record)
    state["selected_atom_id"] = atom_id

    # If attaching to existing atom, create bond
    if attach_to_id is not None and get_atom_by_id(state, attach_to_id):
        bond_id = state["next_bond_id"]
        state["next_bond_id"] += 1
        state["bonds"].append({
            "id": bond_id,
            "atom1": attach_to_id,
            "atom2": atom_id,
            "order": max(1, min(3, int(bond_order)))
        })

    return atom_id


def create_or_update_bond(state: Dict[str, Any], atom1_id: int, atom2_id: int, order: int = 1) -> Optional[int]:
    """
    Create a new bond or update the order of an existing bond between atom1 and atom2.
    """
    if atom1_id == atom2_id:
        return None
    a1 = get_atom_by_id(state, atom1_id)
    a2 = get_atom_by_id(state, atom2_id)
    if not a1 or not a2:
        return None

    push_history(state)
    order = max(1, min(3, int(order)))
    existing = get_bond_between(state, atom1_id, atom2_id)
    if existing:
        existing["order"] = order
        return existing["id"]

    bond_id = state["next_bond_id"]
    state["next_bond_id"] += 1
    state["bonds"].append({
        "id": bond_id,
        "atom1": min(atom1_id, atom2_id),
        "atom2": max(atom1_id, atom2_id),
        "order": order
    })
    return bond_id


def delete_atom(state: Dict[str, Any], atom_id: int) -> bool:
    """Delete an atom and remove all bonds attached to it."""
    atom = get_atom_by_id(state, atom_id)
    if not atom:
        return False

    push_history(state)
    state["atoms"] = [a for a in state["atoms"] if a["id"] != atom_id]
    state["bonds"] = [b for b in state["bonds"] if b["atom1"] != atom_id and b["atom2"] != atom_id]

    if state["selected_atom_id"] == atom_id:
        state["selected_atom_id"] = state["atoms"][0]["id"] if state["atoms"] else None
    return True


def delete_bond(state: Dict[str, Any], atom1_id: int, atom2_id: int) -> bool:
    """Delete the bond connecting two atoms."""
    existing = get_bond_between(state, atom1_id, atom2_id)
    if not existing:
        return False
    push_history(state)
    state["bonds"] = [b for b in state["bonds"] if b["id"] != existing["id"]]
    return True


def delete_bond_by_id(state: Dict[str, Any], bond_id: int) -> bool:
    """Delete a bond by its unique ID."""
    found = any(b["id"] == bond_id for b in state["bonds"])
    if not found:
        return False
    push_history(state)
    state["bonds"] = [b for b in state["bonds"] if b["id"] != bond_id]
    return True


def update_atom_element(state: Dict[str, Any], atom_id: int, new_element: str) -> bool:
    """Change the chemical element of an atom."""
    atom = get_atom_by_id(state, atom_id)
    if not atom:
        return False
    push_history(state)
    atom["element"] = parsers.clean_element_symbol(new_element)
    return True


def update_atom_coordinates(state: Dict[str, Any], atom_id: int, x: float, y: float, z: float) -> bool:
    """Update 3D Cartesian coordinates of an atom."""
    atom = get_atom_by_id(state, atom_id)
    if not atom:
        return False
    push_history(state)
    atom["x"] = round(float(x), 4)
    atom["y"] = round(float(y), 4)
    atom["z"] = round(float(z), 4)
    return True


def update_bond_order(state: Dict[str, Any], atom1_id: int, atom2_id: int, new_order: int) -> bool:
    """Change the bond order (1, 2, or 3) between two atoms."""
    existing = get_bond_between(state, atom1_id, atom2_id)
    if not existing:
        return False
    push_history(state)
    existing["order"] = max(1, min(3, int(new_order)))
    return True


# =========================================================================
# CHEMISTRY & VALENCE VALIDATION
# =========================================================================

def get_atom_bond_order_sum(state: Dict[str, Any], atom_id: int) -> int:
    """Compute total bond order for a given atom (e.g. single=1, double=2, triple=3)."""
    total = 0
    for b in state["bonds"]:
        if b["atom1"] == atom_id or b["atom2"] == atom_id:
            total += b.get("order", 1)
    return total


def get_atom_neighbors(state: Dict[str, Any], atom_id: int) -> List[Tuple[int, int]]:
    """Return list of (neighbor_atom_id, bond_order) for an atom."""
    nbrs = []
    for b in state["bonds"]:
        if b["atom1"] == atom_id:
            nbrs.append((b["atom2"], b.get("order", 1)))
        elif b["atom2"] == atom_id:
            nbrs.append((b["atom1"], b.get("order", 1)))
    return nbrs


def validate_valences(state: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Perform chemistry validation on each atom's valence.
    Returns list of validation feedback dicts:
    {
        'atom_id': int,
        'element': str,
        'current_valence': int,
        'allowed_valences': List[int],
        'status': 'normal' | 'warning' | 'info',
        'message': str
    }
    """
    results = []
    for atom in state["atoms"]:
        aid = atom["id"]
        elem = atom["element"].upper()
        current = get_atom_bond_order_sum(state, aid)
        allowed = ELEMENT_VALENCES.get(elem, [1, 2, 3, 4, 5, 6])
        max_allowed = max(allowed)
        min_allowed = min(allowed)

        if current > max_allowed:
            results.append({
                "atom_id": aid,
                "element": atom["element"],
                "current_valence": current,
                "allowed_valences": allowed,
                "status": "warning",
                "message": f"⚠️ **Valence Warning:** {atom['element']} (Atom #{aid}) currently has {current} bonds, which exceeds its typical valence of {max_allowed}."
            })
        elif current < min_allowed and len(state["atoms"]) > 1:
            diff = min_allowed - current
            h_hint = f" (Needs {diff} more bond{'s' if diff > 1 else ''} or hydrogen{'s' if diff > 1 else ''})"
            results.append({
                "atom_id": aid,
                "element": atom["element"],
                "current_valence": current,
                "allowed_valences": allowed,
                "status": "info",
                "message": f"ℹ️ **Open Valence:** {atom['element']} (Atom #{aid}) has bond order {current} (typical: {allowed[0]}).{h_hint}"
            })
    return results


def compute_molecular_info(state: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate molecular formula, weight, atom & bond counts, and composition."""
    num_atoms = len(state["atoms"])
    num_bonds = len(state["bonds"])
    single_bonds = sum(1 for b in state["bonds"] if b.get("order", 1) == 1)
    double_bonds = sum(1 for b in state["bonds"] if b.get("order", 1) == 2)
    triple_bonds = sum(1 for b in state["bonds"] if b.get("order", 1) == 3)

    element_counts: Dict[str, int] = {}
    for a in state["atoms"]:
        elem = a["element"]
        element_counts[elem] = element_counts.get(elem, 0) + 1

    # Molecular weight calculation
    mw = sum(parsers.ATOMIC_WEIGHTS.get(elem.upper(), 12.0) * count for elem, count in element_counts.items())

    # Hill system formula: C first, then H, then alphabetical
    sorted_elements = sorted(
        element_counts.keys(),
        key=lambda e: (0 if e.upper() == "C" else (1 if e.upper() == "H" else 2), e.upper())
    )
    formula_parts = []
    for elem in sorted_elements:
        c = element_counts[elem]
        formula_parts.append(f"{elem}{c if c > 1 else ''}")
    formula = "".join(formula_parts) if formula_parts else "Empty"

    return {
        "formula": formula,
        "molecular_weight": round(mw, 3),
        "num_atoms": num_atoms,
        "num_bonds": num_bonds,
        "single_bonds": single_bonds,
        "double_bonds": double_bonds,
        "triple_bonds": triple_bonds,
        "element_counts": element_counts,
        "elements_present": list(element_counts.keys())
    }


# =========================================================================
# EXPORT FORMAT GENERATORS (PDB, XYZ, V2000 MOL)
# =========================================================================

def export_to_xyz(state: Dict[str, Any], title: str = "Molecule Builder Model") -> str:
    """
    Generate standard XYZ molecular format string.
    Line 1: Atom count
    Line 2: Title / comments
    Lines 3..: Element X Y Z
    """
    lines = [str(len(state["atoms"])), title]
    for a in state["atoms"]:
        lines.append(f"{a['element']:<3} {a['x']:>12.4f} {a['y']:>12.4f} {a['z']:>12.4f}")
    return "\n".join(lines) + "\n"


def export_to_pdb(state: Dict[str, Any], title: str = "Molecule Builder Model") -> str:
    """
    Generate standard PDB molecular file with ATOM/HETATM and CONECT records.
    Repeats CONECT target entries for double/triple bonds so Py3Dmol renders bond multiplicity.
    """
    lines = [f"HEADER    {title[:40]:<40}", "COMPND    MOL"]
    
    # Map internal IDs to 1-indexed PDB serials
    id_to_serial: Dict[int, int] = {}
    for idx, a in enumerate(state["atoms"], start=1):
        id_to_serial[a["id"]] = idx
        elem = a["element"]
        atom_name = f"{elem[:2]}{idx}"
        # PDB format fixed columns
        lines.append(
            f"HETATM{idx:>5} {atom_name:<4} UNL     1    "
            f"{a['x']:>8.3f}{a['y']:>8.3f}{a['z']:>8.3f}"
            f"  1.00  0.00          {elem:>2}  "
        )

    # CONECT records
    conect_dict: Dict[int, List[int]] = {idx: [] for idx in range(1, len(state["atoms"]) + 1)}
    for b in state["bonds"]:
        s1 = id_to_serial.get(b["atom1"])
        s2 = id_to_serial.get(b["atom2"])
        if s1 and s2:
            order = b.get("order", 1)
            for _ in range(order):
                conect_dict[s1].append(s2)
                conect_dict[s2].append(s1)

    for src in sorted(conect_dict.keys()):
        targets = conect_dict[src]
        if targets:
            # Chunk into lines of up to 4 targets per PDB CONECT specification
            for chunk_start in range(0, len(targets), 4):
                chunk = targets[chunk_start:chunk_start + 4]
                t_str = "".join(f"{t:>5}" for t in chunk)
                lines.append(f"CONECT{src:>5}{t_str}")

    lines.append("END\n")
    return "\n".join(lines)


def export_to_mol_block(state: Dict[str, Any], title: str = "Molecule Builder") -> str:
    """
    Generate MDL Molfile (V2000) format.
    Explicitly encodes bond orders (1, 2, 3) for Py3Dmol.
    """
    atoms = state["atoms"]
    bonds = state["bonds"]
    num_atoms = len(atoms)
    num_bonds = len(bonds)

    lines = [
        title[:80],
        "  MoleculeBuilder 3D",
        "",
        f"{num_atoms:>3}{num_bonds:>3}  0  0  0  0  0  0  0  0999 V2000"
    ]

    id_to_idx = {a["id"]: i + 1 for i, a in enumerate(atoms)}

    for a in atoms:
        elem = a["element"]
        lines.append(f"{a['x']:>10.4f}{a['y']:>10.4f}{a['z']:>10.4f} {elem:<3} 0  0  0  0  0  0  0  0  0  0  0  0")

    for b in bonds:
        i1 = id_to_idx.get(b["atom1"])
        i2 = id_to_idx.get(b["atom2"])
        if i1 and i2:
            order = b.get("order", 1)
            lines.append(f"{i1:>3}{i2:>3}{order:>3}  0  0  0  0")

    lines.append("M  END")
    return "\n".join(lines) + "\n"


def get_builder_mol_data(state: Dict[str, Any]) -> Dict[str, Any]:
    """Convert builder state into the standardized dict format expected by visualizer.py."""
    pdb_content = export_to_pdb(state)
    mol_info = compute_molecular_info(state)

    atoms_list = []
    for idx, a in enumerate(state["atoms"], start=1):
        atoms_list.append({
            "serial": idx,
            "id": a["id"],
            "name": f"{a['element']}{idx}",
            "element": a["element"],
            "x": a["x"],
            "y": a["y"],
            "z": a["z"],
            "resName": "UNL",
            "chain": "A",
            "resSeq": 1
        })

    return {
        "filename": "builder_molecule.pdb",
        "file_type": "PDB",
        "raw_content": pdb_content,
        "title": "Built Molecule",
        "num_atoms": mol_info["num_atoms"],
        "num_bonds": mol_info["num_bonds"],
        "bonds_source": "Builder Connections",
        "element_counts": mol_info["element_counts"],
        "molecular_weight": mol_info["molecular_weight"],
        "formula": mol_info["formula"],
        "atoms": atoms_list,
        "chains": ["A"],
        "residues": ["UNL"]
    }


# =========================================================================
# 3D STRUCTURE OPTIMIZATION (RDKit + Fallback)
# =========================================================================

def fallback_force_directed_3d(state: Dict[str, Any]) -> bool:
    """
    Pure Python force-directed 3D relaxation algorithm (VSEPR-inspired).
    Used as an ultra-robust fallback if RDKit is not available or fails.
    Simulates spring bond forces, steric non-bonded repulsion, and angular separation.
    """
    atoms = state["atoms"]
    bonds = state["bonds"]
    n = len(atoms)
    if n <= 1:
        if n == 1:
            atoms[0]["x"], atoms[0]["y"], atoms[0]["z"] = 0.0, 0.0, 0.0
        return True

    # Build adjacency and target bond lengths
    adj: Dict[int, List[Tuple[int, float]]] = {a["id"]: [] for a in atoms}
    atom_map = {a["id"]: a for a in atoms}

    for b in bonds:
        a1, a2 = b["atom1"], b["atom2"]
        if a1 in atom_map and a2 in atom_map:
            elem1 = atom_map[a1]["element"].upper()
            elem2 = atom_map[a2]["element"].upper()
            r1 = DEFAULT_BOND_LENGTHS.get(elem1, 0.77)
            r2 = DEFAULT_BOND_LENGTHS.get(elem2, 0.77)
            d0 = r1 + r2
            order = b.get("order", 1)
            if order == 2:
                d0 *= 0.88
            elif order == 3:
                d0 *= 0.78
            adj[a1].append((a2, d0))
            adj[a2].append((a1, d0))

    # Add small random jitter if all coordinates are zero or overlapping
    for a in atoms:
        if abs(a["x"]) < 1e-4 and abs(a["y"]) < 1e-4 and abs(a["z"]) < 1e-4:
            a["x"] += random.uniform(-0.3, 0.3)
            a["y"] += random.uniform(-0.3, 0.3)
            a["z"] += random.uniform(-0.3, 0.3)

    # 150 iterations of dampened Verlet/force relaxation
    step_size = 0.12
    damping = 0.85
    velocities = {a["id"]: [0.0, 0.0, 0.0] for a in atoms}

    for iteration in range(160):
        forces = {a["id"]: [0.0, 0.0, 0.0] for a in atoms}

        # 1. Bond Spring Forces (Hooke's Law)
        for a1_id, nbrs in adj.items():
            a1 = atom_map[a1_id]
            for a2_id, d0 in nbrs:
                if a1_id < a2_id:
                    a2 = atom_map[a2_id]
                    dx = a2["x"] - a1["x"]
                    dy = a2["y"] - a1["y"]
                    dz = a2["z"] - a1["z"]
                    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                    if dist < 1e-5:
                        dist = 1e-5
                        dx, dy, dz = random.uniform(0.01, 0.05), random.uniform(0.01, 0.05), 0.02
                    delta = dist - d0
                    # Spring constant k ~ 4.5
                    f_mag = 4.5 * delta
                    fx = f_mag * (dx / dist)
                    fy = f_mag * (dy / dist)
                    fz = f_mag * (dz / dist)
                    forces[a1_id][0] += fx
                    forces[a1_id][1] += fy
                    forces[a1_id][2] += fz
                    forces[a2_id][0] -= fx
                    forces[a2_id][1] -= fy
                    forces[a2_id][2] -= fz

        # 2. Non-bonded Steric Repulsion (Lennard-Jones / Coulomb-like)
        atom_ids = list(atom_map.keys())
        for i in range(len(atom_ids)):
            id1 = atom_ids[i]
            a1 = atom_map[id1]
            elem1 = a1["element"].upper()
            vdw1 = VDW_RADII.get(elem1, 1.6)
            for j in range(i + 1, len(atom_ids)):
                id2 = atom_ids[j]
                a2 = atom_map[id2]
                dx = a1["x"] - a2["x"]
                dy = a1["y"] - a2["y"]
                dz = a1["z"] - a2["z"]
                dist_sq = dx*dx + dy*dy + dz*dz
                elem2 = a2["element"].upper()
                vdw2 = VDW_RADII.get(elem2, 1.6)
                cutoff = (vdw1 + vdw2) * 1.2
                if dist_sq < cutoff * cutoff:
                    dist = math.sqrt(dist_sq)
                    if dist < 1e-4:
                        dist = 1e-4
                    # Repulsion pushes apart
                    rep_mag = 1.8 / (dist * dist + 0.1)
                    fx = rep_mag * (dx / dist)
                    fy = rep_mag * (dy / dist)
                    fz = rep_mag * (dz / dist)
                    forces[id1][0] += fx
                    forces[id1][1] += fy
                    forces[id1][2] += fz
                    forces[id2][0] -= fx
                    forces[id2][1] -= fy
                    forces[id2][2] -= fz

        # 3. VSEPR Angular Repulsion between shared neighbors
        for aid, nbrs in adj.items():
            if len(nbrs) >= 2:
                for k1 in range(len(nbrs)):
                    id_n1, _ = nbrs[k1]
                    n1 = atom_map[id_n1]
                    for k2 in range(k1 + 1, len(nbrs)):
                        id_n2, _ = nbrs[k2]
                        n2 = atom_map[id_n2]
                        dx = n1["x"] - n2["x"]
                        dy = n1["y"] - n2["y"]
                        dz = n1["z"] - n2["z"]
                        dist_sq = dx*dx + dy*dy + dz*dz
                        if dist_sq < 9.0:
                            dist = math.sqrt(dist_sq)
                            if dist < 1e-4:
                                dist = 1e-4
                            ang_f = 2.0 / (dist * dist + 0.2)
                            forces[id_n1][0] += ang_f * (dx / dist)
                            forces[id_n1][1] += ang_f * (dy / dist)
                            forces[id_n1][2] += ang_f * (dz / dist)
                            forces[id_n2][0] -= ang_f * (dx / dist)
                            forces[id_n2][1] -= ang_f * (dy / dist)
                            forces[id_n2][2] -= ang_f * (dz / dist)

        # Update velocities and positions
        for aid, a in atom_map.items():
            vx, vy, vz = velocities[aid]
            fx, fy, fz = forces[aid]
            # Cap force to prevent explosion
            f_norm = math.sqrt(fx*fx + fy*fy + fz*fz)
            if f_norm > 8.0:
                fx = fx * (8.0 / f_norm)
                fy = fy * (8.0 / f_norm)
                fz = fz * (8.0 / f_norm)
            vx = (vx + fx * step_size) * damping
            vy = (vy + fy * step_size) * damping
            vz = (vz + fz * step_size) * damping
            velocities[aid] = [vx, vy, vz]
            a["x"] += vx * step_size
            a["y"] += vy * step_size
            a["z"] += vz * step_size

    # Center molecule around origin
    mean_x = sum(a["x"] for a in atoms) / n
    mean_y = sum(a["y"] for a in atoms) / n
    mean_z = sum(a["z"] for a in atoms) / n
    for a in atoms:
        a["x"] = round(a["x"] - mean_x, 4)
        a["y"] = round(a["y"] - mean_y, 4)
        a["z"] = round(a["z"] - mean_z, 4)

    return True


def optimize_3d_geometry(state: Dict[str, Any], add_hydrogens: bool = False) -> Tuple[bool, str]:
    """
    Generate clean 3D molecular geometry and perform energy minimization.
    Tries RDKit ETKDG + UFF first. If RDKit is unavailable or encounters
    unusual valences/errors, automatically employs the robust force-directed fallback.
    """
    if not state["atoms"]:
        return False, "No atoms present to optimize."

    push_history(state)

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        rw = Chem.RWMol()
        id_to_rd_idx: Dict[int, int] = {}

        for a in state["atoms"]:
            elem = a["element"]
            rd_atom = Chem.Atom(elem)
            idx = rw.AddAtom(rd_atom)
            id_to_rd_idx[a["id"]] = idx

        for b in state["bonds"]:
            idx1 = id_to_rd_idx.get(b["atom1"])
            idx2 = id_to_rd_idx.get(b["atom2"])
            if idx1 is not None and idx2 is not None:
                order = b.get("order", 1)
                btype = Chem.BondType.SINGLE
                if order == 2:
                    btype = Chem.BondType.DOUBLE
                elif order == 3:
                    btype = Chem.BondType.TRIPLE
                rw.AddBond(idx1, idx2, btype)

        mol = rw.GetMol()

        # Sanitize gracefully
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            try:
                Chem.SanitizeMol(
                    mol,
                    sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                )
            except Exception:
                pass

        # If user requested adding hydrogens
        if add_hydrogens:
            try:
                mol = Chem.AddHs(mol)
            except Exception:
                pass

        # 3D Conformer Generation
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        res = AllChem.EmbedMolecule(mol, params)
        if res < 0:
            res = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)

        if res >= 0:
            # UFF force field optimization
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                pass

            conf = mol.GetConformer()

            # Update coordinates
            if add_hydrogens and mol.GetNumAtoms() > len(state["atoms"]):
                # Synchronize newly added hydrogens back into state
                state["atoms"].clear()
                state["bonds"].clear()
                state["next_atom_id"] = 1
                state["next_bond_id"] = 1

                for i in range(mol.GetNumAtoms()):
                    rd_a = mol.GetAtomWithIdx(i)
                    pos = conf.GetAtomPosition(i)
                    state["atoms"].append({
                        "id": i + 1,
                        "element": rd_a.GetSymbol(),
                        "x": round(pos.x, 4),
                        "y": round(pos.y, 4),
                        "z": round(pos.z, 4)
                    })

                for bond in mol.GetBonds():
                    i1 = bond.GetBeginAtomIdx() + 1
                    i2 = bond.GetEndAtomIdx() + 1
                    rd_order = bond.GetBondTypeAsDouble()
                    ord_int = int(rd_order) if rd_order in (1.0, 2.0, 3.0) else 1
                    state["bonds"].append({
                        "id": len(state["bonds"]) + 1,
                        "atom1": i1,
                        "atom2": i2,
                        "order": ord_int
                    })

                state["next_atom_id"] = len(state["atoms"]) + 1
                state["next_bond_id"] = len(state["bonds"]) + 1
                state["selected_atom_id"] = 1
                return True, f"✨ 3D Optimization complete (RDKit UFF with {mol.GetNumAtoms()} atoms including Hydrogens)."
            else:
                for a in state["atoms"]:
                    idx = id_to_rd_idx.get(a["id"])
                    if idx is not None and idx < mol.GetNumAtoms():
                        pos = conf.GetAtomPosition(idx)
                        a["x"] = round(pos.x, 4)
                        a["y"] = round(pos.y, 4)
                        a["z"] = round(pos.z, 4)

                return True, "✨ 3D Optimization complete (RDKit ETKDG + UFF Energy Minimization)."

    except Exception as e:
        # Fallback to pure Python force-directed geometry
        pass

    fallback_force_directed_3d(state)
    return True, "✨ 3D Geometry optimized using built-in force-directed VSEPR relaxation."


def auto_add_hydrogens(state: Dict[str, Any]) -> Tuple[int, str]:
    """
    Saturate open valences by attaching Hydrogen atoms to satisfy standard chemistry.
    """
    if not state["atoms"]:
        return 0, "No atoms present."

    added_count = 0
    # Collect additions first so we don't mutate state during iteration
    additions = []
    for atom in state["atoms"]:
        aid = atom["id"]
        elem = atom["element"].upper()
        current = get_atom_bond_order_sum(state, aid)
        allowed = ELEMENT_VALENCES.get(elem, [1])
        target = allowed[0]
        needed = target - current
        if needed > 0 and elem != "H":
            additions.append((aid, needed))

    if not additions:
        return 0, "All atoms already satisfy their typical valences."

    push_history(state)
    for parent_id, needed in additions:
        for _ in range(needed):
            add_atom(state, element="H", attach_to_id=parent_id, bond_order=1)
            added_count += 1

    # Optimize after adding H
    optimize_3d_geometry(state)
    return added_count, f"✅ Added {added_count} Hydrogen atoms to complete valences."


# =========================================================================
# 3D PY3DMOL BUILDER VIEW GENERATOR
# =========================================================================

def create_builder_viewer(
    state: Dict[str, Any],
    selected_atom_id: Optional[int] = None,
    style: str = "Ball-and-stick",
    bg_color: str = "#0f172a",
    viewer_height: int = 560
) -> str:
    """
    Generate an interactive Py3Dmol WebGL view tailored for the builder workspace.
    Highlights the currently selected atom with a distinct glowing halo and label.
    Supports mouse orbit, zoom, pan, and real-time HUD controls.
    """
    import py3Dmol

    mol_content = export_to_mol_block(state)
    view = py3Dmol.view(width="100%", height=f"{viewer_height}px")
    view.addModel(mol_content, "mol")

    # Set background color
    view.setBackgroundColor(bg_color)

    # Ball-and-stick rendering with multiple bonds support
    stick_spec = {"radius": 0.15, "multipleBonds": True}
    sphere_spec = {"scale": 0.28}
    view.setStyle({}, {"stick": stick_spec, "sphere": sphere_spec})

    # Add persistent atom serial / element labels
    view.addPropertyLabels("elem", {}, {
        "fontColor": "#ffffff",
        "fontSize": 11,
        "showBackground": True,
        "backgroundColor": "#1e293b",
        "backgroundOpacity": 0.65,
        "alignment": "center"
    })

    # Highlight selected atom if one is selected
    if selected_atom_id is not None:
        sel_atom = get_atom_by_id(state, selected_atom_id)
        if sel_atom:
            # 1-indexed index in molfile
            idx_in_mol = None
            for idx, a in enumerate(state["atoms"]):
                if a["id"] == selected_atom_id:
                    idx_in_mol = idx
                    break

            # Add glowing selection sphere around selected atom
            view.addSphere({
                "center": {"x": sel_atom["x"], "y": sel_atom["y"], "z": sel_atom["z"]},
                "radius": 0.58,
                "color": "#facc15",
                "opacity": 0.55
            })

            # Add prominent highlighted label
            view.addLabel(
                f"★ Selected: {sel_atom['element']} #{selected_atom_id}",
                {
                    "position": {"x": sel_atom["x"], "y": sel_atom["y"] + 0.65, "z": sel_atom["z"]},
                    "backgroundColor": "#facc15",
                    "fontColor": "#0f172a",
                    "fontSize": 12,
                    "borderThickness": 1,
                    "borderColor": "#ca8a04",
                    "inFront": True
                }
            )

    # Hover tooltips
    hover_js = """function(atom, viewer) {
        if (!atom.label) {
            var txt = (atom.elem || '') + (atom.serial ? (' #' + atom.serial) : '');
            atom.label = viewer.addLabel(txt, {
                position: atom,
                backgroundColor: '#0f172a',
                backgroundOpacity: 0.9,
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

    view.zoomTo()

    base_html = view._make_html()
    uniqueid = view.uniqueid

    # HUD Toolbar and overlay styles
    is_dark = bg_color.lower() in ("#0f172a", "#18181b", "#0a0f1d", "#000000", "#111827", "#1e1e2f")
    hud_bg = "rgba(15, 23, 42, 0.85)" if is_dark else "rgba(255, 255, 255, 0.9)"
    hud_border = "rgba(56, 189, 248, 0.35)" if is_dark else "rgba(203, 213, 225, 0.85)"
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
        .builder-hud {{
            position: absolute;
            top: 12px;
            right: 12px;
            z-index: 1000;
            display: flex;
            align-items: center;
            gap: 6px;
            background: {hud_bg};
            backdrop-filter: blur(10px);
            -webkit-backdrop-filter: blur(10px);
            padding: 5px 8px;
            border-radius: 10px;
            border: 1px solid {hud_border};
            box-shadow: 0 4px 18px rgba(0, 0, 0, 0.3);
        }}
        .builder-btn {{
            background: transparent;
            border: none;
            color: {btn_color};
            cursor: pointer;
            padding: 5px 9px;
            border-radius: 6px;
            font-size: 12px;
            font-weight: 500;
            display: flex;
            align-items: center;
            gap: 4px;
            transition: all 0.15s ease;
        }}
        .builder-btn:hover {{
            background: {btn_hover};
            color: #38bdf8;
            transform: translateY(-1px);
        }}
        .builder-btn:active {{
            transform: translateY(0);
        }}
        .builder-badge {{
            position: absolute;
            bottom: 12px;
            left: 12px;
            z-index: 1000;
            background: {hud_bg};
            backdrop-filter: blur(10px);
            border: 1px solid {hud_border};
            border-radius: 8px;
            padding: 4px 10px;
            font-size: 11px;
            color: {btn_color};
            opacity: 0.9;
            display: flex;
            gap: 12px;
            pointer-events: none;
        }}
    </style>

    <div class="builder-hud">
        <button class="builder-btn" id="btn-reset" title="Reset View Angle & Zoom" onclick="resetBuilderView()">
            <span>🔄</span> Reset
        </button>
        <button class="builder-btn" id="btn-spin" title="Toggle Auto-Spin Animation" onclick="toggleBuilderSpin()">
            <span>✨</span> Spin
        </button>
        <button class="builder-btn" id="btn-zoomin" title="Zoom In" onclick="zoomInBuilder()">
            <span>➕</span>
        </button>
        <button class="builder-btn" id="btn-zoomout" title="Zoom Out" onclick="zoomOutBuilder()">
            <span>➖</span>
        </button>
        <button class="builder-btn" id="btn-screenshot" title="Download Screenshot PNG" onclick="exportBuilderPNG()">
            <span>📸</span> PNG
        </button>
    </div>

    <div class="builder-badge">
        <span>🖱️ <b>Drag:</b> Rotate</span>
        <span>🔍 <b>Scroll:</b> Zoom</span>
        <span>✋ <b>Right-Drag:</b> Pan</span>
    </div>

    <script>
    var isSpinning = false;
    function resetBuilderView() {{
        if (typeof viewer_{uniqueid} !== 'undefined') {{
            viewer_{uniqueid}.zoomTo(500);
            viewer_{uniqueid}.render();
        }}
    }}
    function toggleBuilderSpin() {{
        if (typeof viewer_{uniqueid} !== 'undefined') {{
            isSpinning = !isSpinning;
            viewer_{uniqueid}.spin(isSpinning, 1);
        }}
    }}
    function zoomInBuilder() {{
        if (typeof viewer_{uniqueid} !== 'undefined') {{
            viewer_{uniqueid}.zoom(1.25);
            viewer_{uniqueid}.render();
        }}
    }}
    function zoomOutBuilder() {{
        if (typeof viewer_{uniqueid} !== 'undefined') {{
            viewer_{uniqueid}.zoom(0.8);
            viewer_{uniqueid}.render();
        }}
    }}
    function exportBuilderPNG() {{
        if (typeof viewer_{uniqueid} !== 'undefined') {{
            var uri = viewer_{uniqueid}.pngURI();
            var link = document.createElement('a');
            link.download = 'molecule_builder.png';
            link.href = uri;
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}
    }}
    </script>
    """

    return f"""<!DOCTYPE html>
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
