"""
Molecular file parser for PDB and XYZ formats.
Extracts atoms, coordinates, bonds, elemental composition, and structural metadata.
"""

from __future__ import annotations
import math
from typing import Dict, List, Tuple, Any, Optional

try:
    import pandas as pd
except Exception:
    pd = None


# Standard covalent radii (in Angstroms) for bond estimation
COVALENT_RADII: Dict[str, float] = {
    'H': 0.31, 'HE': 0.28,
    'LI': 1.28, 'BE': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'NE': 0.58,
    'NA': 1.66, 'MG': 1.41, 'AL': 1.21, 'SI': 1.11, 'P': 1.07, 'S': 1.05, 'CL': 1.02, 'AR': 1.06,
    'K': 2.03, 'CA': 1.76, 'SC': 1.70, 'TI': 1.60, 'V': 1.53, 'CR': 1.39, 'MN': 1.39,
    'FE': 1.32, 'CO': 1.26, 'NI': 1.24, 'CU': 1.32, 'ZN': 1.22, 'GA': 1.22, 'GE': 1.20,
    'AS': 1.19, 'SE': 1.20, 'BR': 1.20, 'KR': 1.16,
    'RB': 2.20, 'SR': 1.95, 'Y': 1.90, 'ZR': 1.75, 'MO': 1.54, 'RU': 1.46, 'RH': 1.42,
    'PD': 1.39, 'AG': 1.45, 'CD': 1.44, 'IN': 1.42, 'SN': 1.39, 'SB': 1.39, 'TE': 1.38,
    'I': 1.39, 'XE': 1.40,
    'CS': 2.44, 'BA': 2.15, 'PT': 1.36, 'AU': 1.36, 'HG': 1.32, 'PB': 1.46, 'U': 1.96
}

# Standard atomic weights (g/mol)
ATOMIC_WEIGHTS: Dict[str, float] = {
    'H': 1.008, 'HE': 4.0026,
    'LI': 6.94, 'BE': 9.0122, 'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'NE': 20.180,
    'NA': 22.990, 'MG': 24.305, 'AL': 26.982, 'SI': 28.085, 'P': 30.974, 'S': 32.06, 'CL': 35.45, 'AR': 39.948,
    'K': 39.098, 'CA': 40.078, 'SC': 44.956, 'TI': 47.867, 'V': 50.942, 'CR': 51.996, 'MN': 54.938,
    'FE': 55.845, 'CO': 58.933, 'NI': 58.693, 'CU': 63.546, 'ZN': 65.38, 'GA': 69.723, 'GE': 72.63,
    'AS': 74.922, 'SE': 78.971, 'BR': 79.904, 'KR': 83.798,
    'RB': 85.468, 'SR': 87.62, 'Y': 88.906, 'ZR': 91.224, 'MO': 95.95, 'RU': 101.07, 'RH': 102.91,
    'PD': 106.42, 'AG': 107.87, 'CD': 112.41, 'IN': 114.82, 'SN': 118.71, 'SB': 121.76, 'TE': 127.60,
    'I': 126.90, 'XE': 131.29,
    'CS': 132.91, 'BA': 137.33, 'PT': 195.08, 'AU': 196.97, 'HG': 200.59, 'PB': 207.2, 'U': 238.03
}

# Standard CPK element colors for UI chips
CPK_COLORS: Dict[str, str] = {
    'H': '#FFFFFF',
    'C': '#909090',
    'N': '#3050F8',
    'O': '#FF0D0D',
    'F': '#90E050',
    'CL': '#1FF01F',
    'BR': '#A62929',
    'I': '#940094',
    'P': '#FF8000',
    'S': '#FFFF30',
    'FE': '#E06633',
    'ZN': '#7D80B0',
    'CA': '#3DFF00',
    'MG': '#2A8000',
    'NA': '#AB5CF2',
    'K': '#8F40D4'
}

def clean_element_symbol(symbol: str) -> str:
    """Normalize and format element symbol (e.g., 'ca' -> 'Ca', 'C' -> 'C')."""
    if not symbol:
        return 'Unknown'
    symbol = symbol.strip().upper()
    # Remove digits or charges (+, -)
    symbol = ''.join([c for c in symbol if c.isalpha()])
    if len(symbol) == 1:
        return symbol
    elif len(symbol) >= 2:
        return symbol[0] + symbol[1:].lower()
    return symbol

def estimate_bonds(atoms: List[Dict[str, Any]], tolerance: float = 0.45) -> int:
    """
    Estimate number of covalent bonds between atoms based on interatomic distances
    and covalent radii heuristics (r1 + r2 + tolerance).
    """
    n = len(atoms)
    if n < 2:
        return 0
    
    # Cap bond computation for massive structures to preserve performance
    if n > 3500:
        return -1
    
    bonds = 0
    coords = [(a['x'], a['y'], a['z']) for a in atoms]
    symbols = [a['element'].upper() for a in atoms]
    
    for i in range(n):
        x1, y1, z1 = coords[i]
        r1 = COVALENT_RADII.get(symbols[i], 0.77)
        for j in range(i + 1, n):
            x2, y2, z2 = coords[j]
            dx = x1 - x2
            dy = y1 - y2
            dz = z1 - z2
            dist_sq = dx*dx + dy*dy + dz*dz
            
            # Fast bounding check: max bond length ~ 3.5 A (dist_sq <= 12.25)
            if dist_sq > 14.0 or dist_sq < 0.16:
                continue
            
            dist = math.sqrt(dist_sq)
            r2 = COVALENT_RADII.get(symbols[j], 0.77)
            if dist <= (r1 + r2 + tolerance):
                bonds += 1
                
    return bonds

def parse_xyz(content: str, filename: str = "molecule.xyz") -> Dict[str, Any]:
    """
    Parse an XYZ molecular file.
    Standard format:
    Line 1: Number of atoms (N)
    Line 2: Title / Comment
    Lines 3..N+2: Element X Y Z
    """
    lines = [line.strip() for line in content.strip().splitlines() if line.strip()]
    if not lines:
        raise ValueError("XYZ file is empty.")
    
    # Check first line: atom count
    try:
        atom_count_decl = int(lines[0].split()[0])
    except Exception:
        raise ValueError("First line of XYZ must specify the integer number of atoms.")
    
    title = lines[1] if len(lines) > 1 else "Unnamed Molecule"
    data_lines = lines[2:]
    
    atoms = []
    element_counts: Dict[str, int] = {}
    
    for idx, line in enumerate(data_lines, start=1):
        parts = line.split()
        if len(parts) < 4:
            continue
        raw_elem = parts[0]
        elem_norm = clean_element_symbol(raw_elem)
        
        try:
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
        except ValueError:
            raise ValueError(f"Invalid coordinate format at line {idx + 2}: '{line}'")
        
        atoms.append({
            'serial': idx,
            'name': f"{elem_norm}{idx}",
            'element': elem_norm,
            'x': x,
            'y': y,
            'z': z,
            'resName': 'MOL',
            'chain': 'A',
            'resSeq': 1
        })
        
        element_counts[elem_norm] = element_counts.get(elem_norm, 0) + 1
        
        if len(atoms) >= atom_count_decl:
            break

    if not atoms:
        raise ValueError("No valid atom coordinate records found in XYZ file.")
    
    # Estimate bonds using distance matrix
    computed_bonds = estimate_bonds(atoms)
    
    # Molecular weight calculation
    mw = sum(ATOMIC_WEIGHTS.get(elem.upper(), 12.0) * count for elem, count in element_counts.items())
    
    # Empirical chemical formula (Hill system: C, H, then alphabetical)
    sorted_elements = sorted(element_counts.keys(), key=lambda e: (0 if e == 'C' else (1 if e == 'H' else 2), e))
    formula_parts = []
    for elem in sorted_elements:
        c = element_counts[elem]
        formula_parts.append(f"{elem}{c if c > 1 else ''}")
    formula = "".join(formula_parts)
    
    df_atoms = pd.DataFrame(atoms) if pd is not None else atoms
    
    return {
        'filename': filename,
        'file_type': 'XYZ',
        'raw_content': content,
        'title': title,
        'num_atoms': len(atoms),
        'declared_atoms': atom_count_decl,
        'num_bonds': computed_bonds,
        'bonds_source': 'Estimated (Covalent radii distances)',
        'element_counts': element_counts,
        'molecular_weight': round(mw, 3),
        'formula': formula,
        'atoms': atoms,
        'dataframe': df_atoms,
        'has_secondary_structure': False,
        'chains': ['A'],
        'residues': ['MOL']
    }

def parse_pdb(content: str, filename: str = "molecule.pdb") -> Dict[str, Any]:
    """
    Parse a PDB molecular file.
    Extracts ATOM and HETATM records, CONECT records, secondary structure info,
    and metadata.
    """
    lines = content.splitlines()
    if not lines:
        raise ValueError("PDB file is empty.")
    
    atoms = []
    element_counts: Dict[str, int] = {}
    chains = set()
    residues = set()
    conect_bonds_set = set()
    
    title = ""
    helix_count = 0
    sheet_count = 0
    
    for line in lines:
        rec_type = line[:6].strip()
        
        if rec_type in ('TITLE', 'HEADER', 'COMPND') and not title:
            title = line[10:70].strip()
            
        elif rec_type == 'HELIX':
            helix_count += 1
            
        elif rec_type == 'SHEET':
            sheet_count += 1
            
        elif rec_type in ('ATOM', 'HETATM'):
            try:
                serial = int(line[6:11].strip() or len(atoms) + 1)
            except ValueError:
                serial = len(atoms) + 1
                
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip() or 'UNK'
            chain_id = line[21:22].strip() or 'A'
            
            try:
                res_seq = int(line[22:26].strip() or 1)
            except ValueError:
                res_seq = 1
                
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                continue
            
            # Element symbol from cols 76-78 or derived from atom name
            elem_raw = line[76:78].strip()
            if not elem_raw:
                elem_raw = ''.join([c for c in atom_name if c.isalpha()])
                if len(elem_raw) > 2:
                    elem_raw = elem_raw[0]
            
            elem_norm = clean_element_symbol(elem_raw)
            if not elem_norm or elem_norm == 'Unknown':
                elem_norm = 'C'
            
            chains.add(chain_id)
            residues.add(f"{res_name}_{chain_id}_{res_seq}")
            element_counts[elem_norm] = element_counts.get(elem_norm, 0) + 1
            
            atoms.append({
                'serial': serial,
                'record': rec_type,
                'name': atom_name,
                'resName': res_name,
                'chain': chain_id,
                'resSeq': res_seq,
                'x': x,
                'y': y,
                'z': z,
                'element': elem_norm
            })
            
        elif rec_type == 'CONECT':
            # CONECT atom_serial bonded1 bonded2 ...
            tokens = line.split()
            if len(tokens) >= 3:
                try:
                    src = int(tokens[1])
                    for tgt_str in tokens[2:]:
                        tgt = int(tgt_str)
                        if src != tgt:
                            pair = tuple(sorted((src, tgt)))
                            conect_bonds_set.add(pair)
                except ValueError:
                    pass
    
    if not atoms:
        raise ValueError("No valid ATOM or HETATM records found in PDB file.")
    
    explicit_bonds = len(conect_bonds_set)
    if explicit_bonds > 0:
        num_bonds = explicit_bonds
        bonds_source = 'CONECT records'
    else:
        # If no explicit CONECT records, estimate from coordinates
        est = estimate_bonds(atoms)
        num_bonds = est if est >= 0 else 'Computed dynamically'
        bonds_source = 'Estimated from coordinates' if est >= 0 else 'Dynamic'
        
    mw = sum(ATOMIC_WEIGHTS.get(elem.upper(), 12.0) * count for elem, count in element_counts.items())
    
    # Formula
    sorted_elements = sorted(element_counts.keys(), key=lambda e: (0 if e == 'C' else (1 if e == 'H' else 2), e))
    formula_parts = []
    for elem in sorted_elements:
        c = element_counts[elem]
        formula_parts.append(f"{elem}{c if c > 1 else ''}")
    formula = "".join(formula_parts)
    
    has_secondary = (helix_count > 0 or sheet_count > 0 or len(residues) > 5)
    df_atoms = pd.DataFrame(atoms) if pd is not None else atoms
    
    return {
        'filename': filename,
        'file_type': 'PDB',
        'raw_content': content,
        'title': title if title else "PDB Molecular Structure",
        'num_atoms': len(atoms),
        'num_bonds': num_bonds,
        'bonds_source': bonds_source,
        'element_counts': element_counts,
        'molecular_weight': round(mw, 3),
        'formula': formula,
        'atoms': atoms,
        'dataframe': df_atoms,
        'has_secondary_structure': has_secondary,
        'helix_count': helix_count,
        'sheet_count': sheet_count,
        'chains': sorted(list(chains)),
        'num_residues': len(residues)
    }

def parse_molecular_file(filename: str, content: str) -> Dict[str, Any]:
    """
    Dispatch parser based on file extension (.pdb or .xyz).
    Handles uppercase/lowercase extensions and auto-detects if missing.
    """
    fn_lower = filename.lower()
    if fn_lower.endswith('.pdb'):
        return parse_pdb(content, filename)
    elif fn_lower.endswith('.xyz'):
        return parse_xyz(content, filename)
    else:
        # Attempt heuristic detection
        stripped = content.strip()
        lines = stripped.splitlines()
        if len(lines) > 0 and lines[0].strip().isdigit():
            return parse_xyz(content, filename)
        elif any(line.startswith(('ATOM  ', 'HETATM', 'HEADER', 'COMPND')) for line in lines[:30]):
            return parse_pdb(content, filename)
        else:
            raise ValueError(
                f"Unsupported file format '{filename}'. Only .pdb and .xyz files are supported."
            )
