import os
import numpy as np
from Bio.PDB import PDBParser

def calculate_grid_center_and_size(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    atoms = [atom.get_coord() for atom in structure.get_atoms() if atom.get_name() == 'CA']
    atoms = np.array(atoms)

    min_coord = np.min(atoms, axis=0)
    max_coord = np.max(atoms, axis=0)
    center = (max_coord + min_coord) / 2
    size = max_coord - min_coord

    return center, size

# Ruta al archivo PDB de la proteína
pdb_file = 'ruta/a/tu/archivo/proteina.pdb'

center, size = calculate_grid_center_and_size(pdb_file)

print("Centro de la cuadrícula:", center)
print("Tamaño de la cuadrícula:", size)
