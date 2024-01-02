from Bio.PDB import PDBParser

def amino_acids_in_grid(pdb_file, center, size):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    amino_acids = []

    x_min, x_max = center[0] - size[0]/2, center[0] + size[0]/2
    y_min, y_max = center[1] - size[1]/2, center[1] + size[1]/2
    z_min, z_max = center[2] - size[2]/2, center[2] + size[2]/2

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if x_min <= atom.coord[0] <= x_max and \
                       y_min <= atom.coord[1] <= y_max and \
                       z_min <= atom.coord[2] <= z_max:
                        amino_acids.append(residue)
                        break

    return amino_acids

# Define el centro y tamaño del grid (ejemplo)
center = (10.0, 10.0, 10.0)  # x, y, z
size = (20.0, 20.0, 20.0)  # ancho, largo, alto

# Archivo PDB de la proteína
pdb_file = 'ruta/al/archivo/proteina.pdb'

# Obtener los aminoácidos dentro del grid
amino_acids = amino_acids_in_grid(pdb_file, center, size)

# Imprimir los aminoácidos encontrados
for residue in amino_acids:
    print("Residuo:", residue.get_resname(), "Cadena:", residue.get_parent().id, "Número de residuo:", residue.get_id()[1])

