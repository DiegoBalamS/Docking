import os
from subprocess import call

# Directorio donde se encuentran tus archivos PDB
pdb_directory = "/ruta/a/tus/archivos/pdb"
# Directorio donde quieres guardar los archivos PDBQT
pdbqt_directory = "/ruta/a/tus/archivos/pdbqt"

# Preparar los archivos PDB para conversión
for filename in os.listdir(pdb_directory):
    if filename.endswith(".pdb"):
        pdb_path = os.path.join(pdb_directory, filename)
        pdbqt_path = os.path.join(pdbqt_directory, filename.replace(".pdb", ".pdbqt"))

        # Comando para preparar el archivo PDB
        prepare_ligand_command = f"prepare_ligand4.py -l {pdb_path} -o {pdbqt_path}"

        # Ejecutar el comando
        call(prepare_ligand_command, shell=True)

print("Conversión completada.")
