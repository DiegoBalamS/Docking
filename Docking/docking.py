import os
from subprocess import call


def dock(protein_file,ligand_file,center,size):
"""
# Rutas a los archivos de la proteína y el ligando
protein_file = '/ruta/a/tu/proteina.pdbqt'
ligand_file = '/ruta/a/tu/ligando.pdbqt'

# Rutas a los archivos de salida
output_file = '/ruta/a/tu/salida_docking.pdbqt'
"""
# Parámetros de la cuadrícula (ajustar según sea necesario)
center_x, center_y, center_z = center  # Centro de la cuadrícula
size_x, size_y, size_z = size          # Tamaño de la cuadrícula

# Comando para ejecutar AutoDock Vina
vina_command = [
    '/Users/diegobalam/Documents/Docking-main-2/Docking/resources/bin/vina_mac_catalina', 
    '--receptor', protein_file,
    '--ligand', ligand_file,
    '--center_x', str(center_x),
    '--center_y', str(center_y),
    '--center_z', str(center_z),
    '--size_x', str(size_x),
    '--size_y', str(size_y),
    '--size_z', str(size_z),
    '--out', output_file
]

# Ejecutar el comando
call(vina_command)

print("Docking completado.")

