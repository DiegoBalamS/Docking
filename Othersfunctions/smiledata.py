pip install pyth3

from pyth.plugins.rtf15.reader import Rtf15Reader
from pyth.plugins.plaintext.writer import PlaintextWriter
import re

# Definir la ruta del archivo subido
file_path = 'ncidb.rtf'

# Leer el contenido del archivo RTF
with open(file_path, 'rb') as file:
    doc = Rtf15Reader.read(file)

# Convertir el contenido a texto plano
text_content = PlaintextWriter.write(doc).getvalue()

# Buscar las líneas que contienen '><SMILES>'
smiles_lines = re.findall(r'> <SMILES>\n\n(.+?)\n\n', text_content, re.DOTALL)

output_file_path = 'smiles_list.txt'

# Escribir los SMILES en un archivo
with open(output_file_path, 'w') as file:
    for smile in smiles_lines:
        file.write(smile + '\n')
