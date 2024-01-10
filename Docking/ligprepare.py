import logging
import os
import pathlib
import re
import subprocess
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
from subprocess import call
 

from rdkit.Chem import AllChem as Chem



from tools import (
    PathType,
    assign_bond_orders,
    assign_stereochemistry,
    canonicalize_smiles,
    check_charges,
    check_mol,
    check_obabel_install,
    check_vina_output,
    convert_mol_file_to_pdbqt,
    convert_pdbqt_to_pdb,
    embed_mol,
    get_resources_dir,
    get_targets_dir,
    get_vina_path,
    parse_affinities_from_output,
    parse_search_box_conf,
    protonate_mol,
    read_mol_from_pdb,
    refine_mol_with_ff,
    sanitize_mol,
    smiles_to_mol,
    verify_docked_ligand,
    write_mol_to_mol_file,
)



def prep(
        smiles: str,
        pH=7.4,
        num_cpus: Optional[int] = None,
        seed=974528263,
        verbose=False,
    ) -> Tuple[Optional[float], Dict[str, Any]]:
        """
        Given a molecule, this method will prepare the ligand for the docking against the target.

        :param smiles: SMILES string of ligand
        :param pH: pH at which the docking should take place (default: 7.4, don't change unless you know what you are doing)
        :param num_cpus: number of CPUs cores available to AutoDock Vina
        :param seed: random seed for conformation generation and docking
        :param verbose: increase verbosity of log messages
        :return: docking score and dictionary containing all poses and binding free energies
        """
        # Auxiliary files
        ligand_mol_file = 'ligand.mol'
        ligand_pdbqt = 'ligand.pdbqt'

        # Make sure user input is standardized
        canonical_smiles = canonicalize_smiles(smiles)

        # Read and check input
        mol = smiles_to_mol(canonical_smiles, verbose=verbose)
        mol = sanitize_mol(mol, verbose=verbose)
        check_mol(mol)
        check_charges(mol)

        # Check that the right Open Babel version is available
        check_obabel_install()

        # Protonate ligand
        protonated_mol = protonate_mol(mol, pH=pH)
        check_mol(protonated_mol)

        # Embed ligand
        embedded_mol = embed_mol(protonated_mol, seed=seed)
        refined_mol = refine_mol_with_ff(embedded_mol)
        assign_stereochemistry(refined_mol)

        # Dock
        write_mol_to_mol_file(refined_mol, ligand_mol_file)
        convert_mol_file_to_pdbqt(ligand_mol_file, ligand_pdbqt)

        return ligand_pdbqt


def pre_ligand_adt(input_ligand, output_ligand_pdbqt):
    # Ruta al script prepare_ligand4.py de AutoDockTools
    prepare_ligand_script = '/ruta/a/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'

    # Comando para preparar el ligando
    command = f'pythonsh {prepare_ligand_script} -l {input_ligand} -o {output_ligand_pdbqt} -A hydrogens'

    # Ejecutar el comando
    call(command, shell=True)

# Ruta al archivo del ligando en formato PDB
#input_ligand_file = 'ruta/al/archivo/ligando.pdb'

# Ruta al archivo de salida en formato PDBQT
#output_ligand_pdbqt_file = 'ruta/al/archivo/ligando.pdbqt'

def prepare_ligand(smiles, output_ligand_pdbqt):
    """
    Prepara un ligando para el docking realizando las siguientes operaciones:
    - Agrega hidrógenos.
    - Realiza merge de los hidrógenos no polares.
    - Agrega cargas de Gasteiger.
    - Convierte el archivo a formato PDBQT.
    
    Args:
    input_ligand (str): Ruta al archivo del ligando en formato PDB.
    output_ligand_pdbqt (str): Ruta al archivo de salida en formato PDBQT.
    """

# Auxiliary files
    ligand_mol_file = 'ligand.mol'
    ligand_pdbqt = 'ligand.pdbqt'
    vina_logfile = 'vina.log'
    vina_outfile = 'vina.out'
    docked_ligand_pdb ='docked_ligand.pdb'


    mol=Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    Chem.MolToMolFile(mol, ligand_mol_file)
    cmd_list = [
        'obabel',
        '-imol', ligand_mol_file,
        '-opdbqt',
        '-O', output_ligand_pdbqt,
        '-h','-p',
        '--partialcharge', 'gasteiger'
    ]
    # yapf: enable
    cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd_return.stdout.decode('utf-8')
    logging.debug(output)
    # Comando para preparar el ligando con obabel
    #command = f'obabel -ipdb {pdb_block} -opdbqt -O {output_ligand_pdbqt} -h -p --partialcharge gasteiger'

    # Ejecutar el comando
    #call(command, shell=True)
