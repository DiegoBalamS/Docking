import logging
import os
import pathlib
import re
from subprocess import call
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from rdkit.Chem import AllChem as Chem

from .errors import DockingError, DockstringError, VinaError
from .tools import (
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


def load_target(name: str, *args, **kwargs) -> 'Target':
    """
    Load target with name <name>.

    :return: a target
    """
    return Target(name, *args, **kwargs)


def list_all_target_names(targets_dir: Optional[PathType] = None) -> List[str]:
    """
    List all available targets in <targets_dir>.

    :return: list of available names targets
    """
    if targets_dir is None:
        targets_dir = get_targets_dir()
    file_names = [f for f in os.listdir(targets_dir) if os.path.isfile(os.path.join(targets_dir, f))]

    target_re = re.compile(r'^(?P<name>\w+)_target\.pdbqt$')
    names = []
    for file_name in file_names:
        match = target_re.match(file_name)
        if match:
            names.append(match.group('name'))

    return names


# Función para preparar el receptor
def prepare_receptor_adt(input_pdb, output_pdbqt):
    # Comando para agregar hidrógenos polares y hacer merge de los no polares
    add_hydrogens_command = f'pythonsh /ruta/a/adt/bin/pythonsh /ruta/a/adt/AutoDockTools/Utilities24/prepare_receptor4.py -r {input_pdb} -o temp.pdb -A hydrogens -U nphs'

    # Ejecutar el comando para agregar hidrógenos
    call(add_hydrogens_command, shell=True)

    # Comando para agregar cargas de Kollman y convertir a PDBQT
    add_charges_command = f'pythonsh /ruta/a/adt/bin/pythonsh /ruta/a/adt/AutoDockTools/Utilities24/prepare_receptor4.py -r temp.pdb -o {output_pdbqt} -A checkhydrogens'

    # Ejecutar el comando para agregar cargas y convertir a PDBQT
    call(add_charges_command, shell=True)

    # Eliminar el archivo temporal
    os.remove('temp.pdb')

# Ruta al archivo PDB del receptor
input_pdb_file = 'ruta/al/archivo/receptor.pdb'

# Ruta al archivo PDBQT de salida
output_pdbqt_file = 'ruta/al/archivo/salida_receptor.pdbqt'

# Llamar a la función para preparar el receptor
prepare_receptor(input_pdb_file, output_pdbqt_file)

print("Preparación del receptor completada.")

import openbabel

def prepare_receptor_babel(input_pdb, output_pdbqt):
    # Crear un objeto de conversión
    conversion = openbabel.OBConversion()

    # Especificar los formatos de entrada y salida
    conversion.SetInAndOutFormats("pdb", "pdbqt")

    # Crear un objeto molécula
    mol = openbabel.OBMol()

    # Leer el archivo de entrada
    conversion.ReadFile(mol, input_pdb)

    # Agregar hidrógenos polares
    mol.AddHydrogens(True, False)

    # Crear un objeto que manejará la adición de cargas
    charge_model = openbabel.OBChargeModel.FindType("gasteiger")

    # Calcular y asignar las cargas
    if charge_model:
        charge_model.ComputeCharges(mol)

    # Escribir el archivo de salida
    conversion.WriteFile(mol, output_pdbqt)

    # Cerrar el archivo de salida
    conversion.CloseOutFile()

# Rutas a los archivos de entrada y salida
input_pdb_file = 'ruta/al/archivo/receptor.pdb'
output_pdbqt_file = 'ruta/al/archivo/salida_receptor.pdbqt'

# Preparar el receptor
prepare_receptor(input_pdb_file, output_pdbqt_file)

print("Preparación del receptor completada.")



class Target:
    def __init__(
        self,
        name: str,
        working_dir: Optional[PathType] = None,
        targets_dir: Optional[PathType] = None,
    ) -> None:
        """
        Target to dock against. Two files are required: <name>_target.pdbqt and <name>_conf.txt.
        <name>_target.pdbqt contains the the protein structure including partial charges.
        <name>_conf.txt contains the coordinates of the search box.

        :param name: target name (e.g,. ABL1)
        :param working_dir: directory for temporary and output files. If None, a temporary directory will be created.
        :param targets_dir: directory in which the required files can be found. If None, a default path will be chosen.
        """
        self.name = name

        # Directory where the ligand and output files will be saved
        self._custom_working_dir = working_dir
        self._tmp_dir_handle: Optional[tempfile.TemporaryDirectory] = None
        self.targets_dir: Path = pathlib.Path(targets_dir) if targets_dir else get_targets_dir()

        # Ensure input files exist
        if not all(p.exists() for p in [self.pdbqt_path, self.conf_path]):
            raise DockstringError(f"'{self.name}' is not a supported target")

    def __repr__(self):
        return f"Target(name='{self.name}', working_dir='{self.working_dir}', targets_dir='{self.targets_dir}')"

    @property
    def pdbqt_path(self) -> Path:
        """
        Path to PDBQT file
        """
        return self.targets_dir / (self.name + '_target.pdbqt')

    @property
    def conf_path(self) -> Path:
        """
        Path to configuration file
        """
        return self.targets_dir / (self.name + '_conf.txt')

    @property
    def working_dir(self) -> Path:
        """
        Path to working directory
        """
        if self._custom_working_dir:
            return Path(self._custom_working_dir).resolve()

        # If no custom working dir is set and the tmp working dir handle is not initialized, initialize it
        if not self._tmp_dir_handle:
            self._tmp_dir_handle = tempfile.TemporaryDirectory()

        return Path(self._tmp_dir_handle.name).resolve()

    def prept(pdbqt_file: PathType
        pH=7.4,
        num_cpus: Optional[int] = None,
        seed=974528263,
        verbose=False,
    ) -> Tuple[Optional[float], Dict[str, Any]]:
        """
        Given a molecule, this method will return a docking score against the current target.

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
        vina_logfile = 'vina.log'
        vina_outfile = 'vina.out'
        docked_ligand_pdb ='docked_ligand.pdb'

    
        return ligand_pdbqt


    def _dock_pdbqt(self, pdbqt_path, log_path, out_path, seed, num_cpus: Optional[int] = None) -> None:
        """
        Run AutoDock Vina.

        :param pdbqt_path: path to PDBQT file
        :param log_path: path to log file
        :param out_path: path to output file
        :param seed: random seed
        :param num_cpus: number of CPU cores available to AutoDock Vina
        """
        # yapf: disable
        cmd_list = [
            get_vina_path(),
            '--receptor', self.pdbqt_path,
            '--config', self.conf_path,
            '--ligand', pdbqt_path,
            '--log', log_path,
            '--out', out_path,
            '--seed', str(seed),
        ]
        # yapf: enable
        if num_cpus is not None:
            cmd_list += ['--cpu', str(num_cpus)]

        cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = cmd_return.stdout.decode('utf-8')
        logging.debug(output)

        # If failure, raise DockingError
        if cmd_return.returncode != 0:
            raise VinaError(f'Docking with Vina failed: {output}')

    def dock(
        self,
        smiles: str,
        pH=7.4,
        num_cpus: Optional[int] = None,
        seed=974528263,
        verbose=False,
    ) -> Tuple[Optional[float], Dict[str, Any]]:
        """
        Given a molecule, this method will return a docking score against the current target.

        :param smiles: SMILES string of ligand
        :param pH: pH at which the docking should take place (default: 7.4, don't change unless you know what you are doing)
        :param num_cpus: number of CPUs cores available to AutoDock Vina
        :param seed: random seed for conformation generation and docking
        :param verbose: increase verbosity of log messages
        :return: docking score and dictionary containing all poses and binding free energies
        """
        # Auxiliary files
        ligand_mol_file = self.working_dir / 'ligand.mol'
        ligand_pdbqt = self.working_dir / 'ligand.pdbqt'
        vina_logfile = self.working_dir / 'vina.log'
        vina_outfile = self.working_dir / 'vina.out'
        docked_ligand_pdb = self.working_dir / 'docked_ligand.pdb'

        
        self._dock_pdbqt(ligand_pdbqt, vina_logfile, vina_outfile, seed=seed, num_cpus=num_cpus)

        # Process docking output
        try:
            check_vina_output(vina_outfile)
        except DockingError:
            return None, {}

        convert_pdbqt_to_pdb(pdbqt_file=vina_outfile, pdb_file=docked_ligand_pdb, disable_bonding=True)
        raw_ligand = read_mol_from_pdb(docked_ligand_pdb)

        # Assign bond orders and stereochemistry
        refined_mol_no_hs = Chem.RemoveHs(refined_mol)  # remove Hs as they are not present in the PDBQT file
        ligand = assign_bond_orders(subject=raw_ligand, ref=refined_mol_no_hs)
        assign_stereochemistry(ligand)

        # Verify docked ligand
        verify_docked_ligand(ref=refined_mol_no_hs, subject=ligand)

        # Parse scores
        affinities = parse_affinities_from_output(docked_ligand_pdb)
        assert len(affinities) == ligand.GetNumConformers()
        score = affinities[0]

        return score, {
            'ligand': ligand,
            'affinities': affinities,
        }

    def view(
        self,
        mol: Union[Chem.Mol, List[Chem.Mol], None] = None,
        include_search_box=True,
    ) -> int:
        """
        Launch PyMol and view the receptor and the search box.

        :param mol: RDKit molecule or list of RDKit molecules containing a conformation
        :param include_search_box: view search box
        :return: return code of PyMol command
        """
        commands: List[Union[str, PathType]] = ['pymol', self.pdbqt_path]

        if include_search_box:
            pymol_script = get_resources_dir() / 'view_search_box.py'
            conf = parse_search_box_conf(self.conf_path)
            # yapf: disable
            commands += [
                pymol_script,
                '-d', 'view_search_box center_x={center_x}, center_y={center_y}, center_z={center_z}, '
                      'size_x={size_x}, size_y={size_y}, size_z={size_z}'.format(**conf)
            ]
            # yapf: enable

        if mol:
            if not isinstance(mol, Iterable):
                mol = [mol]

            tmp_dir_handle = tempfile.TemporaryDirectory()
            tmp_dir = Path(tmp_dir_handle.name).resolve()

            for index, pose in enumerate(mol):
                mol_file = tmp_dir / f'ligand_{index}.mol'
                write_mol_to_mol_file(pose, mol_file)
                commands += [str(mol_file)]

        return subprocess.run(commands).returncode
