import logging
import os
import pathlib
import re
import subprocess
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



def prep(
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

    
