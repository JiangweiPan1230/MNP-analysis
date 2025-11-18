#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import argparse
from tqdm import tqdm
import logging
from typing import Optional

# Set up logging
logging.basicConfig(
    filename='smiles_processing.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def is_valid_molecule(smiles: str) -> bool:
    """Check if SMILES represents a valid molecule (including chemical bond validation)
    
    Args:
        smiles: Input SMILES string
        
    Returns:
        bool: Whether it's a valid molecule
    """
    if not isinstance(smiles, str):
        return False
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    
    try:
        Chem.SanitizeMol(mol)
        return True
    except:
        return False

def get_canonical_kekulized_smiles(smiles: str) -> Optional[str]:
    """Convert SMILES string to canonical form with Kekulé structure
    
    Args:
        smiles: Input SMILES string
        
    Returns:
        Optional[str]: Canonical Kekulé SMILES, returns None if invalid
    """
    if not is_valid_molecule(smiles):
        return None
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    try:
        Chem.SanitizeMol(mol)
        Chem.Kekulize(mol)
        return Chem.MolToSmiles(mol, kekuleSmiles=True, canonical=True)
    except:
        return None

def canonical_smiles(smiles: str) -> Optional[str]:
    """Convert SMILES string to canonical form
    
    Args:
        smiles: Input SMILES string
        
    Returns:
        Optional[str]: Canonical SMILES, returns None if invalid
    """
    if not is_valid_molecule(smiles):
        return None
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except:
        return None

def process_file(file_path: str, output_folder: str, kekule: bool = False, smiles_column: str = 'SMILES') -> None:
    """Read and process a single CSV file
    
    Args:
        file_path: Input file path
        output_folder: Output folder path
        kekule: Whether to use Kekulé structure processing
        smiles_column: SMILES column name
    """
    try:
        df = pd.read_csv(file_path)

        if smiles_column not in df.columns:
            logging.error(f"File {file_path} does not contain '{smiles_column}' column")
            return

        original_count = len(df)
        process_func = get_canonical_kekulized_smiles if kekule else canonical_smiles
        df['Canonical_SMILES'] = df[smiles_column].apply(process_func)
        df = df[df['Canonical_SMILES'].notna()]
        
        processed_count = len(df)
        removed_count = original_count - processed_count
        
        os.makedirs(output_folder, exist_ok=True)
        output_file_path = os.path.join(output_folder, os.path.basename(file_path))
        df.to_csv(output_file_path, index=False)
        
        logging.info(
            f"File {file_path} processing completed: "
            f"Original {original_count} entries, "
            f"Retained {processed_count} entries, "
            f"Removed {removed_count} invalid molecules"
        )
    except Exception as e:
        logging.error(f"Error processing file {file_path}: {e}")

def process_folder(folder_path: str, output_folder: str, kekule: bool = False, smiles_column: str = 'SMILES') -> None:
    """Traverse folder and process all CSV files
    
    Args:
        folder_path: Input folder path
        output_folder: Output folder path
        kekule: Whether to use Kekulé structure processing
        smiles_column: SMILES column name
    """
    os.makedirs(output_folder, exist_ok=True)
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

    for file in tqdm(csv_files, desc="Processing CSV files"):
        file_path = os.path.join(folder_path, file)
        process_file(file_path, output_folder, kekule, smiles_column)

def main():
    parser = argparse.ArgumentParser(
        description="Process SMILES columns in CSV files and generate canonical SMILES"
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help="Input CSV file or folder path"
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help="Output folder path"
    )
    parser.add_argument(
        '-k', '--kekule',
        action='store_true',
        help="Whether to generate canonical SMILES with Kekulé structure (default: False)"
    )
    parser.add_argument(
        '-s', '--smiles_column',
        default='SMILES',
        help="Specify SMILES column name (default: 'SMILES')"
    )
    args = parser.parse_args()

    if os.path.isfile(args.input):
        process_file(args.input, args.output, args.kekule, args.smiles_column)
    elif os.path.isdir(args.input):
        process_folder(args.input, args.output, args.kekule, args.smiles_column)
    else:
        error_msg = f"Input path {args.input} is neither a file nor a folder"
        logging.error(error_msg)
        print(f"Error: {error_msg}")

if __name__ == "__main__":
    main()
