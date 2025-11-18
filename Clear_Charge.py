#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
from tqdm import tqdm

# Disable RDKit warning logs
RDLogger.DisableLog('rdApp.*')

def standardize_smi(smiles, clearCharge=True, clearFrag=False, canonTautomer=False, isomeric=False):
    """Standardize SMILES string (generate Kekulé structure)"""
    try:
        # Skip empty values
        if pd.isna(smiles) or smiles in ['', 'nan', 'NaN']:
            return None
            
        clean_mol = Chem.MolFromSmiles(smiles)
        if clean_mol is None:
            return None
        
        # Mandatory steps (not controlled by parameters)
        clean_mol = rdMolStandardize.Cleanup(clean_mol)  # Basic cleanup
        
        # Parameter-controlled steps
        if clearCharge:
            clean_mol = rdMolStandardize.Uncharger().uncharge(clean_mol)
        if clearFrag:
            clean_mol = rdMolStandardize.FragmentParent(clean_mol)
        if canonTautomer:
            clean_mol = rdMolStandardize.TautomerEnumerator().Canonicalize(clean_mol)
        
        # Generate SMILES with Kekulé structure
        Chem.Kekulize(clean_mol)
        return Chem.MolToSmiles(clean_mol, kekuleSmiles=True, isomericSmiles=isomeric)
    except Exception as e:
        print(f"Error processing '{smiles}': {str(e)}", file=sys.stderr)
        return None

def process_csv(input_path, output_path, smiles_col, **kwargs):
    """Process a single CSV file"""
    df = pd.read_csv(input_path)
    if smiles_col not in df.columns:
        raise ValueError(f"Column '{smiles_col}' not found in {input_path}")
    
    # Add progress bar
    tqdm.pandas(desc=f"Processing {os.path.basename(input_path)}")
    df['standardized_smiles'] = df[smiles_col].progress_apply(
        lambda x: standardize_smi(x, **kwargs)
    )
    
    # Count invalid SMILES
    invalid_count = df['standardized_smiles'].isna().sum()
    if invalid_count > 0:
        print(f"Warning: {invalid_count} invalid SMILES (including empty values) found", file=sys.stderr)
    
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")

def main():
    parser = argparse.ArgumentParser(
        description="SMILES standardization tool (generates Kekulé structure)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, 
                       help="Input path (file or folder)")
    parser.add_argument("-o", "--output", required=True,
                       help="Output path (file or folder)")
    parser.add_argument("-s", "--smiles-col", default="SMILES",
                       help="Specify SMILES column name")
    parser.add_argument("--fragment", action="store_true",
                       help="Enable fragment cleanup")
    parser.add_argument("--tautomer", action="store_true",
                       help="Enable tautomer standardization")
    parser.add_argument("--isomeric", action="store_true",
                       help="Preserve stereochemistry information")
    args = parser.parse_args()

    kwargs = {
        'clearCharge': True,  # Enabled by default
        'clearFrag': args.fragment,
        'canonTautomer': args.tautomer,
        'isomeric': args.isomeric
    }

    # Process input and output
    if os.path.isdir(args.input):
        os.makedirs(args.output, exist_ok=True)
        for fname in tqdm(os.listdir(args.input), desc="Processing files"):
            if fname.endswith('.csv'):
                process_csv(
                    input_path=os.path.join(args.input, fname),
                    output_path=os.path.join(args.output, f"std_{fname}"),
                    smiles_col=args.smiles_col,
                    **kwargs
                )
    else:
        process_csv(args.input, args.output, args.smiles_col, **kwargs)

if __name__ == "__main__":
    main()
