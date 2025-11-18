import os
import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from tqdm import tqdm

def process_smiles(smiles, kekulize=True):
    """Process SMILES string and return scaffold SMILES with optional Kekulization"""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        
        if kekulize:
            # Ensure Kekulé structure
            Chem.Kekulize(scaffold, clearAromaticFlags=True)
            # Generate Kekulé SMILES
            scaffold_smiles = Chem.MolToSmiles(scaffold, kekuleSmiles=True)
        else:
            scaffold_smiles = Chem.MolToSmiles(scaffold)
            
        return scaffold_smiles
    except Exception:
        return None

def process_csv_file(file_path, output_path=None, smiles_column='SMILES', kekulize=True):
    """
    Process a single CSV file to extract scaffolds
    
    Args:
        file_path (str): Path to input CSV file
        output_path (str): Path to output CSV file (optional)
        smiles_column (str): Column name containing SMILES strings
        kekulize (bool): Whether to Kekulize the scaffolds (default: True)
    
    Returns:
        str: Path to output file if successful, None otherwise
    """
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
        
    try:
        df = pd.read_csv(file_path, dtype=str)  # Read all data as strings
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return None
        
    if smiles_column not in df.columns:
        print(f"Column '{smiles_column}' not found in {file_path}. Available columns: {', '.join(df.columns)}")
        return None
        
    if output_path is None:
        base_name = os.path.splitext(file_path)[0]
        output_path = f"{base_name}_scaffolds.csv"
    
    # Process with progress bar
    tqdm.pandas(desc=f"Processing {os.path.basename(file_path)}")
    df['Scaffold'] = df[smiles_column].progress_apply(lambda x: process_smiles(x, kekulize))
    
    # Remove rows where scaffold extraction failed
    df = df[df['Scaffold'].notna()]
    
    if not df.empty:
        df.to_csv(output_path, index=False)
        return output_path
    return None

def process_folder(input_folder, output_folder=None, smiles_column='SMILES', kekulize=True):
    """
    Process all CSV files in a folder
    
    Args:
        input_folder (str): Path to input folder
        output_folder (str): Path to output folder (optional)
        smiles_column (str): Column name containing SMILES strings
        kekulize (bool): Whether to Kekulize the scaffolds
    
    Returns:
        list: List of successfully processed files
    """
    if not os.path.exists(input_folder):
        print(f"Input folder not found: {input_folder}")
        return []
        
    if output_folder is None:
        output_folder = os.path.join(input_folder, "scaffold_output")
    os.makedirs(output_folder, exist_ok=True)
    
    csv_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.csv')]
    processed_files = []
    
    for file in tqdm(csv_files, desc="Processing folder"):
        input_path = os.path.join(input_folder, file)
        output_path = os.path.join(output_folder, file.replace('.csv', '_scaffolds.csv'))
        
        result = process_csv_file(input_path, output_path, smiles_column, kekulize)
        if result:
            processed_files.append(result)
    
    return processed_files

def main():
    parser = argparse.ArgumentParser(
        description='Extract Murcko scaffolds from SMILES in CSV files with optional Kekulization'
    )
    parser.add_argument('-i', '--input', required=True, 
                       help='Input CSV file or folder containing CSV files')
    parser.add_argument('-o', '--output', 
                       help='Output file or folder (optional)')
    parser.add_argument('-s', '--smiles_column', default='SMILES',
                       help='Name of the column containing SMILES (default: SMILES)')
    parser.add_argument('--no-kekulize', action='store_false', dest='kekulize',
                       help='Disable Kekulization of scaffolds')
    
    args = parser.parse_args()
    
    if os.path.isfile(args.input):
        print(f"Processing single file: {args.input}")
        output_file = process_csv_file(
            args.input, 
            args.output, 
            args.smiles_column, 
            args.kekulize
        )
        if output_file:
            print(f"Results saved to: {output_file}")
        else:
            print("Processing failed")
    elif os.path.isdir(args.input):
        print(f"Processing folder: {args.input}")
        output_files = process_folder(
            args.input, 
            args.output, 
            args.smiles_column, 
            args.kekulize
        )
        print(f"\nProcessed {len(output_files)} files. Results saved to: {os.path.dirname(output_files[0]) if output_files else 'N/A'}")
    else:
        print(f"Error: Input path {args.input} is neither a file nor a directory")

if __name__ == '__main__':
    main()
