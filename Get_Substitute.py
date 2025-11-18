import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem
import argparse
from tqdm import tqdm

def get_scaffold_and_substituents(mol):
    """Extract scaffold and substituents from a molecule"""
    if not mol:
        return None, []
    
    # Get Murcko scaffold
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if not scaffold:
        return None, []
    
    # Get substituents
    substituents = []
    try:
        side_chains = Chem.ReplaceCore(mol, scaffold)
        if side_chains:
            for frag in Chem.GetMolFrags(side_chains, asMols=True):
                clean_frag = AllChem.ReplaceSubstructs(
                    frag, 
                    Chem.MolFromSmiles('[*]'), 
                    Chem.MolFromSmiles('[H]'), 
                    replaceAll=True
                )[0]
                substituents.append(Chem.MolToSmiles(clean_frag))
    except:
        pass
    
    return Chem.MolToSmiles(scaffold), substituents

def process_smiles(smiles):
    """Process a single SMILES string"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, None, []
    
    scaffold, substituents = get_scaffold_and_substituents(mol)
    return smiles, scaffold, substituents

def process_file(input_path, output_path, smiles_column='SMILES'):
    """Process a single CSV file"""
    try:
        df = pd.read_csv(input_path)
    except Exception as e:
        print(f"Error: Unable to read file {input_path} - {str(e)}")
        return False
    
    if smiles_column not in df.columns:
        print(f"Error: SMILES column '{smiles_column}' not found in {input_path}")
        return False
    
    results = []
    for smi in tqdm(df[smiles_column], desc=f"Processing {os.path.basename(input_path)}"):
        result = process_smiles(smi)
        if result[0]:  # Skip invalid SMILES
            results.append(result)
    
    if not results:
        print(f"Warning: No valid data found in {input_path}")
        return False
    
    # Convert to DataFrame
    max_subst = max(len(r[2]) for r in results) if results else 0
    data = {
        'SMILES': [r[0] for r in results],
        'Scaffold': [r[1] for r in results]
    }
    
    # Dynamically add substituent columns
    for i in range(1, max_subst + 1):
        data[f'Substituent_{i}'] = [
            r[2][i-1] if i <= len(r[2]) else '' for r in results
        ]
    
    result_df = pd.DataFrame(data)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    result_df.to_csv(output_path, index=False)
    return True

def process_folder(input_dir, output_dir, smiles_column='SMILES'):
    """Process an entire folder"""
    if not os.path.exists(input_dir):
        print(f"Error: Input folder does not exist {input_dir}")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    processed = 0
    for file in os.listdir(input_dir):
        if file.endswith('.csv'):
            input_path = os.path.join(input_dir, file)
            output_path = os.path.join(output_dir, f"processed_{file}")
            if process_file(input_path, output_path, smiles_column):
                processed += 1
    
    print(f"\nProcessing completed. Total files processed: {processed}")

def main():
    parser = argparse.ArgumentParser(description='Molecular Scaffold and Substituent Analysis Tool')
    parser.add_argument('-i', '--input', required=True, 
                       help='Input CSV file or folder path')
    parser.add_argument('-o', '--output', required=True,
                       help='Output CSV file or folder path')
    parser.add_argument('-s', '--smiles_column', default='SMILES',
                       help='Name of the column containing SMILES (default: "SMILES")')
    
    args = parser.parse_args()
    
    # Check if input path exists
    if not os.path.exists(args.input):
        print(f"Error: Input path does not exist {args.input}")
        return
    
    # Determine processing mode
    if os.path.isfile(args.input):
        # File mode
        if os.path.isdir(args.output):
            # If output is a directory, use input filename
            output_path = os.path.join(args.output, f"processed_{os.path.basename(args.input)}")
        else:
            # If output is a file path, use directly
            output_path = args.output
            if not output_path.endswith('.csv'):
                output_path += '.csv'
        
        process_file(args.input, output_path, args.smiles_column)
    elif os.path.isdir(args.input):
        # Folder mode
        if os.path.isfile(args.output):
            print("Error: When input is a folder, output must be a folder")
            return
        
        process_folder(args.input, args.output, args.smiles_column)
    else:
        print(f"Error: Input path is neither a file nor a folder {args.input}")

if __name__ == '__main__':
    main()
