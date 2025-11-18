import pandas as pd
import argparse

def load_smiles(file_path):
    """Load CSV file and return a set of Canonical_SMILES column"""
    try:
        df = pd.read_csv(file_path)
        if 'Canonical_SMILES' not in df.columns:
            raise ValueError(f"Error: {file_path} does not contain a 'Canonical_SMILES' column.")
        return set(df['Canonical_SMILES'].dropna().astype(str).str.strip())
    except Exception as e:
        print(f"Error loading file {file_path}: {e}")
        raise

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Compare two SMILES datasets and output unique SMILES")
    parser.add_argument('-i1', '--input1', required=True, help="First input file path (COCONUT)")
    parser.add_argument('-i2', '--input2', required=True, help="Second input file path (CMNPD)")
    parser.add_argument('-o', '--output', required=True, help="Output file path")
    args = parser.parse_args()

    try:
        # Load SMILES data
        print("Loading SMILES data...")
        coconut_smiles = load_smiles(args.input1)
        mnp_smiles = load_smiles(args.input2)

        # Calculate SMILES unique to the first file
        print("Comparing SMILES sets...")
        unique_smiles = coconut_smiles - mnp_smiles

        # Save results
        print(f"Saving results to {args.output}...")
        df_output = pd.DataFrame(list(unique_smiles), columns=['Canonical_SMILES'])
        df_output.to_csv(args.output, index=False)

        print(f"Success! Saved {len(unique_smiles)} unique SMILES to {args.output}")
    except Exception as e:
        print(f"Error occurred: {e}")

if __name__ == "__main__":
    main()
