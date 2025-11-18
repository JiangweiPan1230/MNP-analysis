import pandas as pd
import os
import argparse
from tqdm import tqdm
import re
from datetime import datetime

# Define list of metal element symbols
metal_elements = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Ti', 'V', 
                 'Cr', 'Mn', 'Fe', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 
                 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 
                 'In',  'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 
                 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 
                 'Ta', 'W', 'Re', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 
                 'Bi', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Pu', 
                 'Am', 'Cm', 'Bk', 'Es', 'Fm', 'Md', 'Lr']

# Create regex pattern to match metal elements
metal_pattern = re.compile(r'(' + '|'.join(metal_elements) + r')', re.IGNORECASE)

def process_file(input_path, output_dir, log_file):
    try:
        # Read CSV file
        df = pd.read_csv(input_path)
        filename = os.path.basename(input_path)
        
        # Check if Canonical_SMILES column exists
        if 'Canonical_SMILES' not in df.columns:
            with open(log_file, 'a') as f:
                f.write(f"{datetime.now()}: {filename} - Missing Canonical_SMILES column\n")
            return None, None
        
        # Record original row count
        original_count = len(df)
        
        # Check if SMILES contain metal elements
        mask = df['Canonical_SMILES'].apply(
            lambda x: bool(metal_pattern.search(str(x))) if pd.notna(x) else False
        )
        
        # Create new column marking if contains metals
        df['contains_metal'] = mask
        
        # Get removed metal compounds
        removed_metals = df[mask].copy()
        
        # Filter out rows containing metals
        df_filtered = df[~mask].copy()
        
        # Record processed row count
        filtered_count = len(df_filtered)
        removed_count = original_count - filtered_count
        
        # Create output filenames
        base_name = os.path.splitext(filename)[0]
        output_path = os.path.join(output_dir, f"{base_name}_filtered.csv")
        removed_path = os.path.join(output_dir, f"{base_name}_removed_metals.csv")
        stats_path = os.path.join(output_dir, f"{base_name}_stats.txt")
        
        # Save processed file
        df_filtered.to_csv(output_path, index=False)
        
        # Save removed metal compounds
        if not removed_metals.empty:
            removed_metals.to_csv(removed_path, index=False)
        
        # Generate statistics
        stats = {
            'Filename': filename,
            'Processing Time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'Original Compounds': original_count,
            'Filtered Compounds': filtered_count,
            'Removed Metal Compounds': removed_count,
            'Removal Rate': f"{(removed_count/original_count)*100:.2f}%" if original_count > 0 else "0%",
            'Metal Elements Found': list(set(
                metal for smi in removed_metals['Canonical_SMILES'].astype(str) 
                for metal in metal_elements 
                if metal.lower() in smi.lower()
            )) if not removed_metals.empty else []
        }
        
        # Save statistics
        with open(stats_path, 'w') as f:
            for key, value in stats.items():
                f.write(f"{key}: {value}\n")
        
        # Log entry
        log_entry = (
            f"{datetime.now()}: {filename}\n"
            f"  Original rows: {original_count}, Filtered rows: {filtered_count}\n"
            f"  Removed rows: {removed_count}, Removal rate: {stats['Removal Rate']}\n"
            f"  Metal elements found: {stats['Metal Elements Found']}\n"
            f"  Output file: {output_path}\n"
            f"  Removed records: {removed_path if not removed_metals.empty else 'None'}\n\n"
        )
        
        with open(log_file, 'a') as f:
            f.write(log_entry)
        
        return stats, removed_metals
    
    except Exception as e:
        error_msg = f"{datetime.now()}: Error processing file {filename}: {str(e)}\n"
        with open(log_file, 'a') as f:
            f.write(error_msg)
        return None, None

def main():
    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Filter compounds containing metal elements from SMILES')
    parser.add_argument('-i', '--input', required=True, help='Input folder path')
    parser.add_argument('-o', '--output', required=True, help='Output folder path')
    args = parser.parse_args()
    
    # Create output folder
    os.makedirs(args.output, exist_ok=True)
    
    # Create log file paths
    log_file = os.path.join(args.output, 'processing_log.txt')
    summary_file = os.path.join(args.output, 'summary_report.txt')
    
    # Write log header
    with open(log_file, 'w') as f:
        f.write(f"Metal Compound Filtering Processing Log\n")
        f.write(f"Start Time: {datetime.now()}\n")
        f.write(f"Input Folder: {args.input}\n")
        f.write(f"Output Folder: {args.output}\n")
        f.write("="*50 + "\n\n")
    
    # Get all CSV files
    csv_files = [f for f in os.listdir(args.input) if f.endswith('.csv')]
    total_files = len(csv_files)
    
    # Initialize summary statistics
    all_stats = []
    total_original = 0
    total_filtered = 0
    
    # Process each file
    for csv_file in tqdm(csv_files, desc="Processing files"):
        input_path = os.path.join(args.input, csv_file)
        stats, _ = process_file(input_path, args.output, log_file)
        if stats:
            all_stats.append(stats)
            total_original += stats['Original Compounds']
            total_filtered += stats['Filtered Compounds']
    
    # Generate summary report
    total_removed = total_original - total_filtered
    removal_rate = (total_removed/total_original)*100 if total_original > 0 else 0
    
    summary = (
        f"\n{'='*50}\n"
        f"Summary Report\n"
        f"{'='*50}\n"
        f"Processing Completed: {datetime.now()}\n"
        f"Total Files Processed: {total_files}\n"
        f"Total Compounds: {total_original}\n"
        f"Filtered Compounds: {total_filtered}\n"
        f"Total Metal Compounds Removed: {total_removed}\n"
        f"Overall Removal Rate: {removal_rate:.2f}%\n"
        f"\nProcessing Details by File:\n"
    )
    
    for stats in all_stats:
        summary += (
            f"\nFile: {stats['Filename']}\n"
            f"  Original: {stats['Original Compounds']}\n"
            f"  Filtered: {stats['Filtered Compounds']}\n"
            f"  Removed: {stats['Removed Metal Compounds']} ({stats['Removal Rate']})\n"
            f"  Metal Elements: {stats['Metal Elements Found']}\n"
        )
    
    with open(summary_file, 'w') as f:
        f.write(summary)
    
    # Complete log
    with open(log_file, 'a') as f:
        f.write(f"\n{'='*50}\n")
        f.write(f"Processing completed at: {datetime.now()}\n")
        f.write(f"Total files processed: {total_files}\n")
        f.write(f"Summary report saved to: {summary_file}\n")

if __name__ == "__main__":
    main()
