#!/usr/bin/env python3
import csv
import os
import argparse
from collections import defaultdict
from tqdm import tqdm

# Rule of Five descriptor list
RULE_OF_FIVE = ['ALOGP', 'HBA', 'HBD', 'PSA', 'ROTB']

def count_rule_violations(row):
    """Count the number of Rule of Five violations for a single molecule"""
    violations = 0
    
    # Lipophilicity (ALOGP <= 5)
    try:
        if float(row.get('ALOGP', 0)) > 5:
            violations += 1
    except (ValueError, TypeError):
        pass
    
    # Molecular Weight (MW <= 500) - Usually needs calculation, assuming MW column exists
    try:
        if float(row.get('MW', 0)) > 500:
            violations += 1
    except (ValueError, TypeError):
        pass
    
    # Hydrogen Bond Acceptors (HBA <= 10)
    try:
        if int(row.get('HBA', 0)) > 10:
            violations += 1
    except (ValueError, TypeError):
        pass
    
    # Hydrogen Bond Donors (HBD <= 5)
    try:
        if int(row.get('HBD', 0)) > 5:
            violations += 1
    except (ValueError, TypeError):
        pass
    
    # Polar Surface Area (PSA <= 140)
    try:
        if float(row.get('PSA', 0)) > 140:
            violations += 1
    except (ValueError, TypeError):
        pass
    
    # Rotatable Bonds (ROTB <= 10)
    try:
        if int(row.get('ROTB', 0)) > 10:
            violations += 1
    except (ValueError, TypeError):
        pass
    
    return violations

def process_csv_file(filepath):
    """Process a single CSV file"""
    violation_counts = defaultdict(int)
    total_molecules = 0
    
    try:
        with open(filepath, 'r', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            # Check if required columns exist
            required_columns = ['SMILES'] + RULE_OF_FIVE
            missing_cols = [col for col in required_columns if col not in reader.fieldnames]
            if missing_cols:
                print(f"Warning: File {os.path.basename(filepath)} is missing required columns: {', '.join(missing_cols)}")
                return None
            
            for row in reader:
                if not row['SMILES'].strip():
                    continue  # Skip empty SMILES
                
                total_molecules += 1
                violations = count_rule_violations(row)
                violation_counts[violations] += 1
    
    except Exception as e:
        print(f"Error processing file {os.path.basename(filepath)}: {str(e)}")
        return None
    
    return {
        'filename': os.path.basename(filepath),
        'total_molecules': total_molecules,
        'violation_counts': violation_counts
    }

def process_folder(input_folder, output_file):
    """Process entire folder"""
    csv_files = [f for f in os.listdir(input_folder) if f.endswith('.csv')]
    if not csv_files:
        print(f"No CSV files found in folder {input_folder}")
        return
    
    results = []
    
    # Use tqdm to show progress bar
    for filename in tqdm(csv_files, desc="Processing files"):
        filepath = os.path.join(input_folder, filename)
        result = process_csv_file(filepath)
        if result:
            results.append(result)
    
    if not results:
        print("No valid statistical results")
        return
    
    # Write summary results
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header row
        header = ['Filename', 'Total Molecules'] + [f'{i}_violations' for i in range(6)]
        writer.writerow(header)
        
        for result in results:
            row = [
                result['filename'],
                result['total_molecules']
            ]
            # Add counts for 0-5 violations
            row.extend([result['violation_counts'].get(i, 0) for i in range(6)])
            writer.writerow(row)
    
    print(f"\nProcessing completed! Results saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Count Rule of Five violations for molecules in CSV files within a folder',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True,
                       help='Input folder path containing CSV files')
    parser.add_argument('-o', '--output', required=True,
                       help='Output CSV file path for statistical results')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input):
        print(f"Error: Input path {args.input} is not a valid folder")
        return
    
    process_folder(args.input, args.output)

if __name__ == '__main__':
    main()
