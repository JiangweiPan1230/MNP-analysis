# Ro5

A Python tool for analyzing molecular datasets and counting violations of the Rule of Five (Lipinski's Rule).

## Overview

The Rule of Five is a rule of thumb to evaluate drug-likeness and determine if a chemical compound has properties that would make it a likely orally active drug in humans. This tool analyzes CSV files containing molecular data and counts violations of these rules.

## Rule of Five Criteria

The tool checks for violations of the following criteria:
- **Lipophilicity**: ALOGP ≤ 5
- **Molecular Weight**: MW ≤ 500
- **Hydrogen Bond Acceptors**: HBA ≤ 10
- **Hydrogen Bond Donors**: HBD ≤ 5
- **Polar Surface Area**: PSA ≤ 140 Å²
- **Rotatable Bonds**: ROTB ≤ 10

## Requirements

- Python 3.6+
- Required packages:
  - tqdm

Install dependencies:
```bash
pip install tqdm
```



# Canonical_Kekulize

A Python tool for processing SMILES strings in CSV files and converting them to canonical forms, with optional Kekulé structure representation.

## Overview

This tool processes molecular data in CSV files by converting SMILES strings to their canonical forms using RDKit. It includes options for generating Kekulé structures and handles invalid molecules gracefully.

## Features

- **Canonical SMILES Generation**: Convert SMILES to standardized canonical form
- **Kekulé Structure Support**: Optional generation of Kekulé-form SMILES
- **Batch Processing**: Process single files or entire folders
- **Validation**: Comprehensive molecule validation including bond checking
- **Error Handling**: Robust error handling with detailed logging
- **Progress Tracking**: Visual progress bar for batch processing

## Requirements

- Python 3.6+
- Required packages:
  - pandas
  - rdkit
  - tqdm

Install dependencies:
```bash
pip install pandas rdkit-pypi tqdm
```



# Clear Charge

A Python tool for standardizing SMILES strings using RDKit's MolStandardize module, with optional generation of Kekulé structures.

## Overview

This tool processes molecular data in CSV files and standardizes SMILES strings using RDKit's molecular standardization capabilities. It generates Kekulé structures by default and provides various standardization options.

## Features

- **Charge Neutralization**: Remove charges from molecules (enabled by default)
- **Fragment Cleanup**: Optional removal of small fragments
- **Tautomer Standardization**: Optional canonical tautomer generation
- **Kekulé Structure**: Generates SMILES with explicit bond orders
- **Stereochemistry Preservation**: Optional retention of stereochemical information
- **Batch Processing**: Process single files or entire folders
- **Progress Tracking**: Visual progress bars for processing
- **Error Handling**: Robust error handling with detailed warnings

## Standardization Steps

### Mandatory Steps (Always Performed)
- Basic molecular cleanup using RDKit's `Cleanup` function

### Optional Steps (Parameter-Controlled)
- **Charge Neutralization**: Removes charges from molecules (default: enabled)
- **Fragment Cleanup**: Keeps only the largest fragment (default: disabled)
- **Tautomer Standardization**: Generates canonical tautomer (default: disabled)
- **Kekulé Structure**: Always generates SMILES with explicit bond orders

## Requirements

- Python 3.6+
- Required packages:
  - pandas
  - rdkit
  - numpy
  - tqdm

Install dependencies:
```bash
pip install pandas rdkit-pypi numpy tqdm
```



# Get_TNP

A Python tool for comparing two SMILES datasets and extracting unique molecules from the first dataset.

## Overview

This tool compares two molecular datasets containing canonical SMILES strings and identifies molecules that are present in the first dataset but not in the second dataset. It's particularly useful for finding unique natural products or compounds in different databases.

## Features

- **Set Operations**: Uses Python set operations for efficient comparison
- **Large Dataset Handling**: Efficiently handles large datasets using set operations
- **Data Validation**: Checks for required columns and handles missing values
- **Simple Output**: Generates clean CSV output with unique SMILES

## Use Cases

- Compare natural product databases (COCONUT vs CMNPD)
- Identify unique compounds in different chemical libraries
- Find novel molecules in one dataset compared to another
- Database deduplication and uniqueness analysis

## Requirements

- Python 3.6+
- pandas

Install dependencies:
```bash
pip install pandas
```



# Get_Scaffold (Murcko Scaffold Extractor)

A Python tool for extracting Murcko scaffolds from SMILES strings in CSV files, with optional Kekulé structure generation.

## Overview

This tool processes molecular data in CSV files and extracts Murcko scaffolds using RDKit. Murcko scaffolds represent the core structure of molecules by removing side chains and retaining the ring systems and linkers. The tool can generate scaffolds with explicit Kekulé bond orders.

## Features

- **Murcko Scaffold Extraction**: Extract core molecular frameworks using RDKit's MurckoScaffold
- **Kekulé Structure Support**: Optional generation of scaffolds with explicit bond orders
- **Batch Processing**: Process single files or entire folders of CSV files
- **Flexible Input**: Customizable SMILES column name
- **Progress Tracking**: Visual progress bars for processing
- **Error Handling**: Robust handling of invalid SMILES and file errors
- **Data Cleaning**: Automatically removes rows where scaffold extraction fails

## What are Murcko Scaffolds?

Murcko scaffolds represent the core framework of a molecule by:
- Retaining all ring systems
- Keeping the linkers between rings
- Removing all side chains and substituents
- Providing a simplified view of molecular architecture

## Requirements

- Python 3.6+
- Required packages:
  - pandas
  - rdkit
  - tqdm

Install dependencies:
```bash
pip install pandas rdkit-pypi tqdm
```



# Get_Ring_System(Ring System Extractor)

A Python tool for extracting individual ring systems from molecules and generating Kekulé SMILES representations.

Reference:[[Natural Product Reports, 2022, DOI: 10.1039/D2NP00001F](https://doi.org/10.1039/D2NP00001F)]

## Overview

This tool processes molecular data in CSV files and extracts all individual ring systems from each molecule. Unlike Murcko scaffolds that provide the core framework, this tool identifies and separates each distinct ring system, making it particularly useful for analyzing complex polycyclic compounds and natural products.

## Features

- **Ring System Extraction**: Identifies and separates individual ring systems from molecules
- **Kekulé Structure Generation**: Outputs ring systems with explicit bond orders
- **Dynamic Column Generation**: Automatically creates columns for the maximum number of ring systems found
- **Batch Processing**: Handles single files or entire directories of CSV files
- **Progress Tracking**: Visual progress bars for both file and molecule processing
- **Error Handling**: Robust error handling with detailed error messages
- **Flexible Input**: Customizable SMILES column name

## What are Ring Systems?

Ring systems are individual cyclic components of a molecule:
- Each connected set of rings is treated as a separate system
- Includes both isolated rings and fused ring systems
- Excludes non-ring atoms and side chains
- Useful for analyzing complex polycyclic structures

## Requirements

- Python 3.6+
- Required packages:
  - pandas
  - rdkit
  - tqdm

Install dependencies:
```bash
pip install pandas rdkit-pypi tqdm
```



# Get Substitute

A Python tool for decomposing molecules into their Murcko scaffolds and substituents using RDKit.

## Overview

This tool processes molecular data in CSV files and systematically separates each molecule into its core scaffold and individual substituents. It's particularly useful for structure-activity relationship (SAR) studies, scaffold-based analysis, and understanding molecular decoration patterns.

## Features

- **Scaffold Extraction**: Identifies Murcko scaffolds (core ring systems with linkers)
- **Substituent Detection**: Automatically detects and separates side chains and functional groups
- **Dynamic Column Generation**: Creates columns for all substituents found across the dataset
- **Batch Processing**: Handles single files or entire directories of CSV files
- **Progress Tracking**: Visual progress bars for processing
- **Error Resilience**: Continues processing even with invalid SMILES strings
- **Flexible Input**: Customizable SMILES column name

## What are Scaffolds and Substituents?

- **Murcko Scaffold**: The core framework of a molecule including ring systems and linkers between them
- **Substituents**: Side chains, functional groups, and decorations attached to the scaffold
- **Application**: Useful for analyzing molecular decoration patterns and scaffold-based drug design

## Requirements

- Python 3.6+
- Required packages:
  - pandas
  - rdkit
  - tqdm

Install dependencies:
```bash
pip install pandas rdkit-pypi tqdm
```



# Metal Compound Filter

A Python tool for filtering out compounds containing metal elements from SMILES datasets.

## Overview

This tool processes CSV files containing molecular data and removes compounds that contain metal elements based on their SMILES representations. It's particularly useful for preparing datasets for organic chemistry applications, drug discovery pipelines, and computational chemistry studies where metal-containing compounds are not desired.

## Features

- **Metal Detection**: Identifies compounds containing any of 70+ metal elements
- **Batch Processing**: Processes entire folders of CSV files automatically
- **Comprehensive Logging**: Detailed processing logs and statistics
- **Data Preservation**: Saves both filtered data and removed metal compounds
- **Progress Tracking**: Visual progress bars for file processing
- **Statistical Reporting**: Generates detailed statistics for each file and overall summary

## Supported Metal Elements

The tool detects the following metal elements:
- Alkali metals: Li, Na, K, Rb, Cs, Fr
- Alkaline earth metals: Be, Mg, Ca, Sr, Ba, Ra
- Transition metals: Ti, V, Cr, Mn, Fe, Ni, Cu, Zn, etc.
- Lanthanides: La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, etc.
- Actinides: Ac, Th, Pa, U, Pu, Am, Cm, etc.
- Post-transition metals: Al, Ga, In, Tl, Pb, Bi
- And many more (70+ elements total)

## Requirements

- Python 3.6+
- Required packages:
  - pandas
  - tqdm

Install dependencies:
```bash
pip install pandas tqdm