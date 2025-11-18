#!/usr/bin/env python
import argparse
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from glob import glob
import sys

class MoleculeWithRings:
    '''
    MoleculeWithRings class contains a molecule and the corresponding ring systems
    '''
    def __init__(self, mol):
        self.mol = mol
        self.ringList = self._get_ring_system_list()

    def get_ring_system_dict(self):
        ringDict = {}
        for ringSystem in self.ringList:
            ringDict[ringSystem] = self._get_smiles_of_ring_system_without_dummy_atoms(ringSystem)  
        return ringDict

    def _get_ring_system_list(self):
        ringAtomsSets = self._get_list_of_ring_atom_sets()
        for ringAtomsSet in ringAtomsSets:
            ringAtomsSet = self._extend_ringSet(ringAtomsSet)
        ringMolList = self._get_fragments_from_mol(ringAtomsSets)
        ringMolListWithoutDummyAtom = self._get_ring_mol_without_dummy_atom(ringMolList)
        return ringMolListWithoutDummyAtom

    def _get_list_of_ring_atom_sets(self):
        ringInfo = self.mol.GetRingInfo()
        listOfRingAtomSets = []
        for atomRing in ringInfo.AtomRings():
            ringAtomSet = set(atomRing)
            newList = []
            for atomSet in listOfRingAtomSets:
                BridgeAtomSet = ringAtomSet.intersection(atomSet)
                if len(BridgeAtomSet):
                    ringAtomSet = (ringAtomSet.union(atomSet))
                else:
                    newList.append(atomSet)
            newList.append(ringAtomSet)
            listOfRingAtomSets = newList
        return listOfRingAtomSets

    def _extend_ringSet(self, ringAtomsSet):
        extendSet = set()
        for moleculeAtom in self.mol.GetAtoms():
            atomIdx = moleculeAtom.GetIdx()
            if atomIdx in ringAtomsSet:
                continue
            else:
                neighborAtoms = moleculeAtom.GetNeighbors()
            for neighbor in neighborAtoms:
                if (neighbor.GetIdx() in ringAtomsSet) and neighbor.IsInRing():
                    if self._is_not_single_bond(moleculeAtom, neighbor):
                        extendSet.add(moleculeAtom.GetIdx())
        ringAtomsSet.update(extendSet)
        return ringAtomsSet

    def _is_not_single_bond(self, atom, neighbor):
        atomIdx = atom.GetIdx()
        neighborIdx = neighbor.GetIdx()
        bond = self.mol.GetBondBetweenAtoms(atomIdx, neighborIdx)
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return True
        else:
            return False

    def _get_fragments_from_mol(self, fragmentAtomSetsList):
        ringMolList = []
        for fragmentSet in fragmentAtomSetsList:
            breakpoints = set()
            markingAtomInRing = self.mol.GetAtomWithIdx(list(fragmentSet)[0])
            Chem.SetSupplementalSmilesLabel(markingAtomInRing, '<xxx>')
            for bond in self.mol.GetBonds():
                begin = bond.GetBeginAtomIdx()
                end = bond.GetEndAtomIdx()
                breakpointIdx = self._is_bond_breakpoint(bond, begin, end, fragmentSet)
                if breakpointIdx is not None:
                    breakpoints.add(breakpointIdx)
            if breakpoints:
                breakpoints = list(breakpoints)
                fragments = Chem.FragmentOnBonds(self.mol, breakpoints)
                Chem.SanitizeMol(fragments, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                ringFragmentList = Chem.GetMolFrags(fragments, asMols=True, sanitizeFrags=False)
                for molFragment in ringFragmentList:
                    for atom in molFragment.GetAtoms():
                        if atom.GetIsAromatic() and not atom.IsInRing():
                            atom.SetIsAromatic(False)
                        if Chem.GetSupplementalSmilesLabel(atom) == '<xxx>':
                            Chem.SetSupplementalSmilesLabel(atom, "")
                            Chem.SetSupplementalSmilesLabel(markingAtomInRing, '')
                            ringMolList.append(molFragment)
            else:
                Chem.SetSupplementalSmilesLabel(markingAtomInRing, '')
                ringMolList.append(self.mol)
                return ringMolList
        return ringMolList

    @staticmethod
    def _is_bond_breakpoint(bond, beginAtom, endAtom, fragmentSet):
        if (beginAtom in fragmentSet and endAtom not in fragmentSet) or (
                beginAtom not in fragmentSet and endAtom in fragmentSet):
            return bond.GetIdx()

    def _get_ring_mol_without_dummy_atom(self, ringMolList): 
        ringMolListWithoutDummyAtom = list()
        for ring in ringMolList:
            ringWithoutDummyAtom = AllChem.ReplaceSubstructs(ring,
                                                         Chem.MolFromSmiles('*'),
                                                         Chem.MolFromSmiles('[H]'),
                                                         True)
            ringWithoutDummyAtom = ringWithoutDummyAtom[0]

            for atom in ringWithoutDummyAtom.GetAtoms():
                if atom.GetIsAromatic() and not atom.IsInRing():
                    atom.SetIsAromatic(False)
                    neighbors = atom.GetNeighbors()
                    for neighbor in neighbors:
                        bond = ringWithoutDummyAtom.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetIsAromatic():
                            bond.SetIsAromatic(False) 
            smiles = Chem.MolToSmiles(ringWithoutDummyAtom)
            ringWithoutDummyAtom = Chem.MolFromSmiles(smiles) 

            ringMolListWithoutDummyAtom.append(ringWithoutDummyAtom)
        return ringMolListWithoutDummyAtom
    
    @staticmethod
    def _get_smiles_of_ring_system_without_dummy_atoms(ringSystem):
        if ringSystem is not None:
            smilesWithoutDummyAtoms = Chem.MolToSmiles(ringSystem)
            return smilesWithoutDummyAtoms
        return None

class MoleculeProcessor:
    def __init__(self):
        pass

    def process_molecule(self, mol):
        if not mol:
            return []
        
        try:
            ring_system_finder = MoleculeWithRings(mol)
            ring_dict = ring_system_finder.get_ring_system_dict()
            return list(ring_dict.values())
        except Exception as e:
            print(f"Error processing molecule: {e}", file=sys.stderr)
            return []

    def kekulize_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            try:
                Chem.Kekulize(mol)
                return Chem.MolToSmiles(mol, kekuleSmiles=True)
            except:
                return smiles
        return smiles

def process_file(input_file, output_file, smiles_column='SMILES'):
    processor = MoleculeProcessor()
    
    try:
        df = pd.read_csv(input_file)
        
        # Validate SMILES column exists
        if smiles_column not in df.columns:
            print(f"Error: Column '{smiles_column}' not found in {input_file}", file=sys.stderr)
            return False

        # Check and clean SMILES data
        if df[smiles_column].isnull().any():
            print(f"Warning: Found empty SMILES in {input_file}, filling with ''", file=sys.stderr)
            df[smiles_column] = df[smiles_column].fillna('')
        
        all_results = []
        max_rings = 0
        
        # First pass to determine maximum number of rings
        for smiles in df[smiles_column]:
            if not isinstance(smiles, str) or not smiles.strip():
                continue
                
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                ring_systems = processor.process_molecule(mol)
                max_rings = max(max_rings, len(ring_systems))
        
        print(f"Maximum ring systems found: {max_rings}", file=sys.stderr)
        
        # Second pass to process all molecules
        for idx, smiles in tqdm(df[smiles_column].items(), total=len(df), desc=f"Processing {os.path.basename(input_file)}"):
            result = {'Original_SMILES': smiles, 'RingSystems_count': 0}
            
            # Initialize all possible ring system columns
            for i in range(1, max_rings + 1):
                result[f'RingSystems{i}'] = ''
                
            if isinstance(smiles, str) and smiles.strip():
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    ring_systems = processor.process_molecule(mol)
                    kekule_ring_systems = [processor.kekulize_smiles(rs) for rs in ring_systems if rs]
                    
                    result['RingSystems_count'] = len(kekule_ring_systems)
                    for i, rs in enumerate(kekule_ring_systems, 1):
                        result[f'RingSystems{i}'] = rs
            
            all_results.append(result)

        # Create DataFrame with consistent columns
        result_df = pd.DataFrame(all_results)
        
        # Reorder columns: Original_SMILES, RingSystems_count, then RingSystems1..N
        columns = ['Original_SMILES', 'RingSystems_count'] + \
                 [f'RingSystems{i}' for i in range(1, max_rings + 1)]
        result_df = result_df[columns]
        
        result_df.to_csv(output_file, index=False)
        print(f"Successfully processed {len(result_df)} molecules to {output_file}", file=sys.stderr)
        return True
        
    except Exception as e:
        print(f"Fatal error processing file {input_file}: {str(e)}", file=sys.stderr)
        return False

def main():
    parser = argparse.ArgumentParser(description='Extract ring systems from molecules and output Kekule SMILES.')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file or directory containing CSV files')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file or directory')
    parser.add_argument('-s', '--smiles_column', default='SMILES', help='Column name containing SMILES (default: SMILES)')
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    if os.path.isdir(args.input) and not os.path.exists(args.output):
        os.makedirs(args.output)

    # Process based on input type
    if os.path.isdir(args.input):
        # Directory mode
        input_files = glob(os.path.join(args.input, '*.csv'))
        processed_count = 0
        
        for input_file in tqdm(input_files, desc="Processing files"):
            output_file = os.path.join(args.output, os.path.basename(input_file))
            if process_file(input_file, output_file, args.smiles_column):
                processed_count += 1
        
        print(f"\nProcessed {processed_count}/{len(input_files)} files successfully.", file=sys.stderr)
    else:
        # Single file mode
        if os.path.isdir(args.output):
            output_file = os.path.join(args.output, os.path.basename(args.input))
        else:
            output_file = args.output
        
        if process_file(args.input, output_file, args.smiles_column):
            print(f"Successfully processed {args.input} -> {output_file}", file=sys.stderr)
        else:
            sys.exit(1)

if __name__ == '__main__':
    main()
