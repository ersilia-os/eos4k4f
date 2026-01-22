# imports
import os
import csv
import sys
from rdkit import Chem
import datamol as dm
from chembl_structure_pipeline import standardizer
from rdkit.Chem.Scaffolds import MurckoScaffold


def get_canonical_smiles_datamol(smiles):
    try:
        mol = dm.to_mol(smiles)
        mol = dm.fix_mol(mol)
        mol = dm.sanitize_mol(mol)
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        return mol, smiles
    except:
        return None, ""

def get_canonical_smiles_rdkit(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        return mol, smiles
    except:
        return None, ""

def get_canonical_smiles(smiles):
    smiles = str(smiles).strip()
    try:
        mol, canonical_smiles = get_canonical_smiles_datamol(smiles)
        return mol, canonical_smiles
    except Exception:
        pass
    try:
        mol, canonical_smiles = get_canonical_smiles_rdkit(smiles)
        return mol, canonical_smiles
    except Exception:
        return None, ""
    
def get_standardized_smiles(mol):
    if mol is None:
        return None, ""
    try:
        mol, _ = standardizer.get_parent_mol(mol)
        mol = standardizer.standardize_mol(mol)
        standardized_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        return mol, standardized_smiles
    except:
        return None, ""
    
def get_flattened_smiles(mol):
    if mol is None:
        return None, ""
    try:
        flattened_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        mol = Chem.MolFromSmiles(flattened_smiles)
        return mol, flattened_smiles
    except:
        return None, ""
    
def get_murcko_scaffold(mol):
    if mol is None:
        return None, ""
    try:
        mol_scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        smiles_scaffold = Chem.MolToSmiles(mol_scaffold, canonical=True)
        return mol_scaffold, smiles_scaffold
    except:
        return None, ""
    
def get_generic_scaffold(mol_scaffold):
    if mol_scaffold is None:
        return None, ""
    try:
        mol_generic_scaffold = MurckoScaffold.MakeScaffoldGeneric(mol_scaffold)
        smiles_generic_scaffold = Chem.MolToSmiles(mol_generic_scaffold, canonical=True)
        return mol_generic_scaffold, smiles_generic_scaffold
    except:
        return None, ""

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# run model
outputs = []

for smiles in smiles_list:

    mol, canonical_smiles = get_canonical_smiles(smiles)
    mol, standardized_smiles = get_standardized_smiles(mol)
    mol, flattened_smiles = get_flattened_smiles(mol)
    mol_scaffold, smiles_scaffold = get_murcko_scaffold(mol)
    mol_generic_scaffold, smiles_generic_scaffold = get_generic_scaffold(mol_scaffold)

    # Store results
    outputs.append([canonical_smiles, standardized_smiles, flattened_smiles, smiles_scaffold, smiles_generic_scaffold])


#check input and output have the same lenght
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["canonical_smiles", "standardized_smiles", "flattened_smiles", "murcko_scaffold", "generic_scaffold"])  # header
    for o in outputs:
        writer.writerow(o)



