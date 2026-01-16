# imports
import os
import csv
import sys
from rdkit import Chem
import datamol as dm
from chembl_structure_pipeline import standardizer
from rdkit.Chem.Scaffolds import MurckoScaffold


def get_canonical_smiles_datamol(smiles):
    mol = dm.to_mol(smiles)
    if not mol:
        return None
    mol = dm.fix_mol(mol)
    mol = dm.sanitize_mol(mol)
    if not mol:
        return None
    smi = dm.to_smiles(mol)
    if not smi:
        return None
    m2 = Chem.MolFromSmiles(smi)
    if not m2:
        return None
    return Chem.MolToSmiles(m2, canonical=True, isomericSmiles=True)

def get_canonical_smiles_rdkit(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

def get_canonical_smiles(smiles):
    smiles = str(smiles).strip()
    try:
        canonical_smiles = get_canonical_smiles_datamol(smiles)
        if canonical_smiles:
            return canonical_smiles
    except Exception:
        pass
    try:
        return get_canonical_smiles_rdkit(smiles)
    except Exception:
        return None
    
def get_standardized_smiles(smiles):
    if not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    if not m:
        return None
    m = standardizer.standardize_mol(m)
    m, _ = standardizer.get_parent_mol(m)
    if not m:
        return None
    m = standardizer.standardize_mol(m)
    if not m:
        return None
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)

def get_flattened_smiles(smiles):
    if not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    if not m:
        return None
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=False)
         
def get_murcko_scaffold(smiles):
    if not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    if not m:
        return None
    scaf = MurckoScaffold.GetScaffoldForMol(m)
    if not scaf:
        return None
    return Chem.MolToSmiles(scaf, canonical=True)
    
def get_generic_scaffold(scaf_smiles):
    if not scaf_smiles:
        return None
    scaf = Chem.MolFromSmiles(scaf_smiles)
    if not scaf:
        return None
    generic = MurckoScaffold.MakeScaffoldGeneric(scaf)
    if not generic:
        return None
    return Chem.MolToSmiles(generic, canonical=True)

def deal_with_nones(string):
    if string is None:
        return ""
    return string


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

    canonical_smiles = get_canonical_smiles(smiles)
    standardized_smiles = get_standardized_smiles(canonical_smiles)
    flattened_smiles = get_flattened_smiles(standardized_smiles)
    murcko_scaffold = get_murcko_scaffold(flattened_smiles)
    generic_scaffold = get_generic_scaffold(murcko_scaffold)

    # Deal with None's
    canonical_smiles = deal_with_nones(canonical_smiles)
    standardized_smiles = deal_with_nones(standardized_smiles)
    flattened_smiles = deal_with_nones(flattened_smiles)
    murcko_scaffold = deal_with_nones(murcko_scaffold)
    generic_scaffold = deal_with_nones(generic_scaffold)

    # Store results
    outputs.append([canonical_smiles, standardized_smiles, flattened_smiles, murcko_scaffold, generic_scaffold])


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



