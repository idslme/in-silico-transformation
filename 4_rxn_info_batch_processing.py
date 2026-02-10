import os
import re
import subprocess
from multiprocessing import Pool
import pandas as pd
import rxn_insight
from typing import Dict
from pyteomics.mass import Composition
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

########### Checking Reaction Stoichiometry Balance of the Reaction ##########
def _normalize_rxn_smiles(rxn_smiles: str) -> str:
    s = (rxn_smiles or "").strip()
    return s.replace("&gt;&gt;", ">>")

def _split_side(smiles: str):
    parts = [p for p in smiles.split(".") if p]
    mols = [Chem.MolFromSmiles(p) for p in parts]
    if not mols or any(m is None for m in mols):
        raise ValueError(f"Invalid SMILES in side: {smiles}")
    return mols

# Precompiled regex for trailing charge annotations: '+', '-', '+2', '-3', etc.
_TRAILING_CHARGE_RE = re.compile(r'([+-]\d*)$')

def _composition_for_side(mols) -> Composition:
    """
    Build a pyteomics Composition for a list of RDKit Mol objects.
    Defensively strips any trailing charge annotation from the formula string
    before passing it to Composition.
    """
    total = Composition()
    for m in mols:
        f = rdMolDescriptors.CalcMolFormula(m)  # RDKit formula (neutral-only in practice)
        # Defensive: remove any trailing charge annotation if present
        # e.g., 'C5H8NO4-', 'Na+', 'C7H6O4-2' -> 'C5H8NO4', 'Na', 'C7H6O4'
        f_neutral = _TRAILING_CHARGE_RE.sub('', f)
        total += Composition(f_neutral)
    return total

def _formal_charge_for_side(mols) -> int:
    return sum(rdmolops.GetFormalCharge(m) for m in mols)

def assess_reaction_balance(rxn_id: str, rxn_smiles: str) -> Dict[str, object]:
    """
    Assess elemental and formal charge balance of a reaction SMILES, and indicate
    which side (reactants/products) is short of atoms/charge.

    Returns a dict:
      {
        "normalized_rxn": "...",
        "reactants": {"composition": {...}, "formal_charge": int},
        "products":  {"composition": {...}, "formal_charge": int},
        "atom_difference": {...},  # products - reactants (signed)
        "mass_status": "balanced | mass_imbalance",
        "missing_atoms_side": "none | reactants | products",
        "missing_atoms_detail": {"C": 1, "H": 2},  # counts on the missing side
        "charge_difference": int,  # products - reactants
        "charge_status": str       # human-readable summary
      }
    """
    s = _normalize_rxn_smiles(rxn_smiles)
    if ">>" not in s:
        raise ValueError("Invalid reaction SMILES format. Expected '>>' delimiter.")

    reactants_str, products_str = s.split(">>")
    reactant_mols = _split_side(reactants_str)
    product_mols  = _split_side(products_str)

    # Elemental compositions
    comp_react = _composition_for_side(reactant_mols)
    comp_prod  = _composition_for_side(product_mols)

    # Convert to plain dicts for readability
    comp_react_dict = dict(comp_react)
    comp_prod_dict  = dict(comp_prod)

    # Signed difference: products - reactants
    diff = comp_prod - comp_react
    diff_dict = dict(diff)

    # Mass status and which side is missing atoms
    if diff_dict:
        mass_status = "mass_imbalance"
        # If products - reactants is positive for an element, products have *more* of that element.
        # Missing side is the one with fewer atoms for elements where the absolute deficit exists.
        # We can aggregate a side classification as follows:
        #   If any element has diff>0 -> reactants missing those counts.
        #   If any element has diff<0 -> products missing those counts.
        # If both signs appear, we report both sides being short (common in truncated SMILES).
        missing_reactants = {el: cnt for el, cnt in diff_dict.items() if cnt > 0}
        missing_products  = {el: -cnt for el, cnt in diff_dict.items() if cnt < 0}

        if missing_reactants and missing_products:
            missing_atoms_side   = "both"
            missing_atoms_detail = {"reactants": missing_reactants, "products": missing_products}
        elif missing_reactants:
            missing_atoms_side   = "reactants"
            missing_atoms_detail = missing_reactants
        elif missing_products:
            missing_atoms_side   = "products"
            missing_atoms_detail = missing_products
        else:
            # Should not occur since diff_dict is non-empty, but keep a safe default.
            missing_atoms_side   = "none"
            missing_atoms_detail = {}
    else:
        mass_status = "balanced"
        missing_atoms_side   = "none"
        missing_atoms_detail = {}

    # Formal charge balance
    charge_react = _formal_charge_for_side(reactant_mols)
    charge_prod  = _formal_charge_for_side(product_mols)
    charge_diff  = charge_prod - charge_react

    if charge_diff == 0:
        charge_status = "charge_balanced"
    elif charge_diff > 0:
        charge_status = f"products have +{charge_diff} more formal charge"
    else:
        charge_status = f"reactants have +{abs(charge_diff)} more formal charge"
    
    is_balanced = (mass_status == "balanced") and (charge_diff == 0)
    reaction_status = "balanced" if is_balanced else "imbalanced"


    return {
        "rxn_id": rxn_id, 
        "rxn": s,
        "reactants": {
            "composition": comp_react_dict,
            "formal_charge": charge_react,
        },
        "products": {
            "composition": comp_prod_dict,
            "formal_charge": charge_prod,
        },
        "atom_difference": diff_dict,
        "mass_status": mass_status,
        "missing_atoms_side": missing_atoms_side,
        "missing_atoms_detail": missing_atoms_detail,
        "charge_difference": charge_diff,
        "charge_status": charge_status,
        "reaction_status": reaction_status,
    }


def asses_reaction_balance_row(row):
    rxn_id = row["RXN_ID"]
    rxn_smi = row["SANITIZED_MAPPED_REACTION"]
    try:
        result= assess_reaction_balance(rxn_id,rxn_smi)
        return result
    except Exception as e:
        return {
            "RXN_ID": rxn_id,
            "normalized_rxn": rxn_smi,
            "reactants": None,
            "products": None,
            "atom_difference": None,
            "mass_status": "error",
            "missing_atoms_side": None,
            "missing_atoms_detail": None,
            "charge_difference": None,
            "charge_status": f"error: {type(e).__name__}: {e}",
            "is_balanced": False,
            "reaction_status": "error",
        }

# Set up directories and read the input CSV
os.makedirs('/home/mani/missing_reaction_processing/rxn_info_batchpro', exist_ok=True)
unprocess_rxns = pd.read_csv('/home/mani/missing_reaction_processing/missing_rxns.csv')
unprocess_rxns.rename(columns={'reaction_cid': 'rxn_id', 'reaction_smiles': 'rxn_smiles'}, inplace=True)
# Define the output directory and subset size
output_dir = '/home/mani/missing_reaction_processing/rxn_info_batchpro'
subset_size = 500

# Split the dataframe into smaller subsets and save them in different folders
for i in range(0, len(unprocess_rxns), subset_size):
    subset = unprocess_rxns.iloc[i:i + subset_size]
    subset_dir = os.path.join(output_dir, f"subset_{i // subset_size + 1}")
    os.makedirs(subset_dir, exist_ok=True)  # Create directory for each subset
    output_file = os.path.join(subset_dir, "data.csv")
    subset.to_csv(output_file, index=False)
    print(f"Saved {output_file}")

# Now, collect paths of all the generated subset files
target_dir = '/home/mani/missing_reaction_processing/rxn_info_batchpro'
sub_dir_paths = []
for sub_dir in os.listdir(target_dir):
    sub_dir_path = os.path.join(target_dir, sub_dir)
    if os.path.isdir(sub_dir_path):  # Ensure we only list directories
        sub_dir_paths.append(sub_dir_path)
       

file_paths = []
for sub_dir_path in sub_dir_paths:
    # Check if 'reaction_info_with_ids.csv' exists in the sub_dir_path
    if 'reaction_info_with_ids.csv' in os.listdir(sub_dir_path):
        print(f"Skipping folder: {sub_dir_path} as it contains reaction_info_with_ids.csv")
        continue  # Skip this folder and move on to the next one
    
    file_names = os.listdir(sub_dir_path)
    for file_name in file_names:
        if file_name.endswith('.csv'):  # Only process .csv files
            file_path = os.path.join(sub_dir_path, file_name)
            file_paths.append(file_path)

# Path to the reaction info generation script
rxn_info_generation = '/home/mani/ManiProject2025/ManiProject2024/pythonscripts-2025/Metabolite_Predicition_ExposureRX/rxn_info_cli.py'

# Iterate over all subset files and process them
for file_path in file_paths:
    input_file = file_path
    max_num_cpus = '16'
    output_dir = os.path.dirname(file_path)
    
    num_reactions = len(pd.read_csv(input_file))
    
    try:
        print(f"Processing file: {input_file}")
        print(f"Number of reactions to be processed: {num_reactions}")
        # Updated to positional arguments
        result = subprocess.run([
            'python', rxn_info_generation,
            input_file,
            max_num_cpus,
            output_dir
        ], capture_output=True, text=True)

        if result.returncode == 0:
            print(f"Successfully processed {input_file}")
            print(f"Output: {result.stdout}")
        else:
            print(f"Error processing {input_file}")
            print(f"Error Output: {result.stderr}")
    except Exception as e:
        print(f"Failed to run the process for {input_file}. Error: {str(e)}")


directory = '/home/mani/missing_reaction_processing/rxn_info_batchpro'
sub_directories = os.listdir(directory)

full_paths = []

for sub_dir in sub_directories:
    sub_dir_path = os.path.join(directory, sub_dir)
    full_paths.append(sub_dir_path)


combined_rxn_info = []

for full_path in full_paths:
    if 'reaction_info_with_ids.csv' in os.listdir(full_path):
        file_path = os.path.join(full_path, 'reaction_info_with_ids.csv')
        df = pd.read_csv(file_path)
        combined_rxn_info.append(df)
rxn_info_combined = pd.concat(combined_rxn_info, ignore_index=True)
# Save the combined DataFrame to a CSV file

results= rxn_info_combined.apply(asses_reaction_balance_row, axis= 1)
reaction_quality_assesment= pd.DataFrame(results.to_list()) 

rxn_quality_subset = rxn_quality_assesment[[
    'reactants', 'products', 'atom_difference',
    'mass_status', 'missing_atoms_side', 'missing_atoms_detail',
    'charge_difference', 'charge_status', 'reaction_status', 'rxn_id'
]]



combined_rxn_hyperlink_checked= pd.read_csv('/home/mani/ista_revision_final/Zenodo_Table_1.csv', low_memory= False)
combined_rxn_hyperlink_checked_unique = (
    combined_rxn_hyperlink_checked
    .groupby('rxn_text', as_index=False)
    .apply(aggregate_with_priority, include_groups=False)
    .reset_index(drop=True)
)

rxns_hyperlinked= combined_rxn_hyperlink_checked_unique[combined_rxn_hyperlink_checked_unique['hyperlink_status'] == 'Reaction CID is Hyperlinked'].reset_index(drop= True)
rxns_hyperlinked['rxn_id']= rxncid(rxns_hyperlinked)

rxns_hyperlinked['is_reactant_product_identical']= rxns_hyperlinked.apply(
    lambda row: set(row['cidsreactant'].split('|')) == set(row['cidsproduct'].split('|')),
    axis=1
    )

rxns_for_downstream_analysis= rxns_hyperlinked[rxns_hyperlinked['is_reactant_product_identical'] == False].reset_index(drop= True)
rxns_for_downstream_analysis['rxn_id'] = rxns_for_downstream_analysis['rxn_id'].astype(str).str.strip()


def _agg_join(series: pd.Series) -> str:
    values = [str(v).strip() for v in series.dropna().astype(str) if str(v).strip()]
    return ';'.join(sorted(set(values)))

_join_cols = ['enzyme', 'biosystem', 'source', 'file_name']
_agg_dict = {
    col: _agg_join
    for col in rxns_for_downstream_analysis.columns
    if col != 'rxn_id'
}

rxns_for_downstream_analysis_agg = (
    rxns_for_downstream_analysis
    .groupby('rxn_id', as_index=False)
    .agg(_agg_dict)
)

rxns_for_downstream_analysis_agg = rxns_for_downstream_analysis_agg.rename(columns={"rxn_id": "rxn_cid"})


# Select columns from aggregated dataframe
aggregated_df_subset = rxns_for_downstream_analysis_agg[['rxn_id', 'rxn_text', 'hyperlink_status', 'reactant_counts',
       'product_counts', 'reactant_cid_missing_info',
       'product_cid_missing_info', 'source', 'biosystem', 'file_name', 'is_reactant_product_identical']]


# Merge all three dataframes
merged_processed_rxns = (
    rxn_info_combined
    .merge(rxn_quality_subset, left_on='RXN_ID', right_on='rxn_id', how='left', suffixes=('', '_quality'))
    .merge(aggregated_df_subset, left_on='RXN_ID', right_on='rxn_id', how='inner', suffixes=('', '_proc'))
    .drop(columns=['rxn_id', 'rxn_id_quality', 'rxn_id_proc'], errors='ignore')
)


merged_processed_rxns.to_csv('/home/mani/missing_reaction_processing/rxn_info_combined.csv', index=False)
print("Combined reaction info saved to /home/mani/missing_reaction_processing/rxn_info_combined.csv")

