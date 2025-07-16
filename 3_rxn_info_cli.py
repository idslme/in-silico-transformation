"""
Generating Reaction Info from Reaction SMILES using rxn_insight workflow
"""
import os
import logging
import argparse
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from rxn_insight.classification import ReactionClassifier
from rxn_insight.reaction import Reaction, Molecule
from tqdm import tqdm

logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

def get_fg(rxn: str):
    """Get Functional Group for Each SMILES in Reaction SMILES using rxn_insight Molecule class"""
    try:
        reactants, products = rxn.split(">>")
    except ValueError:
        raise ValueError(f'Invalid Reaction Equation Format for:{rxn}')

    def get_functional_groups_dict(smiles_list):
        """Function to get a dictionary of SMILES and their functional groups"""
        return {smiles: Molecule(smiles).get_functional_groups() for smiles in smiles_list}

    reactant_dict = get_functional_groups_dict(reactants.split("."))
    product_dict = get_functional_groups_dict(products.split("."))

    return reactant_dict, product_dict

def get_rxn_info(rxn_smiles: str, rxn_id: str):
    """Generate Reaction Info For a Single Reaction Equation"""

    if '-->' in rxn_smiles:
        rxn_smiles = rxn_smiles.replace('-->', '>>')
    elif '>>' not in rxn_smiles:
        raise ValueError(f'Reaction Equation Format is Invalid for {rxn_smiles}')
    try:
        rxn = Reaction(rxn_smiles)
        info_rxn = rxn.get_reaction_info()
        info_rxn.pop('TAG', None)
       
        classifier = ReactionClassifier(rxn_smiles)
        info_rxn['NOS_REACTION_CENTER'] = classifier.nos_reaction_center
        info_rxn['RING_CHANGING'] = classifier.ring_changing()
        info_rxn['IS_AROMATIC_HETEROCYCLE'] = classifier.is_aromatic_heterocycle()
        info_rxn['IS_ACYLATION'] = classifier.is_acylation()
        info_rxn['IS_CC_COUPLING'] = classifier.is_cc_coupling()
        info_rxn['IS_DEPROTECTION'] = classifier.is_deprotection()
        info_rxn['IS_FGA'] = classifier.is_fga()
        info_rxn['IS_FGI'] = classifier.is_fgi()
        info_rxn['IS_HEREROATOM_ALKYLATION'] = classifier.is_heteroatom_alkylation()
        info_rxn['IS_OXIDATION'] = classifier.is_oxidation()
        info_rxn['IS_PROTECTION'] = classifier.is_protection()
        info_rxn['IS_REDUCTION'] = classifier.is_reduction()
        info_rxn['SANITIZED_TRANSFORMATION_MAPPING'] = classifier.sanitized_transformation_mapping
        info_rxn['REACTION_CENTER_ATOMS'] = classifier.reaction_center_atoms
        info_rxn['REACTION_CENTER_IDX'] = classifier.reaction_center_idx
       
        react_fg_dict, prod_fg_dict =  get_fg(rxn_smiles)
       
        info_rxn['FG_REACTANTS'] = react_fg_dict
        info_rxn['FG_PRODUCTS'] = prod_fg_dict
       
        order_of_dictionary = ['REACTION', 'MAPPED_REACTION', 'NAME', 'CLASS',
 'N_REACTANTS', 'N_PRODUCTS', 'FG_REACTANTS', 'FG_PRODUCTS', 'BY-PRODUCTS',
 'PARTICIPATING_RINGS_REACTANTS', 'PARTICIPATING_RINGS_PRODUCTS',
 'ALL_RINGS_PRODUCTS', 'SOLVENT', 'REAGENT', 'CATALYST',
 'IS_CC_COUPLING', 'IS_AROMATIC_HETEROCYCLE', 'IS_ACYLATION',
 'IS_DEPROTECTION', 'IS_FGA', 'IS_FGI', 'IS_OXIDATION',
 'IS_PROTECTION', 'IS_REDUCTION', 'SCAFFOLD',
 'NOS_REACTION_CENTER', 'RING_CHANGING', 'SANITIZED_TRANSFORMATION_MAPPING','REF']

        rxn_info = {key: info_rxn[key] for key in order_of_dictionary}
         
        return rxn_info
         
    except Exception as e:
        error_message = f"Error processing reaction {rxn_id}: {rxn_smiles}. Error: {str(e)}"
        logging.error(error_message)
        return None
       
def get_rxn_info_batch(df: pd.DataFrame, batch_size: int, max_num_cpus: int, output_dir: str, min_num_cpus: int=8, error_log_path: str="error_log.txt"):
    """Generate Reaction Info From a Dataframe Containing Reaction SMILES and Reaction Id in Batch Mode"""

    workers = min(max_num_cpus, max(min_num_cpus, os.cpu_count()))
    combined_rxn_info = []
    num_batches = len(df) // batch_size + (1 if len(df) % batch_size != 0 else 0)

    batch_error_log = []
   
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = []
        for batch_num in range(num_batches):
            start_index = batch_num * batch_size
            end_index = min((batch_num + 1) * batch_size, len(df))
            batch = df.iloc[start_index:end_index]
            logging.info(f'Processing Batch {batch_num + 1}/{num_batches}')
           
            for _, row in batch.iterrows():
                futures.append(executor.submit(get_rxn_info, row['rxn_smiles'], row['rxn_id']))
       
        for future in tqdm(as_completed(futures), total=len(futures)):
            rxn_info = future.result()
            if rxn_info:
                combined_rxn_info.append(rxn_info)
            else:
                row = batch.iloc[futures.index(future)]
                batch_error_log.append(f"Error processing reaction {row['rxn_id']}: {row['rxn_smiles']}")

        if batch_error_log:
            error_log_file = os.path.join(output_dir, "batch_error_log.txt")
            with open(error_log_file, 'w') as error_log:
                error_log.write("\n".join(batch_error_log))
            logging.info(f"Errors encountered during batch processing. See {error_log_file} for details.")

        info_df = pd.DataFrame(combined_rxn_info)
        info_df.to_csv(os.path.join(output_dir, "reaction_info.csv"), index=False)
        logging.info(f'The Result File is Written in the {output_dir}')

if __name__ == '__main__':
   
    parser = argparse.ArgumentParser(description='Reaction Info Generation From Reaction Smiles')
    parser.add_argument('input_file', type=str, help='CSV file with Reaction SMILES and its corresponding Reaction ID')
    parser.add_argument('batch_size', type=int, help='Batch Size for Parallel Processing')
    parser.add_argument('max_num_cpus', type=int, help='Maximum Number of Workers for Parallel Processing', default=8)
    parser.add_argument('output_dir', type=str, help='Output Directory Path')
    args = parser.parse_args()
    
    required_columns = ['rxn_smiles', 'rxn_id']
    input_file = pd.read_csv(args.input_file)
    for col in required_columns:
        if col not in input_file.columns:
            raise ValueError(f"Column {col} not found in the input file.")
    get_rxn_info_batch(df=input_file, batch_size=args.batch_size, max_num_cpus=args.max_num_cpus, output_dir=args.output_dir)

