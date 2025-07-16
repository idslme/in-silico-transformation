"""Transformation template generation for reactions using rdchiral."""
import logging
import os
import pandas as pd
from rdchiral import template_extractor
from rdkit import Chem
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

def configure_logging(input_file_path: str):
    """Configure logging to save log file in the same directory as the input file."""
    log_dir = os.path.dirname(input_file_path)  
    log_file_path = os.path.join(log_dir, 'reac_temp_generation.log')    

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    logging.basicConfig(level=logging.INFO,  
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=[
                            logging.StreamHandler(), 
                            logging.FileHandler(log_file_path, mode='w')  
                        ])
    
    
def mapped_rxn_to_rdchiral_input(rxn: str, rxn_id: str):
    """Converts a reaction string into the format required by rdchiral."""
    try:
        # Split reaction into reactants and products
        reaction = {
            'reactants': rxn.split('>>')[1],
            'products': rxn.split('>>')[0],
            '_id': rxn_id
        }
        
        # Convert reactants and products into RDKit molecules
        reactant_mols = [Chem.MolFromSmiles(reactant) for reactant in reaction['reactants'].split('.')]
        product_mols = [Chem.MolFromSmiles(product) for product in reaction['products'].split('.')]
                
        if None in reactant_mols or None in product_mols:
            raise ValueError(f"Error in converting SMILES to molecules for reaction {rxn_id}.")
        
        logging.info(f"Successfully processed reaction {rxn_id}.")
        return reaction, reactant_mols, product_mols
    except Exception as e:
        logging.error(f"Error processing reaction {rxn_id}: {e}")
        return None, [], []
    
    
# Extract the transformation template for a given reaction
def get_transformation_template(rxn: str, rxn_id: str):
    """Get the transformation template for a given reaction."""
    try:        
        reaction, reactant_mols, product_mols = mapped_rxn_to_rdchiral_input(rxn, rxn_id)        
        temp = template_extractor.extract_from_reaction(reaction)               
        return temp
    except Exception as e:
        logging.error(f"Error getting transformation template for reaction {rxn_id}: {e}")
        return None

def get_transformation_template_parallel(df, rxn_col='SANITIZED_MAPPED_REACTION', id_col='RXN_ID', n_workers=None):
    """
    Apply get_transformation_template in parallel to a DataFrame of reactions, with progress bar.
    """
    if n_workers is None:
        n_workers = max(1, cpu_count() // 2)
    args = [(row[rxn_col], row[id_col]) for _, row in df.iterrows()]
    results = [None] * len(args)
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_idx = {executor.submit(get_transformation_template, *arg): idx for idx, arg in enumerate(args)}
        for i, future in enumerate(tqdm(as_completed(future_to_idx), total=len(future_to_idx), desc="Processing reactions")):
            idx = future_to_idx[future]
            try:
                results[idx] = future.result()
            except Exception as e:
                results[idx] = None
    return results



if __name__ == "__main__":
    
    input_file_path = '/home/mani/pubchem_rxns_2025_processing/combined_info_15may_with_templates.csv'
    configure_logging(input_file_path) 
    
    info_df = pd.read_csv('/home/mani/pubchem_rxns_2025_processing/combined_info_15may_with_templates.csv')
    reaction_template = get_transformation_template_parallel(info_df, rxn_col='SANITIZED_MAPPED_REACTION', id_col='RXN_ID', n_workers=16)
    
    rxn_temp_nan_removed = [mem for mem in reaction_template if isinstance(mem, dict)]
    reaction_template_forward = pd.DataFrame(rxn_temp_nan_removed)
    reaction_template_forward.to_csv('/home/mani/pubchem_rxns_2025_processing/reaction_template_forward.csv', index=False)
    
    info_df.drop(columns=['reac_temp'], inplace=True, errors='ignore')

    filtered_df = info_df[info_df['RXN_ID'].isin(reaction_template_forward['reaction_id'])]  
    
       
    merged_df = pd.merge(
        filtered_df,
        reaction_template_forward[['reaction_id', 'reaction_smarts']],
        left_on='RXN_ID',
        right_on='reaction_id',
        how='left'
    )
    
    merged_df.rename(columns={'reaction_smarts': 'reac_temp'}, inplace=True)
    merged_df.drop(columns=['reaction_id'], inplace=True, axis=1)
    merged_df.to_csv('/home/mani/pubchem_rxns_2025_processing/combined_info_28may_with_templates.csv', index=False)

