import logging
import os
import time
from typing import List, Dict, Any, Optional
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit import DataStructs 
from rxnmapper import RXNMapper
import pandas as pd
from tqdm.auto import tqdm
import pubchempy as pcp
from admet_ai import ADMETModel


log_file = "/home/mani/caffeine_trans_products/transformation_log_file.log"
os.makedirs(os.path.dirname(log_file), exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()  
    ]
)

# Set up logging to console
logger = logging.getLogger(__name__)

def get_matched_substructures_from_library(
    query_smiles: str, 
    sub_struc_lib_df: pd.DataFrame
) -> Optional[pd.DataFrame]:
    """
    Identifies substructures from a library that are present in a query molecule,
    and returns a subsetted DataFrame of the matched library entries.

    Args:
        query_smiles (str): The SMILES string of the molecule to query against the library.
        sub_struc_lib_df (pd.DataFrame): The DataFrame containing the substructure library,
                                         expected to have a 'Fragmants Smarts' column.

    Returns:
        Optional[pd.DataFrame]: A DataFrame containing only the rows from sub_struc_lib_df
                                that matched a pattern in the query_smiles,
                                with an additional 'query_smiles' column.
                                Returns None if the query SMILES is invalid,
                                the 'Fragmants Smarts' column is missing, or no matches are found.
    """
    
    if 'Fragmants Smarts' not in sub_struc_lib_df.columns:
        logger.error("Input DataFrame must contain a 'Fragmants Smarts' column.")
        return None
   
    query_molecule = Chem.MolFromSmiles(query_smiles)
    if query_molecule is None:
        logger.error(f"Invalid query SMILES string provided: '{query_smiles}'. Cannot proceed with matching.")
        return None

    #Prepare the list of SMARTS patterns from the library
    smarts_patterns: List[str] = sub_struc_lib_df['Fragmants Smarts'].apply(lambda x: str(x).strip('()')).tolist()
    num_patterns = len(smarts_patterns)

    # Generate the custom fingerprint for the query molecule
    # Initialize a bit vector with a length equal to the number of SMARTS patterns
    custom_fp = DataStructs.ExplicitBitVect(num_patterns)

    for i, smarts in enumerate(smarts_patterns):
        try:
            
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                logger.warning(f"Could not parse SMARTS pattern at index {i}: '{smarts}'. Skipping this pattern.")
                continue 
                        
            if query_molecule.HasSubstructMatch(pattern):                
                custom_fp.SetBit(i)
        except Exception as e:            
            logger.error(f"Error processing SMARTS pattern '{smarts}' at index {i}: {e}", exc_info=True)
            continue            
    
    fp_list = custom_fp.ToList()
    matched_indices = [i for i, bit_value in enumerate(fp_list) if bit_value == 1]

    if not matched_indices:
        logger.info(f"No patterns from the library matched '{query_smiles}'.")
        return None
    
    matched_substrucs_df = sub_struc_lib_df.iloc[matched_indices].copy()
    
    matched_substrucs_df.rename(columns={'Smiles': 'matched_smiles'}, inplace=True)
    
    matched_substrucs_df['query_smiles'] = query_smiles

    logger.info(f"Found {len(matched_indices)} matching patterns for '{query_smiles}'.")
    return matched_substrucs_df


rxn_mapper = None

def init_rxn_mapper_worker():
    """
    Initializer function for multiprocessing pool workers that need RXNMapper.
    This ensures RXNMapper is initialized only once per worker process.
    """
    global rxn_mapper
    try:
        rxn_mapper = RXNMapper()
        logger.info(f"RXNMapper initialized in worker process {os.getpid()}")
    except Exception as e:
        logger.error(f"Failed to initialize RXNMapper in worker process {os.getpid()}: {e}")
        rxn_mapper = None

def _get_mapped_reaction(constructed_reaction_smiles: str, rxn_id: Any) -> str | None:
    """
    Applies RXNMapper to get atom-mapped reaction SMILES.
    This function expects `rxn_mapper` to be globally available and initialized
    in the process where this function is called.

    Args:
        constructed_reaction_smiles (str): The reaction SMILES string to be mapped.
        rxn_id (Any): The reaction ID for logging purposes.

    Returns:
        str | None: The atom-mapped reaction SMILES, or None if mapping fails or rxn_mapper is not initialized.
    """
    global rxn_mapper
    if rxn_mapper is None:
        logger.warning(f"RXNMapper not initialized in this process. Skipping atom mapping for rxn_id: {rxn_id}.")
        return None
    try:
        mapping_result = rxn_mapper.get_attention_guided_atom_maps([constructed_reaction_smiles])[0]
        return mapping_result['mapped_rxn']
    except Exception as e:
        logger.error(f"RxnMapper error for rxn_id {rxn_id} in process {os.getpid()}: {e}")
        return None
    
def _rank_products(
    mapped_reaction: str,
    query_smiles: str,
    rxn_id: Any,
    ) -> List[Dict[str, Any]]:
    """
    Ranks products based on the overlap of atom maps with the query compound.

    Args:
        mapped_reaction (str): The atom-mapped reaction SMILES.
        query_smiles (str): The original SMILES of the query compound.
        unique_smiles_tuples (List[Tuple[str, ...]]): List of unique product SMILES tuples.
        rxn_id (Any): The reaction ID for logging purposes.

    Returns:
        List[Dict[str, Any]]: A list of dictionaries, each containing ranked product information.
                               Returns an empty list if ranking fails.
                               
    """
      
    ranked_results: List[Dict[str, Any]] = []
    try:
        mapped_parts = mapped_reaction.split(">>")
        mapped_reactants = mapped_parts[0].split(".")
        mapped_products = mapped_parts[1].split(".")

        query_cmpd_mapped = None
        for index, reactant_mapped_smiles in enumerate(mapped_reactants):
            reactant_mol_mapped = Chem.MolFromSmiles(reactant_mapped_smiles)
            if reactant_mol_mapped:
                for atom in reactant_mol_mapped.GetAtoms():
                    atom.SetAtomMapNum(0)
                reactant_canonical = Chem.MolToSmiles(reactant_mol_mapped, canonical=True, allHsExplicit=True)
                if reactant_canonical == Chem.MolToSmiles(Chem.MolFromSmiles(query_smiles), canonical=True, allHsExplicit=True):
                    query_cmpd_mapped = mapped_reactants[index]
                    break

        if not query_cmpd_mapped:
            logger.warning(f"Query compound not found in mapped reactants for rxn_id: {rxn_id}.")
            return []

        map_index = []
        mol = Chem.MolFromSmiles(query_cmpd_mapped)
        if mol:
            for atom in mol.GetAtoms():
                map_number = atom.GetAtomMapNum()
                if map_number != 0:
                    map_index.append(map_number)

        overlap_count = []
        for prod in mapped_products:
            
            mol = Chem.MolFromSmiles(prod)
            if mol:
                prod_map = [atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() != 0]
                matched_atoms = [i for i in prod_map if i in map_index]
                matched_count = len(matched_atoms)
                overlap_count.append((prod, matched_count))
            else:
                logger.warning(f"Could not parse product SMILES '{prod}' for ranking for rxn_id: {rxn_id}. Skipping.")

        overlap_count.sort(key=lambda x: x[1], reverse=True)

        current_rank = 1
        prev_count = None

        for idx, (product, count) in enumerate(overlap_count):


            if count != prev_count:
                rank_label = current_rank
                prev_count = count
            ranked_results.append({
                "rank": rank_label,
                "product_smiles": product,
                "no_matched_atom_maps": count
            })
            if idx + 1 < len(overlap_count):
                next_count = overlap_count[idx + 1][1]
                if next_count != count:
                    current_rank += 1

    except Exception as e:
        logger.error(f"Post-mapping processing error for rxn_id {rxn_id}: {e}")
        return []
    return ranked_results



def generate_transformation_product(trans_info: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Generates transformed products for a single reaction entry,
    without performing atom mapping or product ranking.
    This function is designed to be run in parallel.

    Args:
        trans_info (dict): A dictionary containing reaction details:
                            'query_smiles', 'matched_smiles', 'reac_temp', 'rxn_id',
                            'SANITIZED_MAPPED_REACTION'.

    Returns:
        List[Dict[str, Any]]: A list of dictionaries with transformed product information,
                              or an empty list if an error occurs.
    """
    results: List[Dict[str, Any]] = []

    try:
        query_smiles = trans_info['query_smiles']
        match_smiles = trans_info['matched_smiles']
        reac_template = trans_info['reac_temp']
        reaction_smiles = trans_info['SANITIZED_MAPPED_REACTION']
        rxn_id = trans_info['rxn_id']

        updated_reactants = [
            query_smiles if r == match_smiles else r
            for r in reaction_smiles.split('>>')[0].split('.')
        ]

        reactant_mols = [Chem.MolFromSmiles(r) for r in updated_reactants]
        if not all(reactant_mols):
            logger.warning(f"Invalid reactant SMILES found for rxn_id: {rxn_id}. Skipping reaction.")
            return []

        reaction = rdChemReactions.ReactionFromSmarts(reac_template)
        if not reaction:
            logger.error(f"Failed to parse reaction SMARTS for rxn_id: {rxn_id}.")
            return []

        products = reaction.RunReactants(reactant_mols)
        unique_smiles_tuples = list({
            tuple(Chem.MolToSmiles(mol, canonical=True) for mol in product_set if mol)
            for product_set in products
        })

        if not unique_smiles_tuples:
            logger.info(f"No unique products generated for rxn_id: {rxn_id}.")
            return []

        for product_tuple in unique_smiles_tuples:
            results.append({
                'rxn_id': rxn_id,
                'query_smiles': query_smiles,
                'matched_smiles': match_smiles,
                'sanitized_reaction': reaction_smiles,
                'reaction_template': reac_template,
                'transformed_products': product_tuple,
                'rxn_smiles_with_query_cmpd': f"{'.'.join(updated_reactants)}>>{'.'.join(product_tuple)}"
            })

    except Exception as e:
        logger.error(f"Initial transformation failed for rxn_id {trans_info.get('rxn_id', 'unknown')}: {e}")
        return []

    return results

def generate_transformation_products_parallel(merged_info: pd.DataFrame, max_workers: int = 16) -> pd.DataFrame:
    """
    Generates transformed products in parallel for a DataFrame of reaction information.

    Args:
        merged_info (pd.DataFrame): DataFrame containing reaction details for transformation.
        max_workers (int): Maximum number of worker processes.

    Returns:
        pd.DataFrame: DataFrame containing transformed product information.
    """
    all_results: List[Dict[str, Any]] = []
    reaction_rows: List[Dict[str, Any]] = merged_info.to_dict(orient="records")

    logger.info(f"Starting parallel product generation with {max_workers} workers...")
    with ThreadPoolExecutor(max_workers=16) as executor:
        futures = [executor.submit(generate_transformation_product, trans_info) for trans_info in reaction_rows]

        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing reactions"):
            try:
                result = future.result()
                if result:
                    all_results.append(result)
            except Exception as e:
                logger.error(f"Error retrieving result from worker: {e}", exc_info=True)

    all_results_flattened = [item for sublist in all_results for item in sublist]

    return pd.DataFrame(all_results_flattened)


def get_canonical_smiles(smiles: str):
    """
    Convert a SMILES string to its canonical form with all atom map numbers removed.

    Parameters:
        smiles (str): Input SMILES string.

    Returns:
        str: Canonical SMILES string without atom map numbers,
             or an error message if the input is invalid.
    """
    if smiles is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f"Error: Invalid SMILES string '{smiles}'."

    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, allHsExplicit=False)
    return canonical_smiles


def get_cid_from_smiles_with_retry(
    smiles: str,
    max_retries: int = 5,
    initial_delay: float = 1.0,
    cache: Optional[Dict[str, Optional[str]]] = None
) -> Optional[str]:
    """
    Retrieve PubChem CID for a SMILES string with retry logic and exponential backoff.
    
    Args:
        smiles (str): The SMILES string to look up.
        max_retries (int): Maximum number of retry attempts.
        initial_delay (float): Initial delay in seconds between retries.
    
    Returns:
        Optional[str]: The CID as a string, or None if the lookup fails.
    """
    if smiles is None:
        return None
    if cache is not None and smiles in cache:
        return cache[smiles]
    
    delay = initial_delay
    for attempt in range(max_retries):
        try:
            compounds = pcp.get_compounds(smiles, 'smiles')
            if compounds:
                cid = str(compounds[0].cid)
                if cache is not None:
                    cache[smiles] = cid
                return cid
            else:
                logger.warning(f"No compound found for SMILES: {smiles}")
                if cache is not None:
                    cache[smiles] = None
                return None
        except pcp.PubChemHTTPError as e:
            if 'ServerBusy' in str(e) or 'PUGREST.ServerBusy' in str(e):
                if attempt < max_retries - 1:
                    logger.warning(f"PubChem server busy for SMILES {smiles}. Retrying in {delay:.1f}s (attempt {attempt + 1}/{max_retries})...")
                    time.sleep(delay)
                    delay *= 2  # Exponential backoff
                else:
                    logger.error(f"PubChem server busy after {max_retries} attempts for SMILES: {smiles}")
                    if cache is not None:
                        cache[smiles] = None
                    return None
            else:
                logger.error(f"PubChem HTTP error for SMILES {smiles}: {e}")
                if cache is not None:
                    cache[smiles] = None
                return None
        except IndexError:
            logger.warning(f"No compounds found for SMILES: {smiles}")
            if cache is not None:
                cache[smiles] = None
            return None
        except Exception as e:
            logger.error(f"Unexpected error getting CID for SMILES {smiles}: {e}")
            if cache is not None:
                cache[smiles] = None
            return None
    
    return None


if __name__ == "__main__":

    try:
        sub_struc_lib = pd.read_csv('/home/mani/missing_reaction_processing/transformed_atom_substructure.csv')
        info_df = pd.read_csv('/home/mani/missing_reaction_processing/combined_info_with_templates.csv')
        
        pubchem_toc= pd.read_csv('/home/mani/PubChem-TOC-Pesticides-Insecticides-KEGG.csv')
        query_compounds = pubchem_toc['SMILES'].to_list()

        
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        exit()

    results_dir = "/home/mani/missing_reaction_processing/transformation_product_generation"
    os.makedirs(results_dir, exist_ok=True)
    
    rxn_mapper = RXNMapper()
    model = ADMETModel()
    cid_cache: Dict[str, Optional[str]] = {}

    for index, query_compound in enumerate(tqdm(query_compounds, desc="Query compounds")):

        logger.info(f"Processing query compound {index + 1}/{len(query_compounds)}: {query_compound}")

        matched_substruc_df = get_matched_substructures_from_library(query_compound, sub_struc_lib)
        if matched_substruc_df is None or matched_substruc_df.empty:
            logger.warning(f"No matched substructures for compound index {index + 1}")
            continue

        merged_info = pd.merge(
            matched_substruc_df,
            info_df[['RXN_ID', 'SANITIZED_MAPPED_REACTION', 'reac_temp']],
            left_on='rxn_id',
            right_on='RXN_ID',
            how='inner'
        )

        merged_info.drop(columns=['RXN_ID'], inplace=True)

        result = pd.DataFrame(generate_transformation_products_parallel(merged_info, max_workers=16))
        
        if result.empty:
            logger.warning(f"No transformation products generated for query compound {index + 1}. Skipping.")
            continue

        result['mapped_rxn_query_cmpd'] = result.apply(
            lambda row: _get_mapped_reaction(row['rxn_smiles_with_query_cmpd'], row['rxn_id']),
            axis=1
        )

        result['query_cmpd_to_prdt_correspondence'] = result.apply(
            lambda row: _rank_products(row['mapped_rxn_query_cmpd'], row['query_smiles'], row['rxn_id']),
            axis=1
        )

        rank_1_products = result['query_cmpd_to_prdt_correspondence'].apply(lambda x: x[0]['product_smiles'] if x and isinstance(x, list) and 'product_smiles' in x[0] else None)
        
        result['rank_1_product_canonical_smiles'] = rank_1_products.apply(lambda x: get_canonical_smiles(x))
 
        # Get unique SMILES for CID fetching (to avoid redundant API calls)
        unique_smiles = result['rank_1_product_canonical_smiles'].dropna().unique()
        
        # Retrieve CIDs with retry logic and rate limiting
        logger.info(f"Fetching CIDs for {len(unique_smiles)} unique SMILES...")
        cid_data = []
        for idx, smi in enumerate(unique_smiles):
            cid = get_cid_from_smiles_with_retry(smi, cache=cid_cache)
            cid_data.append({'smiles': smi, 'cid': cid})
            # Add a small delay between requests to avoid overwhelming the server
            if idx < len(unique_smiles) - 1:
                time.sleep(0.3)  # 300ms delay between requests
        
        smiles_to_cid_df = pd.DataFrame(cid_data)
             
        # Merge with your result DataFrame
        result = result.merge(
            smiles_to_cid_df,
            left_on='rank_1_product_canonical_smiles',
            right_on='smiles',
            how='left'
        ).drop(columns=['smiles'], axis=1)
        
        subdir = os.path.join(results_dir, f'query_cmpd_{index + 1}')
        os.makedirs(subdir, exist_ok=True)
        
        result.to_csv(os.path.join(subdir, f'query_compound_{index + 1}_transformation_products.csv'), index=False)

        # For ADME properties, use only unique SMILES from smiles_to_cid_df
        unique_product_smiles = smiles_to_cid_df['smiles'].dropna().tolist()
        adme_prop = []
        for smi in unique_product_smiles:
            try:
                pred = model.predict(smi)
                pred['smiles'] = smi
                adme_prop.append(pred)
            except Exception as e:
                logger.error(f"Error predicting ADME properties for SMILES {smi}: {e}")
        
        adme_prop_df = pd.DataFrame(adme_prop)
        adme_prop_df['Parent Compound/Product'] = 'Product'
        adme_prop_query = pd.DataFrame([{**model.predict(query_compound), 'smiles': query_compound, 'Parent Compound/Product': 'Parent Compound'}])
        
        adme_table = pd.concat([adme_prop_query, adme_prop_df], ignore_index=True)
     
        column_order = ['Parent Compound/Product', 'smiles']
        remaining_cols = [col for col in adme_table.columns if col not in column_order]
        final_order = column_order + remaining_cols
        adme_table = adme_table[final_order]
        
        adme_table.to_csv(
            os.path.join(subdir, f'adme_properties_query_compound_{index + 1}.csv'),
            index=False
        )





############### Thermodynamic stability ##################################
from thermo import Joback
import rdkit
from rdkit import Chem

def compute_gibbs_free_energy(smiles: str, temp_k: float = 298.15):
    """
    Calculate thermodynamic properties from SMILES using Joback group contribution method.
    
    Args:
        smiles: SMILES string of the molecule
        temp_k: Temperature in Kelvin (default: 298.15 K)
    
    Returns:
        Tuple of (enthalpy, entropy, gibbs_free_energy) in kJ/mol, kJ/(mol·K), kJ/mol
        Returns (None, None, None) if invalid SMILES
    """
    if smiles is None or str(smiles).strip() == '':
        return 'Invalid SMILES'    
   
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 'Failed To Generate Rdkit Mol Object From Smiles'
    
    try:
        data = Joback(mol)
        result= data.estimate()
        # 1. Get Enthalpy (Hf) in kJ/mol
        enthalpy = result['Hf'] / 1000  # Convert J to kJ
        # 2. Get Gibbs free energy (Gf) in kJ/mol  
        gibbs_free_energy = result['Gf'] / 1000  # Convert J to kJ
        # 3. Derive Entropy (S) from relationship: G = H - TS
        # Rearranging: S = (H - G) / T
        entropy = (enthalpy - gibbs_free_energy) / temp_k  # kJ/(mol·K)

        # 4. Determine Stability Rating
        if gibbs_free_energy < -100:
            stability = "Highly Stable"
        elif gibbs_free_energy < 0:
            stability = "Stable"
        else:
            stability = "Unstable / Highly Reactive"

        return {
            'Enthalpy': round(enthalpy, 2), 
            'Entropy': round(entropy, 4), 
            'Gibbs Free Energy': round(gibbs_free_energy, 2),
            'Stability': stability
        }
        
    except Exception as e:
                
        return {
            'Enthalpy': None, 
            'Entropy': None, 
            'Gibbs Free Energy': None, 
            'Stability': f'Error: {str(e)}'
        }


target_directory = "/home/mani/missing_reaction_processing/transformation_product_generation"
table2= pd.read_csv('/home/mani/ista_revision_final/Zenodo_Table_2.csv', low_memory= False)

full_file_paths = [
    os.path.join(root, name)
    for root, _, files in os.walk(target_directory)
    for name in files
    if "transformation" in name
]

df_list = [pd.read_csv(p) for p in full_file_paths]
final_trans_products = pd.concat(df_list, ignore_index=True)

insec_trans_dict = {smiles: group for smiles, group in final_trans_products.groupby('query_smiles')}

gibbs_free_energy_calculation_results = []

for key in insec_trans_dict.keys():
    df = insec_trans_dict.get(key)
    gibbs_results = []
    
    # Get unique SMILES and compute Gibbs free energy for each
    for smil in set(df['rank_1_product_canonical_smiles']):
        try:
            # Compute Gibbs free energy
            energy_result = compute_gibbs_free_energy(smil)
            
            # Handle different return types from compute_gibbs_free_energy
            if isinstance(energy_result, dict):
                gibbs_results.append({'smiles': smil, **energy_result})
            elif isinstance(energy_result, (int, float)):
                # If function returns just a number, store it
                gibbs_results.append({'smiles': smil, 'Gibbs_Free_Energy': energy_result})
            else:
                # If unknown type, convert to string and log
                print(f"Warning: Unexpected return type for {smil}: {type(energy_result)}")
                gibbs_results.append({'smiles': smil, 'error': str(energy_result)})
        except Exception as e:
            # Handle errors gracefully without breaking the loop
            print(f"Error computing Gibbs free energy for {smil}: {str(e)}")
            gibbs_results.append({'smiles': smil, 'error': str(e)})
    
    # Create DataFrame from results
    if gibbs_results:
        gibbs_df = pd.DataFrame(gibbs_results)
        
        # Rank by Gibbs Free Energy if column exists
        if 'Gibbs_Free_Energy' in gibbs_df.columns:
            gibbs_df['rank'] = gibbs_df['Gibbs_Free_Energy'].rank(ascending=True, method='min')
        elif 'Gibbs Free Energy' in gibbs_df.columns:
            gibbs_df['rank'] = gibbs_df['Gibbs Free Energy'].rank(ascending=True, method='min')
        
        gibbs_df['query_smiles'] = key  # Add query SMILES for reference
        gibbs_free_energy_calculation_results.append(gibbs_df)
    else:
        print(f"No results for query SMILES: {key}")


gibbs_free_energy_results = pd.concat(gibbs_free_energy_calculation_results).drop(columns=['Stability'], errors='ignore')

merged_gibbs_trans = final_trans_products.merge(
    gibbs_free_energy_results,
    left_on='rank_1_product_canonical_smiles',
    right_on='smiles',
    how='left'    
).drop(columns= ['smiles', 'query_smiles_y', 'error'])

# Map REACTION_STATUS and TAX_ID from final_combined_rxns
merged_gibbs_trans = merged_gibbs_trans.merge(
    table2[['RXN_ID', 'ECS', 'TAX_ID', 'REACTION_STATUS']],
    left_on='rxn_id',
    right_on='RXN_ID',
    how='left'
).drop(columns=['RXN_ID'])

merged_gibbs_trans.to_csv('/home/mani/ista_revision_final/Zenodo_Table5.csv', index= False)


######## Zenodo Table 6 ###########################################

target_directory = "/home/mani/missing_reaction_processing/transformation_product_generation"

full_file_paths = [
    os.path.join(root, name)
    for root, _, files in os.walk(target_directory)
    for name in files
    if "adme_properties" in name
]

df_list = [pd.read_csv(p) for p in full_file_paths]
combined_adme= pd.concat(df_list).reset_index(drop= True)

parent_cmpd_index = list(combined_adme[combined_adme['Parent Compound/Product'] == 'Parent Compound'].index)

splitted_df = []
for index, value in enumerate(parent_cmpd_index):
  start_index = value
  if index < len(parent_cmpd_index) - 1:
    end_index = parent_cmpd_index[index + 1]
  else:
    end_index = len(combined_adme)
  df_subset = combined_adme.iloc[start_index:end_index].copy()
  splitted_df.append(df_subset)

top_rank_computed_dfs= []
for df in splitted_df:
    percentiles = df.filter(like='_drugbank_approved_percentile')
    top_1_mask = percentiles >= 95
    df['top_1_percentile_count'] = top_1_mask.sum(axis=1)
    top_rank_computed_dfs.append(df)


final_adme_table= pd.concat(top_rank_computed_dfs, ignore_index=True)
final_adme_table.to_csv('/home/mani/ista_revision_final/Zenodo_Table6.csv', index= False)


############# SI Table 3 Generation ###################################
import re
import pandas as pd

table2= pd.read_csv('/home/mani/ista_revision_final/Zenodo_Table_2.csv', low_memory= False)
trans_prducts= pd.read_csv('/home/mani/ista_revision_final/Zenodo_Table5.csv', low_memory= False)
adme_table= pd.read_csv('/home/mani/ista_revision_final/Zenodo_Table6.csv', low_memory= False)
pubchem_toc = pd.read_csv('/home/mani/PubChem-TOC-Pesticides-Insecticides-KEGG.csv')

query_cids_unique= list(set(adme_table['query_cid']))

splitted_df= []
summary_rows = []

for cid in query_cids_unique:
    df= adme_table[adme_table['query_cid'] == cid]
    splitted_df.append(df)
    
for df in splitted_df:
    query_cid= list(set(df['query_cid']))[0]
    target_cids = []
    for _cid in df['cid'].dropna():
        try:
            target_cids.append(int(float(_cid)))
        except (ValueError, TypeError):
            continue
    target_cids = list(set(target_cids))
    no_of_trabsformation_products= len(list(set(df['smiles'])))
    
    target_compound_rxn_id_list= []
    cids_associted_with_any_reaction= []
    cids_not_associated_with_any_reaction= []
    for tar_cid in target_cids:
        if pd.isna(tar_cid):
            continue
        target_cid_str = str(tar_cid).strip()
        if not target_cid_str:
            continue
        token_pattern = rf'(^|\||>>)\s*{re.escape(target_cid_str)}\s*($|\||>>)'
        table2_exact = table2[table2['RXN_ID'].astype(str).str.contains(token_pattern, na=False, regex=True)]
        reaction_id_list = table2_exact['RXN_ID'].tolist()
        if reaction_id_list:
            cids_associted_with_any_reaction.append(tar_cid)
        else:
            cids_not_associated_with_any_reaction.append(tar_cid)
        target_compound_rxn_id_list.extend(reaction_id_list)

    target_compound_rxn_id_list = sorted(set(target_compound_rxn_id_list))

    query_compound_rxn_id_list = []
   
    query_cid_str = str(query_cid).strip()
    
    if not query_cid_str:
        continue
    token_pattern = rf'(^|\||>>)\s*{re.escape(query_cid_str)}\s*($|\||>>)'
    table2_exact = table2[table2['RXN_ID'].astype(str).str.contains(token_pattern, na=False, regex=True)]
    reaction_id_list = table2_exact['RXN_ID'].tolist()
    query_compound_rxn_id_list.extend(reaction_id_list)

    query_compound_rxn_id_list = sorted(set(query_compound_rxn_id_list))

    overlapping_reaction_id = sorted(set(query_compound_rxn_id_list).intersection(target_compound_rxn_id_list))

    compound_match = pubchem_toc[pubchem_toc['Compound CID'] == query_cid]
    compound_name = compound_match['Name'].iloc[0] if not compound_match.empty else None
    compound_class = compound_match['Compound Class'].iloc[0] if not compound_match.empty else None


    summary_rows.append({
        'Compound Name': compound_name,
        'Query CIDs': query_cid,
        'No of Transformation Product Metabolites Enlisted in Pubchem': len(target_cids),
        'Transformation Product CIDs': target_cids,
        'CIDs Associated with Any Reaction PubChem': cids_associted_with_any_reaction,
        'CIDs Not Associated With Any Reaction PubChem': cids_not_associated_with_any_reaction,
        'Overlapping Reaction ID Between Query CID and Transformation Product': overlapping_reaction_id,
        'Overlapping Reaction Count Between Query CID and Transformation Product': len(overlapping_reaction_id),
        'Compound Class': compound_class,
    })

summary_df = pd.DataFrame(summary_rows)

summary_df['No of Unique Transformation Product Not Enlisted in PubChem'] = (
    summary_df['No of Unique Transformation Products']
    - summary_df['No of Transformation Product Metabolites Enlisted in Pubchem']
)

summary_df.to_csv('/home/mani/ista_revision_final/tables3.csv', index= False)

############# SI Table 4 Generation ###################################

adme_table_2 = final_adme_table.loc[:, ~adme_table.columns.str.contains('_drugbank_approved_percentile', na=False)]
adme_table_3 = adme_table_2[adme_table_2['Parent Compound/Product'] != 'Parent Compound'].reset_index(drop=True)
adme_4= adme_table_3.loc[:, 'molecular_weight': 'VDss_Lombardo']
adme_4_stats = adme_4.describe()
adme_4_stats.to_csv('/home/mani/ista_revision_final/tableS4.csv', index= False)
        