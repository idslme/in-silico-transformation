import logging
import os
from typing import List, Dict, Any, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit import DataStructs 
from rxnmapper import RXNMapper
import pandas as pd
from tqdm.auto import tqdm
import pubchempy as pcp
from admet_ai import ADMETModel


log_file = "/home/mani/pubchem_rxns_2025_processing/transformation_products_generation/transformation_log_file.log"
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
    with ProcessPoolExecutor(max_workers=16) as executor:
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

if __name__ == "__main__":

    try:
        sub_struc_lib = pd.read_csv('/home/mani/pubchem_rxns_2025_processing/transformed_atom_substructure.csv')
        info_df = pd.read_csv('/home/mani/pubchem_rxns_2025_processing/combined_info_28may_with_templates.csv')
        insecticides = pd.read_csv('/home/mani/PubChem-TOC-Pesticides-Insecticides-KEGG.csv')
        query_compounds = insecticides['SMILES'].dropna().tolist()
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        exit()

    results_dir = "/home/mani/pubchem_rxns_2025_processing/transformation_products_generation"
    os.makedirs(results_dir, exist_ok=True)
    
    rxn_mapper = RXNMapper()
    model = ADMETModel()

    for index, query_compound in enumerate(query_compounds):

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
 
        smiles_series = pd.Series(list(set(result['rank_1_product_canonical_smiles'])), name='smiles')
        smiles_to_cid_df = smiles_series.apply(lambda smi: pd.Series({'smiles':smi, 'cid':str(pcp.get_compounds(smi, 'smiles')[0].cid)}))
             
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

        adme_prop = smiles_series.apply(lambda smi: {**model.predict(smi), 'smiles': smi})
        adme_prop_df = pd.DataFrame(adme_prop.tolist())
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
        