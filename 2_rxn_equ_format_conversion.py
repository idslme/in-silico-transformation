"""Converting rxn_cid to different rxn equations format"""
import os
from concurrent.futures import  ThreadPoolExecutor, as_completed
import logging
import requests
import gzip
import pandas as pd
import math

logging.basicConfig(
    filename='missing_compounds.log',
    level=logging.WARNING,
    format='%(asctime)s - %(message)s'
)


def log_missing_cid(cid: int):
    logging.warning(f"Missing CID: {cid}. No information available.")

def process_cids(cid_list: str, compound_lookup: dict):
    """Get SMILES, formula and name for a given list of CIDs"""    
    cmpd_smiles, cmpd_formula, cmpd_name = [], [], []
    incomplete_flag = False
    
    cids = [int(cid.strip()) for cid in cid_list.split('|')] if '|' in cid_list else[int(cid_list.strip())]
    
    for cid in cids:
        cmpd_info = compound_lookup.get(cid)
        if cmpd_info:
            cmpd_smiles.append(cmpd_info.get('smiles', 'Missing_Info'))
            cmpd_formula.append(cmpd_info.get('mf', 'Missing_Info'))
            cmpd_name.append(cmpd_info.get('cmpdname', 'Missing_Info'))
        else:
            log_missing_cid(cid)
            cmpd_smiles.append('Missing_Info')
            cmpd_formula.append('Mising_Info')
            cmpd_name.append('Missing_Info')
            incomplete_flag = True
    return cmpd_smiles, cmpd_formula, cmpd_name, incomplete_flag

def convert_rxn_equ_format(rxn_cid: str, compound_lookup: dict):
    """Convert rxn_cid to different reaction equation formats"""
    if '>>' not in rxn_cid:
        print(f'Invalid Rxn Format for Reaction: {rxn_cid}')
        return None

    reactants, products = rxn_cid.split('>>')
    react_smiles, react_formula, react_names, react_incomplte = process_cids(reactants, compound_lookup)
    prdt_smiles, prdt_formula, prdt_names, prdt_incomplete = process_cids(products, compound_lookup)
    
    incomplete_flag = react_incomplte or prdt_incomplete
    
    if not incomplete_flag:
        reaction_smiles = str(".".join(react_smiles) + ">>" + ".".join(prdt_smiles))
        reaction_formula = str(" + ".join(react_formula) + ">>" + " + ".join(prdt_formula))
        reaction_name = str(" + ".join(react_names) +  ">>" + " + ".join(prdt_names))
    else:
        reaction_smiles = reaction_formula = reaction_name = "Incomplete"
        
    return {
        "reaction_cid": rxn_cid,
        "reaction_smiles": reaction_smiles,
        "reaction_formula": reaction_formula,
        "reaction_name": reaction_name,
        "incomplete_flag": incomplete_flag
    }
    
def convert_rxn_equ_batchmode(input_data: pd.DataFrame, compound_lookup: dict, batch_size: int):
    """Batch processing of conversion of diffrent reaction equation formats"""
    batch_results = []
    
    for i in range(0, len(input_data['rxn_cid']), batch_size):
        batch = input_data['rxn_cid'][i:i+batch_size]
        print(f"Processing batch {i // batch_size + 1} of size {len(batch)}")
        with ThreadPoolExecutor(max_workers=16) as executor:
            futures = {executor.submit(convert_rxn_equ_format, rxn, compound_lookup): (rxn, idx)
                       for idx, rxn in enumerate(batch)}

            for future in as_completed(futures):
                result = future.result()
                if result:
                    batch_results.append(result)
            futures.clear()
        
    rxn_with_diff_format = pd.DataFrame(batch_results)
    return rxn_with_diff_format    

# Download the PubChem TOC file
url = 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22compound%22,%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_compound_cache_k9Q1vLg73YfqqVWw18gcmZo6zFr4EFvzIdZAvzrHUr463m4%22,%22where%22:{%22ands%22:[{%22input%22:{%22type%22:%22netcachekey%22,%22idtype%22:%22cid%22,%22key%22:%22k9Q1vLg73YfqqVWw18gcmZo6zFr4EFvzIdZAvzrHUr463m4%22}}]}}'
response = requests.get(url)

if response.status_code == 200:
    compressed_file_path = '/home/mani/pubchem-reactions/pubchem_toc.csv.gz'
    with gzip.open(compressed_file_path, 'wb') as f_out:
        f_out.write(response.content)
    print(f"File successfully downloaded and saved as {compressed_file_path}")
else:
    print(f"Failed to download the file. Status code: {response.status_code}")

# read pharmacology and biochemistry toc file downloaded from PubChem
toc = pd.read_csv('/home/mani/pubchem-reactions/pubchem_toc.csv.gz', compression='gzip', low_memory = False)
compound_lookup = toc.set_index('cid').to_dict(orient='index')
unprocessed_overall = pd.read_csv('/home/mani/pubchem_rxns_2025_processing/unprocessed_overall.csv', low_memory = False)
reaction_id = pd.DataFrame(unprocessed_overall['rxn_cid'], columns=['rxn_cid'])

rxn_with_diff_format = convert_rxn_equ_batchmode(reaction_id, compound_lookup, 5000)
rxn_diff_euq_complete = rxn_with_diff_format[rxn_with_diff_format['incomplete_flag'] == False]
rxn_diff_euq_complete.to_csv('/home/mani/pubchem_rxns_2025_processing/rxn_diff_euq_complete_unprocessed.csv', index=False)

# Get the CIDs for which no information is available
cid_no_available_info = []

with open('missing_compounds.log', 'r') as file:
    for line in file:
        parts = line.split(' ')
        cid = [parts[index + 1].rstrip('.') for index, part in enumerate(parts) if 'CID:' in part]
        cid_no_available_info.extend(cid)
        
if cid_no_available_info:
    with open('/home/mani/pubchem_rxns_2025_processing/missing_cids.txt', 'w') as output_file:
        for cid in cid_no_available_info:
            output_file.write(f"{cid}\n")
    print(f"Written {len(cid_no_available_info)} missing CIDs to 'missing_cids.txt'.")
else:
    print("No missing CIDs found.")   

