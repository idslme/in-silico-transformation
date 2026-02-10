"""Converting rxn_cid to different rxn equations format"""
import time
import re
from concurrent.futures import  ThreadPoolExecutor, as_completed
import logging
import requests
from multiprocessing import Pool
import gzip
import html
import pandas as pd

logging.basicConfig(
    filename='missing_compounds.log',
    level=logging.WARNING,
    format='%(asctime)s - %(message)s'
)


def log_missing_cid(cid: int):
    logging.warning(f"Missing CID: {cid}. No information available.")

# Fromatting HTML format of Reaction Equation in the text format

# 1) Helpers to clean HTML and arrows
HTML_TAG = re.compile(r'<[^>]+>')

ARROW_FIXES = {
    'âŸ¶': '⟶',   # mojibake for U+27F6 (Long Rightwards Arrow)
    'â†’': '→',   # mojibake for U+2192 (Rightwards Arrow)
    '->': '⟶',    # normalize ASCII arrow to unicode
    '→': '⟶',
}

def normalize_arrows(s: str) -> str:
    for bad, good in ARROW_FIXES.items():
        s = s.replace(bad, good)
    return s

def strip_html_formatting(s: str) -> str:
    """
    - HTML-unescape entities (&lt;sup&gt;+&lt;/sup&gt; -> <sup>+</sup>)
    - Strip tags (<i>R</i> -> R, <sup>+</sup> -> +)
    - Normalize arrows (âŸ¶, ->, → -> ⟶)
    - Condense whitespace
    """
    if s is None:
        return ''
    s = html.unescape(str(s))              # unescape &lt; &gt; &amp; etc.
    s = HTML_TAG.sub('', s)                # remove HTML tags
    s = normalize_arrows(s)                # fix mojibake/variants
    s = re.sub(r'\s+', ' ', s).strip()     # tidy spaces
    return s


## Aggregate the dataframe based on the rxn_string
def aggregate_with_priority(group):
    # Find hyperlinked row if exists, otherwise use first row
    hyperlinked = group[group['hyperlink_status'] == 'Reaction CID is Hyperlinked']
    priority_row = hyperlinked.iloc[0] if len(hyperlinked) > 0 else group.iloc[0]
    definition_clean_value = (
        priority_row['definition_clean']
        if 'definition_clean' in priority_row.index
        else strip_html_formatting(priority_row.get('definition', ''))
    )
    
    return pd.Series({
        'cids': priority_row['cids'],
        'cidsreactant': priority_row['cidsreactant'],
        'cidsproduct': priority_row['cidsproduct'],
        'definition': priority_row['definition'],
        'definition_clean': definition_clean_value,
        'hyperlink_status': priority_row['hyperlink_status'],
        'reactant_counts': priority_row['reactant_counts'],
        'product_counts': priority_row['product_counts'],
        'reactant_cid_missing_info': priority_row['reactant_cid_missing_info'],
        'product_cid_missing_info': priority_row['product_cid_missing_info'],
        # Aggregate these columns from all rows
        'enzyme': ';'.join([str(x) for x in group['enzyme'].dropna().unique() if str(x).strip()]),
        'biosystem': ';'.join([str(x) for x in group['biosystem'].dropna().unique() if str(x).strip()]),
        'source': ';'.join([str(x) for x in group['source'].dropna().unique() if str(x).strip()]),
        'file_name': ';'.join([str(x) for x in group['file_name'].dropna().unique() if str(x).strip()])
    })


# Function to convert reactant and product cids to reaction equation id
def react_prod_cid_to_rxncid(reactant_cids: str, product_cids: str):
    """Function to convert reactant and product cids to reaction equation id"""
    reactant_list = sorted(reactant_cids.split('|')) if '|' in reactant_cids else[reactant_cids]
    product_list = sorted(product_cids.split('|')) if '|' in product_cids else [product_cids]
    reactant_product = '|'.join(reactant_list) + ">>" + '|'.join(product_list)
    return reactant_product

def rxncid(df:pd.DataFrame):
    """Function to generate reaction equation id for each reaction with parallel processing"""
    
    rows = [(row['cidsreactant'], row['cidsproduct']) for index, row in df.iterrows()]
    with Pool() as pool:
        result = pool.starmap(react_prod_cid_to_rxncid, rows)
    return result


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


def get_smiles_from_cid(cid, timeout=10, max_retries=3, sleep_time=0.2):
    """
    Get SMILES from PubChem CID with retry logic.
    
    Args:
        cid: PubChem Compound ID
        timeout: Request timeout in seconds
        max_retries: Maximum number of retry attempts
        sleep_time: Base sleep time between requests (seconds)
    
    Returns:
        SMILES string or None if failed
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/TXT"
    
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=timeout)
            if response.status_code == 200:
                time.sleep(sleep_time)  # Rate limiting
                return response.text.strip()
            elif response.status_code == 404:
                # CID not found, no point retrying
                return None
            else:
                # Other errors, retry
                time.sleep(sleep_time * (2 ** attempt))  # Exponential backoff
        except (requests.RequestException, requests.Timeout) as e:
            if attempt == max_retries - 1:
                print(f"Failed to fetch CID {cid} after {max_retries} attempts: {e}")
                return None
            time.sleep(sleep_time * (2 ** attempt))  # Exponential backoff
    
    return None
    

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

rxns_for_downstream_analysis_agg["rxn_cid"] = (
    rxns_for_downstream_analysis_agg["rxn_id"].astype(str).str.replace("'", "", regex=False)
)

reaction_id = pd.DataFrame(rxns_for_downstream_analysis_agg['rxn_cid'], columns=['rxn_cid'])

rxn_with_diff_format = convert_rxn_equ_batchmode(reaction_id, compound_lookup, 5000)
rxn_diff_euq_complete = rxn_with_diff_format[rxn_with_diff_format['incomplete_flag'] == False]

rxn_diff_euq_incomplete= rxn_with_diff_format[rxn_with_diff_format['incomplete_flag'] == True]['reaction_cid'].to_list()

missing_cids_combined= []
for rxn_id in rxn_diff_euq_incomplete:
    for rxn_cid_list in rxn_id.split('>>'):
        for cid in rxn_cid_list.split('|'):
            if cid.strip():  # Only add non-empty CIDs
                missing_cids_combined.append(cid.strip())
            
missing_cids= list(set(missing_cids_combined))

cid_to_smiles_mapping_dict= {}

if missing_cids:
    print(f"Fetching SMILES for {len(missing_cids)} missing CIDs...")
    for cid in missing_cids:
        smiles = get_smiles_from_cid(cid)
        if smiles:
            cid_to_smiles_mapping_dict[str(cid).strip()] = smiles

rxn_cid_to_rxn_smiles_converstion= []

for index, rxn_cid in enumerate(rxn_diff_euq_incomplete):
    print(f'Processing Index: {index+ 1}')
    reactants, products= rxn_cid.split('>>')
    reactant_cids = [int(cid.strip()) for cid in reactants.split('|')] if '|' in reactants else[int(reactants.strip())]
    product_cids= [int(cid.strip()) for cid in products.split('|')] if '|' in products else[int(products.strip())]
    reactant_smiles= [cid_to_smiles_mapping_dict.get(str(cid)) for cid in reactant_cids]
    product_smiles= [cid_to_smiles_mapping_dict.get(str(cid)) for cid in product_cids]
    if any(smi is None for smi in reactant_smiles + product_smiles):
        reaction_smiles = None
    else:
        reaction_smiles = str(".".join(reactant_smiles) + ">>" + ".".join(product_smiles))
    rxn_cid_to_rxn_smiles_converstion.append({rxn_cid: reaction_smiles})


rxn_smiles_data = [
    {'rxn_cid': list(d.keys())[0], 'reaction_smiles': list(d.values())[0]}
    for d in rxn_cid_to_rxn_smiles_converstion
]
rxn_smiles_df = pd.DataFrame(rxn_smiles_data)
(rxn_smiles_df['reaction_smiles'].isna()).sum()

missing_rxn_smiles= rxn_smiles_df[rxn_smiles_df['reaction_smiles'].notna()]
missing_rxns = missing_rxn_smiles.rename(columns={'rxn_cid': 'reaction_cid'})

combined = pd.concat([rxn_diff_euq_complete, missing_rxns], ignore_index=True, sort=False)
combined.to_csv('/home/mani/pubchem_rxns_2025_processing/rxn_diff_euq_complete_unprocessed.csv', index=False)
