
import os
from multiprocessing import Pool
import json
import pandas as pd

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

target_directory = '/home/mani/pubchem-rxn-2025'

file_names = os.listdir(target_directory)

# compiling RHEA reactions
rhea_rxns_compiled = []

for file_name in file_names:
    if '_rhea' in file_name:
        with open(os.path.join(target_directory, file_name)) as f:
            data = json.load(f)
            df = pd.DataFrame(data['SDQOutputSet'][0]['rows'])
            df['file_name'] = file_name
            rhea_rxns_compiled.append(df)
            
rhea_rxns_compiled = pd.concat(rhea_rxns_compiled, ignore_index=True)
rhea_rxns_compiled.rename(columns={'ecs': 'enzyme'}, inplace=True)
rhea_rxns_compiled['source'] = 'RHEA:' + rhea_rxns_compiled['rhid'].astype(str)

# compiling transformation reactions
trans_rxns_compiled = []

for file_name in file_names:
    if '_trans' in file_name:
        with open(os.path.join(target_directory, file_name)) as f:
            data = json.load(f)
            df = pd.DataFrame(data['SDQOutputSet'][0]['rows'])
            df['file_name'] = file_name
            trans_rxns_compiled.append(df)

trans_rxns_compiled = pd.concat(trans_rxns_compiled, ignore_index=True)
trans_rxns_compiled.rename(columns={'predecessorcid': 'cidsreactant', 'successorcid': 'cidsproduct', 'sourcecomment': 'source'}, inplace=True)

# compiling pathway reactions
pwys_rxns_compiled = []

for file_name in file_names:
    if '_path' in file_name:
        with open(os.path.join(target_directory, file_name)) as f:
            data = json.load(f)
            df = pd.DataFrame(data['SDQOutputSet'][0]['rows'])
            df['file_name'] = file_name
            pwys_rxns_compiled.append(df)
            
pwys_rxns_compiled = pd.concat(pwys_rxns_compiled, ignore_index=True)
pwys_rxns_compiled.rename(columns={'ecs': 'enzyme', 'taxname': 'biosystem'}, inplace=True)


# Combine the three DataFrames
combined_rxns = pd.concat(
    [rhea_rxns_compiled[['cids', 'cidsreactant', 'cidsproduct', 'enzyme', 'source', 'file_name']],
     pwys_rxns_compiled[['cids', 'cidsreactant', 'cidsproduct', 'enzyme', 'biosystem', 'source', 'file_name']],
     trans_rxns_compiled[['cids', 'cidsreactant', 'cidsproduct', 'enzyme', 'biosystem', 'source', 'file_name']]],
    ignore_index=True,
    join='outer'  
)


# Reaction with Missing cid in Reactants or Products
combined_rxns_missing_cid = combined_rxns[
    combined_rxns['cidsreactant'].isna() | (combined_rxns['cidsreactant'] == '') |
    combined_rxns['cidsproduct'].isna() | (combined_rxns['cidsproduct'] == '') |
    combined_rxns['cids'].isna() | (combined_rxns['cids'] == '')
].copy()

# Removing Reactions with Missing cid in Reactants or Products
combined_rxns_missing_cid_rm = combined_rxns[
    combined_rxns['cidsreactant'].notna() & (combined_rxns['cidsreactant'] != '') &
    combined_rxns['cidsproduct'].notna() & (combined_rxns['cidsproduct'] != '') &
    combined_rxns['cids'].notna() & (combined_rxns['cids'] != '')
].copy()

# Generating Unique Reaction cid for Each Reaction
combined_rxns_missing_cid_rm['rxn_cid'] = rxncid(combined_rxns_missing_cid_rm)

# Filtering Reactions with Same Reactant and Product
combined_agg = combined_rxns_missing_cid_rm.groupby('rxn_cid')[['cids', 'cidsreactant', 'cidsproduct', 'enzyme', 'source', 'biosystem', 'file_name']].agg(
    lambda x: ';'.join([str(i) for i in set(x) if pd.notnull(i)])).reset_index()

combined_agg_same_react_product = combined_agg[combined_agg['cidsreactant'] == combined_agg['cidsproduct']]

# Removing Reactions with Same Reactant and Product
combined_agg_filtered = combined_agg[combined_agg['cidsreactant'] != combined_agg['cidsproduct']].copy()

# save the filtered DataFrame to a CSV file
combined_agg_filtered.to_csv('/home/mani/pubchem_rxns_2025_processing/combined_rxns_15may.csv', index=False)

