
import os
import re
import html
import time
from typing import Dict, Tuple, Optional
import json
import pandas as pd
from flask import request
from httpx import get
import requests
from urllib.parse import quote
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
from streamlit import columns
from urllib.parse import quote

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

###### Fromat the Cleaned Reaction Definition ################################
SEP_PLUS = re.compile(r'\s+\+\s+')

def format_reaction_string(definition_clean: str):
    if '⟶' not in definition_clean:
        return definition_clean
    reactants, products = definition_clean.split('⟶')
    reactants_list = sorted([r.strip() for r in SEP_PLUS.split(reactants.strip()) if r.strip()])
    products_list = sorted([p.strip() for p in SEP_PLUS.split(products.strip()) if p.strip()])
    return ' + '.join(reactants_list) + ' ⟶ ' + ' + '.join(products_list)

## Checking Reaction CIDs are hyperlinked or Not
def safe_str(val):
    return '' if pd.isna(val) else str(val).strip()

def check_cid_hyperlink_connection(cids_reactant: str,
                                   cids_product: str,
                                   rxn_definition: str) -> Dict[str, str]:
    """
    Compare the number of chemical tokens on each side of a reaction definition
    with the number of CIDs provided for reactants/products.

    Args:
        cids_reactant: pipe-separated CIDs for reactants, e.g. "962|12345"
        cids_product: pipe-separated CIDs for products, e.g. "44237165|67890"
        rxn_definition: reaction text, e.g. "A + B ⟶ C + D"

    Returns:
        Dict with keys:
           - reactant_cid_missing_info
           - product_cid_missing_info
           - hyperlink_status
           - reactant_counts
           - product_counts
    """

    # --- Basic validation ---
    arrow = '⟶'
    if arrow not in rxn_definition:
        return {
            "reactant_cid_missing_info": "Reaction text missing arrow (⟶)",
            "product_cid_missing_info": "Reaction text missing arrow (⟶)",
            "hyperlink_status": "Unknown",
            "reactant_counts": "0/0",
            "product_counts": "0/0",
        }

    # --- Split reaction into sides ---
    reactants_text, products_text = rxn_definition.split(arrow, 1)

    # --- Tokenize chemicals: split on separator '+' but not charges ---
    reactant_tokens = [t for t in SEP_PLUS.split(reactants_text.strip()) if t.strip()]
    product_tokens  = [t for t in SEP_PLUS.split(products_text.strip()) if t.strip()]
    reactants_definition_len = len(reactant_tokens)
    products_definition_len  = len(product_tokens)

    # --- Count provided CIDs (handle empty strings safely) ---
    reactant_cids = [c for c in (cids_reactant or '').split('|') if c.strip()]
    product_cids  = [c for c in (cids_product  or '').split('|') if c.strip()]
    reactant_cids_len = len(reactant_cids)
    product_cids_len  = len(product_cids)

    # --- Side-by-side messages ---
    reactant_cid_missing_info = ("No Missing CIDs in Reactant Side"
                                 if reactants_definition_len == reactant_cids_len
                                 else "Missing CIDs in Reactant Side")

    product_cid_missing_info = ("No Missing CIDs in Product Side"
                                if products_definition_len == product_cids_len
                                else "Missing CIDs in Product Side")

    # --- Overall hyperlink status ---
    if (reactant_cid_missing_info == "No Missing CIDs in Reactant Side"
        and product_cid_missing_info == "No Missing CIDs in Product Side"):
        hyperlink_status = "Reaction CID is Hyperlinked"
    else:
        hyperlink_status = "Reaction CID is Not Hyperlinked"

    return {
        "reactant_cid_missing_info": reactant_cid_missing_info,
        "product_cid_missing_info": product_cid_missing_info,
        "hyperlink_status": hyperlink_status,
        "reactant_counts": f"{reactant_cids_len}/{reactants_definition_len}",
        "product_counts": f"{product_cids_len}/{products_definition_len}",
    }
    
    
# Counting Number of Unique EC for Each Reaction Source
def count_unique_ecs_for_source(combined_rxns, source_name):
    # Filter for the given source
    df = combined_rxns[combined_rxns['source'] == source_name]
    # Drop NaN in 'enzyme' column
    ecs_series = df['enzyme'].dropna()
    # Split by '|' and flatten
    ecs_split = ecs_series.str.split('|').explode()
    # Remove empty strings and strip whitespace
    ecs_clean = ecs_split[ecs_split.str.strip() != ''].str.strip()
    # Remove NaN again (in case splitting introduced any)
    ecs_clean = ecs_clean.dropna()
    # Get unique EC numbers
    unique_ecs = set(ecs_clean)
    return len(unique_ecs) 


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
rhea_rxns_compiled.rename(columns={'ecs': 'enzyme', 'equation': 'definition'}, inplace=True)
rhea_rxns_compiled['definition'] = rhea_rxns_compiled['definition'].str.replace(' = ', ' ⟶ ', regex=False)
rhea_rxns_compiled['source'] = 'RHEA'

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
trans_rxns_compiled['definition'] = (
    trans_rxns_compiled['predecessor'].fillna('') + 
    ' ⟶ ' + 
    trans_rxns_compiled['successor'].fillna('')
)

trans_rxns_compiled.rename(columns={'predecessorcid': 'cidsreactant', 'successorcid': 'cidsproduct'}, inplace=True)
trans_rxns_compiled['source'] = trans_rxns_compiled['sourcecomment'].apply(
    lambda x: x.split('|')[0].strip() if 'MetXBioDB |' in x else x
)
trans_rxns_compiled['source'] = trans_rxns_compiled['source'].apply(
    lambda x: x.split('_')[0] if 'PPDB_ID_' in x else x
)
mask = trans_rxns_compiled['source'] == ''
trans_rxns_compiled.loc[mask, 'source'] = trans_rxns_compiled.loc[mask, 'evidenceref']


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

combined_rxns= pd.concat([rhea_rxns_compiled, pwys_rxns_compiled, trans_rxns_compiled], ignore_index= True, join= 'outer')
combined_rxns['rxn_text']= combined_rxns['definition'].apply(strip_html_formatting)

hyperlink_check_results= [
                check_cid_hyperlink_connection(
                    safe_str(row['cidsreactant']),
                    safe_str(row['cidsproduct']),
                    safe_str(row['rxn_text'])
                )
                for _, row in combined_rxns.iterrows()
            ]

# Combine results with original data
combined_rxn_hyperlink_checked= pd.concat(
    [combined_rxns.reset_index(drop=True), 
        pd.DataFrame(hyperlink_check_results)], 
    axis=1)

column_reordering_list= ['taxid', 'cids', 'name', 'definition', 'reaction', 'control', 
                            'cidsreactant', 'cidsproduct', 'url', 'source', 
                            'externalid', 'file_name', 'geneids', 'protacxns', 'dois', 
                            'pmids', 'pclids', 'citations', 'pmcids', 'enzyme', 'gid', 'srcid', 
                            'predecessor', 'transformation', 'successor', 'evidencedoi', 
                            'evidenceref', 'sourcecomment', 'sourcecommentfull', 'datasetdoi', 
                            'datasetref', 'biosystem', 'rhid', 'htmlequation', 'direction', 'otherdirections',
                            'rxn_text', 'reactant_counts', 'product_counts', 
                            'reactant_cid_missing_info', 'product_cid_missing_info', 'hyperlink_status']

# Select columns and create a copy for formatting
table_1_final = combined_rxn_hyperlink_checked[column_reordering_list].copy()

# Clean HTML formatting from text columns
for col in ['definition', 'rxn_text', 'enzyme', 'source']:
    if col in table_1_final.columns:
        table_1_final[col] = table_1_final[col].apply(
            lambda x: strip_html_formatting(str(x)) if pd.notna(x) else ''
        )

# Format text columns to force Excel text mode (single quote prefix for text recognition)
# This applies to CID counts and text that should not be auto-formatted
for col in ['reactant_counts', 'product_counts', 'rxn_text', 'cidsreactant', 'cidsproduct']:
    if col in table_1_final.columns:
        table_1_final[col] = "'" + table_1_final[col].astype(str)

# Save to CSV with proper encoding
table_1_final.to_csv('/home/mani/ista_revision_final/Zenodo_Table_1.csv', index=False, quoting=csv.QUOTE_ALL, encoding='utf-8-sig')
print("Table 1 saved with HTML cleaned and text columns properly formatted for Excel")


############# SI Table 1_Generation #################################

table_S1 = (table_1_final['source']
    .value_counts()
    .rename_axis('source')
    .reset_index(name='reaction_count')
)

unique_rxn_counts = (
    table_1_final.groupby('source')['rxn_text']
    .nunique()
    .reset_index(name='unique_rxn_count')
    .sort_values('unique_rxn_count', ascending=False)
)

# Filter for human biosystem
human_mask = (
    table_1_final['biosystem'].fillna('').str.contains('Homo sapiens (human)', case=False, na=False, regex=False)
    | (table_1_final['biosystem'].fillna('').str.lower() == 'human')
)
human_reactions = table_1_final[human_mask]
human_reaction_count= (human_reactions['source']
                       .value_counts()
                       .rename_axis('source')
                       .reset_index(name= 'human_unique_rxn_counts')
                       )


# Count unique reactions per source for human biosystem
human_unique_rxn_counts = (
    human_reactions.groupby('source')['rxn_text']
    .nunique()
    .reset_index(name='unique_human_rxn_count')
    .sort_values('unique_human_rxn_count', ascending=False)
)


# Combine everything into table_S1
# Merge all counts into a single DataFrame with desired column names

table_S1 = (
    table_S1
    .merge(unique_rxn_counts.rename(columns={'unique_rxn_count': 'Unique Reactions by Source'}), on='source', how='left')
    .merge(human_reaction_count.rename(columns={'human_unique_rxn_counts': 'Reaction in Homo Sapiens'}), on='source', how='left')
    .merge(human_unique_rxn_counts.rename(columns={'unique_human_rxn_count': 'Unique Reactions in Homo Sapiens'}), on='source', how='left')
)

# Add unique EC counts per source (using combined_rxns to keep all enzyme entries)
table_S1['No of Unique EC Numbers'] = table_S1['source'].apply(
    lambda src: count_unique_ecs_for_source(combined_rxns, src)
)

# Rename count column to match requested output
if 'reaction_count' in table_S1.columns:
    table_S1 = table_S1.rename(columns={'reaction_count': 'count'})


rxs_not_hyperlinked = table_1_final[table_1_final['hyperlink_status'] == 'Reaction CID is Not Hyperlinked']

missing_name_to_cid_counts = pd.DataFrame({
    'No of Reactions': table_1_final['source'].value_counts(dropna=False),
    'No of Raections Missing Name to CID Map': rxs_not_hyperlinked['source'].value_counts(dropna=False)
}).fillna(0).astype(int).reset_index()
missing_name_to_cid_counts.columns = ['source', 'No of Reactions', 'No of Raections Missing Name to CID Map']

missing_name_to_cid_counts['Percentage Missing (%)'] = (missing_name_to_cid_counts['No of Raections Missing Name to CID Map'] / missing_name_to_cid_counts['No of Reactions'] * 100).round(2)

# Merge missing counts into table_S1
table_S1 = table_S1.merge(
    missing_name_to_cid_counts[['source', 'No of Raections Missing Name to CID Map', 'Percentage Missing (%)']],
    on='source',
    how='left'
)

# Reorder columns for final output
table_S1 = table_S1[[
    'source', 'count', 'Unique Reactions by Source',
    'Reaction in Homo Sapiens', 'Unique Reactions in Homo Sapiens',
    'No of Unique EC Numbers', 'No of Raections Missing Name to CID Map', 'Percentage Missing (%)'
]]

table_S1.to_csv('/home/mani/ista_revision_final/table_S1.csv', index= False)

#### SI Table 2 Generation #######################################
biosystem_counts = table_1_final['biosystem'].value_counts().reset_index()

biosystem_sources = (
    table_1_final.groupby('biosystem')['source']
    .apply(lambda x: '|'.join(sorted(set(x.dropna()))))
    .reset_index(name='sources')
)

# Merge with biosystem_counts for a full summary
biosystem_counts = biosystem_counts.merge(biosystem_sources, on='biosystem', how='left')

biosystem_unique_rxn_counts = (table_1_final.groupby('biosystem')['rxn_text']
    .nunique()
    .reset_index(name='unique_rxn_count_per_biosystem')
    .sort_values('unique_rxn_count_per_biosystem', ascending=False)
)


missing_name_to_cid_counts_biosys = pd.DataFrame({
    'No of Reactions': table_1_final['biosystem'].value_counts(dropna=False),
    'No of Raections Missing Name to CID Map': rxs_not_hyperlinked['biosystem'].value_counts(dropna=False)
}).fillna(0).astype(int).reset_index()
missing_name_to_cid_counts_biosys.columns = ['biosystem', 'No of Reactions', 'No of Raections Missing Name to CID Map']

missing_name_to_cid_counts_biosys['Percentage Missing (%)'] = (missing_name_to_cid_counts_biosys['No of Raections Missing Name to CID Map'] / missing_name_to_cid_counts_biosys['No of Reactions'] * 100).round(2)

# Merge all biosystem data for Table S2
table_S2 = (
    biosystem_counts
    .merge(biosystem_unique_rxn_counts, on='biosystem', how='left')
    .merge(missing_name_to_cid_counts_biosys[['biosystem', 'No of Raections Missing Name to CID Map', 'Percentage Missing (%)']], 
           on='biosystem', how='left')
)

# Reorder columns for final output
table_S2 = table_S2[[
    'biosystem', 'count', 'unique_rxn_count_per_biosystem', 'sources',
    'No of Raections Missing Name to CID Map', 'Percentage Missing (%)'
]]

# Rename columns for clarity
table_S2 = table_S2.rename(columns={
    'count': 'No of Reactions',
    'unique_rxn_count_per_biosystem': 'Unique Reactions per Biosystem',
    'sources': 'Reaction Sources'
})

table_S2.to_csv('/home/mani/ista_revision_final/table_S2.csv', index=False)
print(f"Table S2 saved with {len(table_S2)} biosystems")
