"""Generating Substructure of Transformation Site of Reactants Involving in the Chemical Reaction"""
import logging
import os
import re
import pandas as pd
from rdchiral import template_extractor
from rdkit import Chem
from multiprocessing import Pool

def configure_logging(input_file_path: str):
    """Configure logging to save log file in the same directory as the input file."""
    log_dir = os.path.dirname(input_file_path)  
    log_file_path = os.path.join(log_dir, 'reaction_processing.log')    
    
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
        reaction = {
            'reactants': rxn.split('>>')[0],
            'products': rxn.split('>>')[1],
            '_id': rxn_id
        }
                
        reactant_mols = [Chem.MolFromSmiles(reactant) for reactant in reaction['reactants'].split('.')]
        product_mols = [Chem.MolFromSmiles(product) for product in reaction['products'].split('.')]      
        
        if None in reactant_mols or None in product_mols:
            raise ValueError(f"Error in converting SMILES to molecules for reaction {rxn_id}.")
        
        logging.info(f"Successfully processed reaction {rxn_id}.")
        return reaction, reactant_mols, product_mols
    except Exception as e:
        logging.error(f"Error processing reaction {rxn_id}: {e}")
        return None, [], []

def extract_extended_map_numbers(smarts):
    """Extract atom map numbers from a SMARTS pattern."""
    try:        
        map_numbers = re.findall(r':(\d+)', smarts)        
        return list(map(int, map_numbers))
    except Exception as e:
        print(f"Error extracting map numbers: {e}")
        return None
    
def generate_fragment_smiles(mol_smiles, extended_atom_tags):
    """
    Generate fragment SMILES using MolFragmentToSmiles based on extended atom tags.

    Args:
        mol_smiles (str): The SMILES string of the molecule.
        extended_atom_tags (list): List of extended atom map numbers to include in the fragment.

    Returns:
        str: The fragment SMILES.
    """
    try:        
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            raise ValueError("Invalid molecule SMILES.")
        
        atom_indices = []
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() in extended_atom_tags:
                atom_indices.append(atom.GetIdx())

        if not atom_indices:
            raise ValueError("No atoms found with the specified atom map numbers.")
        
        fragment_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=atom_indices, canonical=True, allHsExplicit=True, allBondsExplicit=True)
        return fragment_smiles
    except Exception as e:
        print(f"Error generating fragment SMILES: {e}")
        return None


def get_transformed_atoms_substructure(rxn: str, rxn_id: str):
    """Get the substructure of the transformed atoms in the reactants."""
    reactants_trans_info = []
    reaction, reactant_mols, product_mols = mapped_rxn_to_rdchiral_input(rxn, rxn_id)
    
    if not reaction:
        return None
    
    try:
        
        
        changed_atoms, changed_atom_tags, err = template_extractor.get_changed_atoms(reactant_mols, product_mols)
        if changed_atom_tags is not None:
            for mol in reactant_mols:
                atoms, atom_tags = template_extractor.get_tagged_atoms_from_mol(mol)
                changed_mapnum = [tag for tag in atom_tags if tag in changed_atom_tags]
                changed_atoms_smarts = [atom.GetSmarts() for atom in atoms if str(atom.GetAtomMapNum()) in changed_mapnum]
                changed_atoms_index = [atom.GetIdx() for atom in atoms if str(atom.GetAtomMapNum()) in changed_mapnum]
                fragments, intra_only, dimer_only = template_extractor.get_fragments_for_changed_atoms([mol], changed_mapnum, radius=1, category='reactants')
                extended_map_numbers = extract_extended_map_numbers(fragments)
                extended_atom_indices = [atom.GetIdx() for atom in atoms if atom.GetAtomMapNum() in extended_map_numbers]
                smiles_fragment = Chem.MolFragmentToSmiles(mol, atomsToUse=extended_atom_indices, canonical=True)

                reactant_info = {
                    'rxn_id': rxn_id,
                    'Smiles': Chem.MolToSmiles(mol),
                    'Changed Atom SMARTS': changed_atoms_smarts,
                    'Changed Atom Tags': changed_mapnum,
                    'Extended Changed Atom Tags': extended_map_numbers,
                    'Changed Atom Index': changed_atoms_index,
                    'Extended Atom Index': extended_atom_indices,
                    'Fragmants Smarts': fragments,
                    'Fragment SMILES': smiles_fragment,
                }
                             
                reactants_trans_info.append(reactant_info)
            return reactants_trans_info
    except Exception as e:
        logging.error(f"Error processing reaction {rxn_id}: {e}")
        return None


def process_reactions_in_parallel(input_file_path: str, output_file_path: str, max_workers=None):
    """
    Process reactions in parallel to extract transformed atom substructures.

    Args:
        input_file_path (str): Path to the input CSV file containing reactions.
        output_file_path (str): Path to save the output CSV file with results.

    Returns:
        None
    """
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    try:        
        input_data = pd.read_csv(input_file_path)
                
        required_columns = ['SANITIZED_MAPPED_REACTION', 'RXN_ID']
        if not all(col in input_data.columns for col in required_columns):
            raise ValueError(f"Input file must contain the following columns: {required_columns}")
        
        args = [(row['SANITIZED_MAPPED_REACTION'], row['RXN_ID']) for _, row in input_data.iterrows()]
        
        if max_workers is None:
            max_workers = max(1, os.cpu_count() // 2)            
        with Pool(processes=max_workers) as pool:
            results = pool.starmap(get_transformed_atoms_substructure, args)
        
        processed_results = [result for result in results if result is not None]
        
        flattened_results = [item for sublist in processed_results for item in sublist]
        
        output_df = pd.DataFrame(flattened_results)
        output_df.to_csv(output_file_path, index=False)
        logging.info(f"Successfully processed reactions. Results saved to {output_file_path}")

    except Exception as e:
        logging.error(f"Error processing reactions: {e}")


if __name__ == "__main__":
    
    input_file_path = '/home/mani/missing_reaction_processing/combined_info_with_templates.csv'
    output_file_path =  '/home/mani/missing_reaction_processing/transformed_atom_substructure.csv'
    process_reactions_in_parallel(input_file_path, output_file_path, max_workers=16)
