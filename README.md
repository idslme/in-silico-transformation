# In-Silico Transformation Analysis (ISTA) Workflow

This repository provides a comprehensive pipeline for reaction data compilation, preprocessing, transformation modeling, and ADME-Tox evaluation of compounds.
---

## 📁 Contents

### 1. `pubchem_rxns_compilation.py`

Compiles chemical/biochemical reactions catalogued in PubChem from various sources:
- Generates unique reaction identifier
- Aggregates reaction information
- Removes incomplete entries
- Prepares the data for downstream processing

---

### 2. `rxn_equ_format_conversion.py`

Converts reaction CIDs into multiple reaction equation formats:
- **SMILES**
- **Text**
- **InChIKey**

This is a foundational step for downstream analysis and compatibility with cheminformatics tools.

---

### 3. `rxn_info_cli.py`

A command-line tool for extracting detailed information from reactions in **SMILES** format. Outputs include:
- Reaction class and name
- Functional groups (reactants and products)
- Reaction centers
- Transformation matrix

Supports **parallel processing** for large datasets.

---

### 4. `rxn_info_batch_processing.py`

Automates batch processing of large reaction datasets using `rxn_info_cli.py`. Features:
- Splits data into batches (default: 500 reactions)
- Utilizes up to **16 CPUs** for parallel execution

---

### 5. `reaction_template_generation.py`

Generates **reaction templates** from sanitized atom-mapped reactions. These templates capture the core chemical transformation pattern of reactions and are used in in-silico transformation modeling.

---
### 6. `substructure_generation.py`

Identifies and extracts **SMARTS substructures** of transformation sites:
- Detects atoms altered during the reaction
- Includes neighboring atoms for context
- Requires atom-mapped reactions as input

---

### 7. `biotransformation.py`

Simulates **in-silico biotransformations** of a query compound (input as SMILES):
- Matches reactive substructures using SMARTS
- Applies corresponding transformation templates
- Computes **ADME-Tox** properties for each transformation product

---
