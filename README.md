# In Silico Transformation

A comprehensive workflow for compiling reactions from PubChem, preprocessing reaction datasets, generating chemical transformation templates, identifying metabolic hotspots, and applying these templates to generate in silico biotransformation products.

---

## Overview

This repository provides tools for:

- Compiling chemical reactions from PubChem  
- Preprocessing and standardizing reaction datasets  
- Generating chemical transformation templates  
- Identifying sites of reactivity (metabolic hotspots) in reactants  
- Substructure analysis of reactants at transformation sites  
- Applying reaction templates to unknown compounds to generate predicted in silico transformation products  

---

## Features

- Automated reaction curation from PubChem  
- Batch reaction data preprocessing  
- Extraction and management of transformation templates  
- Site-of-reactivity (hotspot) identification  
- Fragmentation/substructure generation at metabolic sites  
- High-throughput application of transformation templates  
- Parallelization support for large-scale datasets  

---

## Requirements

- Python 3.7+
- [RDKit](https://www.rdkit.org/)
- pandas  
- numpy  
- requests  
- rdchiral  
- tqdm  

Install Python packages:
```bash
pip install pandas numpy requests tqdm
# Install RDKit: see https://www.rdkit.org/docs/Install.html
# Install rdchiral: https://github.com/connorcoley/rdchiral
```

---

## Installation

```bash
git clone https://github.com/idslme/in-silico-transformation.git
cd in-silico-transformation
```

---

## Workflow & Usage

Each major step is implemented as a standalone script. Run the scripts sequentially as needed:

1. **Compile Reactions from PubChem**  
   ```bash
   python 1_pubchem_rxns_compilation.py --input <input_parameters_file> --output <compiled_reactions_file>
   ```

2. **Convert rxn_cid to Reaction Equation Formats**  
   ```bash
   python 2_rxn_equ_format_conversion.py
   ```

3. **Extract Reaction Information (CLI)**  
   ```bash
   python 3_rxn_info_cli.py --input <input_reaction_file.csv> --output <output_info_file.csv>
   ```

4. **Batch Reaction Information Processing**  
   ```bash
   python 4_rxn_info_batch_processing.py --input <input_reaction_file.csv> --output <output_info_file.csv> --batch_size <N>
   ```

5. **Generate Reaction Templates**  
   ```bash
   python 5_reaction_template_generation.py --input <preprocessed_reactions.csv> --output <templates_output.csv>
   ```

6. **Substructure Generation at Transformation Sites**  
   ```bash
   python 6_substructure_generation.py
   ```

7. **In Silico Biotransformation (Product Prediction)**  
   ```bash
   python 7_biotransformation.py --input <input_compounds.csv> --templates <transformation_templates.csv> --output <biotransformation_products.csv>
   ```

---

## Notes

- Adjust file names/paths and arguments to match your data and workflow.
- Use `--help` with any CLI script to check for additional options:
  ```bash
  python <script>.py --help
  ```
- For input/output column and formatting details, see inline comments or function docstrings in each script.

---

## Directory Structure

```
in-silico-transformation/
├── 1_pubchem_rxns_compilation.py
├── 2_rxn_equ_format_conversion.py
├── 3_rxn_info_cli.py
├── 4_rxn_info_batch_processing.py
├── 5_reaction_template_generation.py
├── 6_substructure_generation.py
├── 7_biotransformation.py
├── README.md
```

---

## License

No license specified. Please add a license file if you wish to specify usage permissions.

---

## Citation

If you use this workflow in research or publication, please cite this repository:

```bibtex
@software{in_silico_transformation,
  author = {IDSLME Team},
  title = {In Silico Transformation: Reaction and Template Compilation Workflow},
  url = {https://github.com/idslme/in-silico-transformation},
  year = 2025
}
```

---

## Authors

Maintained by the **IDSLME group**.
