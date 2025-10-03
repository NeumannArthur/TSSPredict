# TSSprediction from RNAseq
## Date
2025-08-04
## Project Description
TSS prediction based solely on RNAseq data. 
Using ML methods (NN+gLM) to embed and predict data.
## Project Organization
The project is organized in the following directories:
- src/: scripts that are connected to a GitHub repo
- data/: raw data, cleaned data, etc.
- analysis/: scripts that are not connected to a GitHub repo
- figures/: plots, graphs, etc.
- playground/: for testing purposes
- results/: results of the project
## Project Team
- Student: Arthur Neumann
- Supervisor: Maik Wolfram-Schauerte + Mathias Witte Paz


# CLI Documentation

### Top-Level Commands

- `data`  
    Data-related operations. Contains subcommands for cleaning, viewing, and windowing datasets.

- `train`
    No current implementation, planned for future ML implementations


#### Subcommands under `data`

##### 1. `clean`
- Cleans the raw dataset and outputs processed files.
- Usage Example:
    ```bash
    python3 main.py data clean
    ```

##### 2. `view`
- Previews cleaned or raw datasets in the terminal, with options to save a preview.
- Options:
    - `--dataset {Master,Clean}`: Choose which dataset to view (`Master` for raw, `Clean` for processed). Default: `Clean`.
    - `--rows <int>`: Number of rows to preview. Default: `5000`.
    - `--save <path>`: Save the preview to a `.tsv` file at the specified path.
- Usage Example:
    ```bash
    python3 main.py data view --dataset Clean --rows 1000 --save /path/to/save
    ```

##### 3. `window`
- Performs window-based analysis on the cleaned dataset.
- Options:
    - `--size <int>`: Window size in each direction from TSS. Default: `1000`.
    - `--conds <list>`: List of conditions to filter by. Default: `Control`.
                            Choices: Control, nov, rif, tet
    - `--strand {+, -, master, all}`: Specify strand for windowing. Default: all.
    - `--preview <path>`: [ Optional ] Path to save preview output.
    - `--wiggle_f <path>`: [ Optional ] Path to forward wiggle file.
    - `--wiggle_b <path>`: [ Optional ] Path to backward wiggle file.
    - `--output_dir <path>`: Output directory for results.
- Usage Example:
    ```bash
    python3 main.py data window --size 1000 --conds Control nov --strand + --preview /path/to/preview
    ```




