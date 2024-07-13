import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os
import json
import csv
import numpy as np
from functools import reduce
import glob
import math

  # Define file paths for easy access and modification
json_file = '/home/greg/ms-pred/quickstart/scarf/out/tree_preds_inten/pred_output.json'
tsv_file = '/home/greg/ms-pred/data/spec_datasets/sample_labels.tsv'


def process_smiles(smiles):
    # Create a molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

    # Check the correctness of the entered formula
    if mol is None:
        print("Error: Invalid SMILES formula entered.")
        return "Error: Invalid SMILES formula entered."

    # Convert the molecule to InChI
    inchi = Chem.MolToInchi(mol)

    # Convert InChI to InChIKey
    inchi_key = Chem.InchiToInchiKey(inchi)

    # Write the InChI string to the file
    with open(os.path.expanduser('~/rassp-public/rassp/sample_data/in.txt'), 'w') as file:
        file.write(inchi)

    print("SMILES:", smiles)
    print("InChI:", inchi)
    print("InChIKey:", inchi_key)
    print("The InChI string has been written to the file rassp_input.txt")

    # Read the existing example file to understand the format
    example_file_path = '~/massformer/predictions/example_smiles.csv'
    example_df = pd.read_csv(example_file_path)

    # Determine the next mol_id based on the last entry in the example file
    next_mol_id = example_df['mol_id'].max() + 1

    # Create a new DataFrame with the entered SMILES formula
    new_data = pd.DataFrame({'mol_id': [next_mol_id], 'smiles': [smiles]})

    # Write the new data to massformer_input.csv in the same format
    massformer_file_path = 'massformer_input.csv'
    new_data.to_csv(massformer_file_path, index=False)

    print(f"The SMILES formula has been written to the file {massformer_file_path} in the same format as example_smiles.csv.")

    # Update sample_labes.tsv
    sample_labes_file_path = '~/ms-pred/data/spec_datasets/sample_labels.tsv'
    sample_df = pd.read_csv(sample_labes_file_path, sep='\t')

    # Insert the SMILES formula into the specified columns and row
    sample_df.iloc[0, 1] = "output"  # Column 2 (index 1)
    sample_df.iloc[0, 5] = smiles  # Column 6 (index 5)

    # Save the updated DataFrame back to the file
    sample_df.to_csv(sample_labes_file_path, sep='\t', index=False)

    print(f"The SMILES formula has been inserted into the specified positions in the file {sample_labes_file_path}.")

    mass = calculate_molecular_mass(smiles)
    print(f"Molecular Mass of {smiles}: {mass}")

    # Create massformer input file
    create_massformer_input(smiles)

    # Run all processes
    flags = {'cmfid': True, 'massformer': True, 'rassp': True, 'scarf': True, 'iceberg': True}
    run_all(flags, smiles, mass)

    # Merge spectra files into a master CSV file
    merge_spects('/home/greg/MDU_outputs/output_*.csv', '/home/greg/MDU_outputs/master.csv')
    
    return "Process completed successfully."

def calculate_molecular_mass(smiles):
    # Convert SMILES to a molecule object
    molecule = Chem.MolFromSmiles(smiles)
    
    # Calculate the molecular mass
    mass = Descriptors.MolWt(molecule)
    
    return mass

def create_massformer_input(smiles):
    filename = "/home/greg/MDU_inputs/massformer_input1.csv"

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["mol_id", "smiles"])
        writer.writerow(["0", smiles])

    print(f"File '{filename}' created successfully.")

def run_conda_script(env_MDU, working_folder, MDU):
    os.system(f"conda run --no-capture-output -n {env_MDU} --cwd {working_folder} {MDU}")

def run_script(script):
    os.system(f"{script}")

def convert_json_to_csv():
    json_file_path = '/home/greg/ms-pred/quickstart/iceberg/out/tree_preds_inten/pred_output.json'
    tsv_file_path = '/home/greg/ms-pred/data/spec_datasets/sample_labels.tsv'
    
    with open(json_file_path, 'r') as file:
        data = json.load(file)

    rows = []
    for frag_hash, frag_info in data['frags'].items():
        row = {
            'frag': frag_info['frag'],
            'frag_hash': frag_info['frag_hash'],
            'parents': frag_info['parents'],
            'atoms_pulled': frag_info['atoms_pulled'],
            'left_pred': frag_info['left_pred'],
            'max_broken': frag_info['max_broken'],
            'tree_depth': frag_info['tree_depth'],
            'id': frag_info['id'],
            'prob_gen': frag_info['prob_gen'],
            'score': frag_info['score'],
            'form': frag_info['form'],
            'base_mass': frag_info['base_mass'],
            'frag_hs': frag_info['frag_hs'],
            'max_remove_hs': frag_info['max_remove_hs'],
            'max_add_hs': frag_info['max_add_hs'],
            'intens': frag_info['intens'],
            'mz_no_charge': frag_info['mz_no_charge'],
            'mz_charge': frag_info['mz_charge'],
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    root_inchi = data['root_inchi']
    name = data['name']
    df['root_inchi'] = root_inchi
    df['name'] = name

    cols = ['root_inchi', 'name'] + [col for col in df.columns if col not in ['root_inchi', 'name']]
    df = df[cols]

    df_selected = df[['name', 'prob_gen', 'base_mass']]

    tsv_data = pd.read_csv(tsv_file_path, sep='\t')

    merged_df = df_selected.merge(tsv_data[['spec', 'smiles']], left_on='name', right_on='spec', how='left')

    merged_df = merged_df.drop(columns=['name', 'spec'])
    merged_df = merged_df.rename(columns={'smiles': 'name'})

    merged_df['peaks'] = merged_df.apply(lambda row: f"({row['base_mass']}, {row['prob_gen']})", axis=1)

    peaks_combined = '[' + ', '.join(merged_df['peaks'].tolist()) + ']'

    final_df = pd.DataFrame({
        'mol_id': [1],
        'smiles': merged_df['name'].iloc[0],
        'prec_mz': merged_df['base_mass'].iloc[0],
        'nce': [''],
        'tool': ['ms-pred iceberg'],
        'peaks': [peaks_combined]
    })

    final_df.to_csv('/home/greg/MDU_outputs/output_iceberg.csv', index=False)

def convert_json_to_csv_with_titles(json_file, labels_df):
    with open(json_file) as f:
        data = json.load(f)
    
    root_inchi = data.get('cand_form', 'Unknown')
    spec_name = data.get('spec_name', 'Unknown')
    df = pd.DataFrame(data['output_tbl'])
    
    smiles_row = labels_df.loc[labels_df['spec'] == spec_name, 'smiles']
    if not smiles_row.empty:
        smiles = smiles_row.values[0]
    else:
        raise ValueError(f"Spec name {spec_name} not found in the TSV file.")
    
    df['Title'] = ''
    df.loc[0, 'Title'] = f'Root InChI: {root_inchi}'
    df.loc[1, 'Title'] = f'Name: {spec_name}'
    
    return df, smiles

def process_data(df, smiles):
    df['peaks'] = list(zip(df['mz'], df['rel_inten']))
    
    combined_peaks = [tuple(x) for x in df['peaks']]
    combined_peaks_str = str(combined_peaks)
    
    new_data = pd.DataFrame({
        'smiles': [smiles],
        'prec_mz': [df.loc[0, 'formula_mass_no_adduct']],
        'peaks': [combined_peaks_str]
    })
    
    new_data['mol_id'] = 1
    new_data['nce'] = np.nan
    new_data['tool'] = 'ms-pred scarf'
    
    new_data = new_data[['mol_id', 'smiles', 'prec_mz', 'nce', 'tool', 'peaks']]
    
    return new_data

def convert_json_to_final_csv(json_file, tsv_file):
    labels_df = pd.read_csv(tsv_file, sep='\t')
    
    df, smiles = convert_json_to_csv_with_titles(json_file, labels_df)
    
    final_df = process_data(df, smiles)
    
    output_file = os.path.join(os.getcwd(), '/home/greg/MDU_outputs/output_scarf.csv')
    final_df.to_csv(output_file, index=False)
    
    print(f"Final CSV file saved as {output_file}")

def combine_csv(input_mol_file, output_spectrum_file, output_file):
    mol_data = {}
    with open(input_mol_file, 'r', newline='') as input_mol_csv:
        reader = csv.DictReader(input_mol_csv)
        for row in reader:
            mol_data[row['mol_id']] = row['smiles']

    spectrum_data = []
    with open(output_spectrum_file, 'r', newline='') as output_spectrum_csv:
        reader = csv.DictReader(output_spectrum_csv)
        for row in reader:
            spectrum_data.append(row)

    with open(output_file, 'w+', newline='') as combined_csv:
        fieldnames = ['mol_id', 'smiles', 'prec_mz', 'nce', 'tool', 'peaks']
        writer = csv.DictWriter(combined_csv, fieldnames=fieldnames)
        writer.writeheader()

        for row in spectrum_data:
            mol_id = row['mol_id']
            smiles = mol_data.get(mol_id, '')
            writer.writerow({
                'mol_id': mol_id,
                'smiles': smiles,
                'prec_mz': row['prec_mz'],
                'nce': row['nce'],
                'tool': 'massformer',
                'peaks': row['peaks']
            })

import subprocess, csv, os
from pathlib import Path
from matchms.importing import load_from_mgf

ENERGY_LEVELS = {"energy0": 10, "energy1": 20, "energy2": 40}
def convert_from_mgf(mgffile, outfile):
    file = load_from_mgf(mgffile)
    with open(outfile, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["mol_id", "smiles", "prec_mz", "nce", "tool", "peaks"])
        for spect in file:
            title = spect.metadata['title'].split(";")
            line = [
                title[0] if title[0].lower() != 'nullid' else 0,
                title[3],
                spect.metadata["precursor_mz"],
                ENERGY_LEVELS[title[1].lower()],
                "CFM-ID",
                [(a[0], a[1]) for a in spect.peaks.to_numpy]
            ]
            writer.writerow(line)

def merge_spects(specs_path, output):
    file_dfs = [pd.read_csv(file, converters={'peaks': lambda s: {t[0]:t[1] for t in pd.eval(s)}}) for file in glob.glob(specs_path)]
    full_df = pd.concat([pd.DataFrame(columns=["smiles","tool","nce","peaks"]), *file_dfs], ignore_index=True)
    full_df["col_name"] = full_df["tool"] + full_df["nce"].map(lambda f:"" if math.isnan(f) else f" {f:.0f}")
    
    spec_dfs = full_df.apply(lambda row: pd.DataFrame.from_dict(row["peaks"], orient="index", columns=[row["col_name"]]), 
                             axis=1, result_type="reduce")
    
    spec_dfs = [(df / df.max() * 100).round(decimals=2) for df in spec_dfs]

    master_df = reduce(lambda l, r: l.join(r, how="outer"), spec_dfs, pd.DataFrame()).fillna(0).sort_index()
    master_df.rename_axis(full_df["smiles"].dropna().get(0, "m/z"), inplace=True)
    master_df.to_csv(output)

   

def run_all(flags, smiles, mass):
    os.chdir('/home/greg')

    if flags['cmfid']:
        run_script(f'docker run --rm=true -v .:/cfmid/public/ -i wishartlab/cfmid:latest sh -c "cd /cfmid/public/; cfm-predict \'{smiles}\' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 out.mgf"')
        import shutil
        shutil.move('./out.mgf','/home/greg/MDU_outputs/out.mgf')
        convert_from_mgf(f'/home/greg/MDU_outputs/out.mgf', '/home/greg/MDU_outputs/output_mgf.csv')
 
    if flags['massformer']:
        print("Running Massformer...")
        run_conda_script('MF-CPU','/home/greg/massformer','python scripts/run_inference.py -c config/demo/demo_eval.yml -s ~/MDU_inputs/massformer_input1.csv -o ~/MDU_outputs/massformer_output.csv -d -1')
        print("Massformer run completed.")
        combine_csv('/home/greg/MDU_inputs/massformer_input1.csv', '/home/greg/MDU_outputs/massformer_output.csv', '/home/greg/MDU_outputs/output_massformer.csv')
        print("Massformer output combined.")
   
    if flags['rassp']:
        run_conda_script('rassp2','~/rassp-public','USE_CUDA=0 python rassp/run_rassp.py rassp/sample_data/in.txt rassp/sample_data/out2.json --file-type inchi --model-name "FormulaNet" --no-gpu --output-file-type json')
        dfRaw = pd.read_json('/home/greg/rassp-public/rassp/sample_data/out2.json')

        dfResult = pd.DataFrame()    
        dfResult['mol_id'] = 1
        dfResult.loc[0] = [1]
        dfResult['peaks'] = repr(list(dfRaw['predictions'][0]['spect']))
        dfResult['nce'] = ''
        dfResult['tool'] = 'rassp'
        dfResult['smiles'] = smiles
        dfResult['prec_mz'] = mass

        dfResult = dfResult[['mol_id','smiles','prec_mz','nce','tool','peaks']]
        dfResult.to_csv('/home/greg/MDU_outputs/output_rassp.csv', index=False)

    if flags['scarf']:
        run_conda_script('ms-genNEW','~/ms-pred-new','. quickstart/scarf/run_model.sh')
        convert_json_to_final_csv(json_file, tsv_file)

    if flags['iceberg']:
        run_conda_script('ms-genNEW','~/ms-pred-new','. quickstart/iceberg/run_model.sh')
        convert_json_to_csv()
