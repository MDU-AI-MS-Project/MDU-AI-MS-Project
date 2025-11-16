#!/usr//bin/env python3
import csv
import importlib.util
from typing import Union
from pathlib import Path
import datetime
import time
import shutil
import subprocess
import json
from rdkit.Chem import MolFromSmiles, MolToSmiles, Descriptors
import pandas as pd
import msp_utils
from table_aggregator import TableAggregator

rootdir = Path("/mnt/c/MDU")
csv_path = rootdir / 'GNPSnew/GNPSnew_canonical_interim.csv'
output_path = rootdir / 'GNPSnew/GNPSnew_canonical_output.csv'
output_path_split = rootdir / 'GNPSnew/GNPSnew_canonical_output_split.csv'
output_path_mod1 = rootdir / 'GNPSnew/GNPSnew_canonical_output_mod1.csv'
output_path_mod2 = rootdir / 'GNPSnew/GNPSnew_canonical_output_mod2.csv'
descriptor_keys = set()

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def get_smiles_based_folder(base_path, smiles):
    #timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    #folder_name = f"{stem}_{timestamp}"
    folder_name = f"{smiles}"
    new_folder = base_path / folder_name
    if new_folder.exists():
        print(f"Folder already exists {new_folder}")
        return new_folder
    else:
        print(f"Creating folder: {new_folder}")
        new_folder.mkdir(parents=True, exist_ok=True)            
        try:
            new_folder.mkdir(parents=True, exist_ok=True)
            print(f"Created folder: {new_folder}")
            return new_folder
        except Exception as e:
            print(f"Error creating folder: {e}")
            exit(1)

def move_files(source_path, destination_folder, pattern="*"):
    try:
        for file in source_path.glob(pattern):
            if file.is_file():
                shutil.move(str(file), str(destination_folder))
                print(f"Moved {file.name} to {destination_folder}")
    except Exception as e:
        print(f"Error moving files: {e}")
        exit(2)

def process_massbank_reference_spectra():
    msp_spectra = msp_utils.read_msp_file(msp_archive)
    massbank = msp_utils.SpectraLookup(msp_spectra)
    for smiles in massbank.smiles_set:
        target_folder = get_smiles_based_folder(rootdir, smiles)
        descriptors_output = target_folder / ("descriptors_" + smiles)
        ref_spectrum_output = target_folder / ("ref_spectrum_" + smiles)
        #run_models_module.process_smiles(smiles)
        #move_files(rootdir, target_folder,"output_*.csv")
        #move_files(rootdir, target_folder,"master.csv")
        subprocess.run([
            'conda', 'run', '-n', 'rdkit-dev',
            'python3', descriptors_script,
            smiles,
            descriptors_output
        ])
        with open(ref_spectrum_output, 'w') as f:
            json.dump(["0",smiles,"500","10","REF",massbank.get_merged_normalized_cutoff_binned_spectra(smiles)], f, indent=2)

def process_referense_spectra(input_file, aggregator: TableAggregator):
    input_file = Path(input_file)
    spectra_list = []
    indices_to_remove = []

    if 'reference_spectra' in aggregator.df.columns:
        return

    extractor = msp_utils.SpectraLookup(msp_utils.read_reference_file(input_file, aggregator))
    
    not_in_dictionary = 0
    not_primary = 0
    for index, row in aggregator.df.iterrows():
        smiles = str(row['canonical_smiles']).strip()
        if not smiles in extractor.smiles_set:
            print(f"Smiles not found in dictionary {smiles}")
            not_in_dictionary += 1
            #spectra_list.append([0.0 for i in range(2800)])
            indices_to_remove.append(index)
        elif not row['is_primary']:
            print(f"Smiles is not primary {smiles}")
            not_primary += 1
            indices_to_remove.append(index)
        else:
            spectra_list.append(extractor.get_merged_normalized_cutoff_binned_spectra(smiles))
    aggregator.df = aggregator.df.drop(indices_to_remove)
    print(f"Removed {not_in_dictionary} + {not_primary} rows for smiles not found in reference spectra")
    aggregator.df['reference_spectra'] = spectra_list


def process_cfmid_predictions(input_folder, aggregator: TableAggregator):
    input_folder = Path(input_folder)
    #output_root_folder = Path(output_root_folder)
    #smiles_csv = Path(smiles_ordered_csv)
    #if 'cfmid_spectra' not in aggregator.df.columns:
    #    aggregator.df['cfmid_spectra'] = None
    
    if 'cfmid_spectra' in aggregator.df.columns:
        return
    
    spectra_list = []
    spectra_dict = {}
    for index, row in aggregator.df.iterrows():
        if not row['is_primary']:
            spectra_list.append([0.0 for i in range(2800)])
            continue
        #with open(smiles_csv, 'r') as csvfile:
        #found = 0
        smiles = str(row['smiles']).strip()
        canonical_smiles = str(row['canonical_smiles']).strip()
        #csv_reader = csv.reader(csvfile)
        #for smiles_index, row in enumerate(csv_reader):
            # Assume the string is in the second field (index 1)
            #if len(row) < 2:
                #print(f"Warning: Row {row} has insufficient fields")
                #continue
            #if "mol_id" in row:
                #continue
            #smiles = row[1].strip()
        input_file = input_folder / f"GNPSnew_cfmpredict_{index + 1}.txt"
        print (f"Processing {input_file}")
        if not input_file.exists():
            print(f"Warning: File {input_file} does not exist. Skipping SMILES {smiles}")
            spectra_list.append([0.0 for i in range(2800)])
            continue
        else:
            cfmid_spectra = msp_utils.read_cfmid_file(input_file, smiles)
            if cfmid_spectra:
                #found += 1
                #output_folder = output_root_folder / str(smiles_index)
                #output_folder.mkdir(parents=True, exist_ok=True)
                #output_file = output_folder / "output_cfmid.csv"
                extractor = msp_utils.SpectraLookup(cfmid_spectra)
                #with open(output_file, 'w') as f:
                #    json.dump(["0",smiles,"500","10","CFMID",extractor.get_merged_normalized_cutoff_binned_spectra(smiles)], f, indent=2)
                #aggregator.df.at[index, 'cfmid_spectra'] = extractor.get_merged_normalized_cutoff_binned_spectra(canonical_smiles)
                #aggregator.add_binned_spect_as_row(extractor.get_merged_normalized_cutoff_binned_spectra(canonical_smiles), index)
                
                #spectra_dict[index] = extractor.get_merged_normalized_cutoff_binned_spectra(canonical_smiles)
                spectra_list.append(extractor.get_merged_normalized_cutoff_binned_spectra(canonical_smiles))
            else:
                print(f"Warning: No CFMID data found for {smiles}")
                spectra_list.append([0.0 for i in range(2800)])
    #column_names = [f"cfmid_{i}" for i in range(len(spectra_dict[1]))]
    #aggregator.df = pd.concat([aggregator.df, pd.DataFrame.from_dict(spectra_dict, orient='index', columns=column_names)], axis=1)
    aggregator.df['cfmid_spectra'] = spectra_list

def process_scarf_predictions(input_folder, aggregator: TableAggregator):
    input_folder = Path(input_folder)
    
    #if 'scarf_spectra' not in aggregator.df.columns:
    #    aggregator.df['scarf_spectra'] = None

    if 'scarf_spectra' in aggregator.df.columns:
        return

    spectra_list = []
    for index, row in aggregator.df.iterrows():
        if not row['is_primary']:
            spectra_list.append([0.0 for i in range(2800)])
            continue

        smiles = str(row['canonical_smiles']).strip()
        input_file = input_folder / f"tree_preds_inten_chunk-{(index // 49):04d}" / f"pred_{index + 1}.json"
        print (f"Processing {input_file}")
        if not input_file.exists():
            print(f"Warning: File {input_file} does not exist. Skipping SMILES {smiles}")
            spectra_list.append([0.0 for i in range(2800)])
            continue
        else:
            scarf_spectra = msp_utils.read_scarf_file(input_file)
            extractor = msp_utils.SpectraLookup(scarf_spectra)
            if scarf_spectra and smiles in extractor.smiles_set:
                spectra_list.append(extractor.get_merged_normalized_cutoff_binned_spectra(smiles))
            else:
                print(f"Warning: No SCARF data found for {smiles}")
                spectra_list.append([0.0 for i in range(2800)])

    aggregator.df['scarf_spectra'] = spectra_list


def process_iceberg_predictions(input_folder, aggregator: TableAggregator):
    input_folder = Path(input_folder)
    
    #if 'iceberg_spectra' not in aggregator.df.columns:
    #    aggregator.df['iceberg_spectra'] = None
    
    if 'iceberg_spectra' in aggregator.df.columns:
        return

    spectra_list = []
    probability_list = []
    for index, row in aggregator.df.iterrows():
        if not row['is_primary']:
            spectra_list.append([0.0 for i in range(2800)])
            probability_list.append([0.0 for i in range(2800)])
            continue

        smiles = str(row['canonical_smiles']).strip()
        input_file = input_folder / f"tree_preds_inten_chunk-{(index // 49):04d}" / f"pred_{index + 1}.json"
        print (f"Processing {input_file}")
        if not input_file.exists():
            print(f"Warning: File {input_file} does not exist. Skipping SMILES {smiles}")
            spectra_list.append([0.0 for i in range(2800)])
            probability_list.append([0.0 for i in range(2800)])
            continue
        else:
            iceberg_spectra, iceberg_probability = msp_utils.read_iceberg_file(input_file)
            extractor = msp_utils.SpectraLookup(iceberg_spectra)
            probability_extractor = msp_utils.SpectraLookup(iceberg_probability)
            if iceberg_spectra and smiles in extractor.smiles_set:
                spectra_list.append(extractor.get_merged_normalized_cutoff_binned_spectra(smiles))
                probability_list.append(probability_extractor.get_merged_normalized_cutoff_binned_spectra(smiles))
            else:
                print(f"Warning: No ICEBERG data found for {smiles}")
                spectra_list.append([0.0 for i in range(2800)])
                probability_list.append([0.0 for i in range(2800)])

    aggregator.df['iceberg_spectra'] = spectra_list
    aggregator.df['iceberg_probability'] = probability_list

def process_rassp_predictions(input_folder: Union[str, Path], aggregator: TableAggregator):
    input_folder = Path(input_folder)
    rassp_spectra = []
    spectra_list = []
    hit_map = {}

    if 'rassp_spectra' in aggregator.df.columns:
        return
    
    for input_file in input_folder.glob("GNPSnew_rassp.chunk-*.json"):
        print (f"Processing {input_file}")
        rassp_spectra.extend(msp_utils.read_rassp_file(input_file))
 

    extractor = msp_utils.SpectraLookup(rassp_spectra)
    
    for index, row in aggregator.df.iterrows():
        if not row['is_primary']:
            spectra_list.append([0.0 for i in range(2800)])
            continue
        smiles = str(row['canonical_smiles']).strip()
        print (f"Processing row: {index} smiles: {smiles}")
        hit_map[smiles] = 0
        if not smiles in extractor.smiles_set:
            print(f"Smiles not found in dictionary {smiles}")
            spectra_list.append([0.0 for i in range(2800)])
        else:
            hit_map[smiles] = extractor.smiles_set.count(smiles)
            spectra_list.append(extractor.get_merged_normalized_cutoff_binned_spectra(smiles))
    aggregator.df['rassp_spectra'] = spectra_list

    hits = 0
    multihits = 0
    for smiles, count in hit_map.items():
        if count == 0:
            print(f"Warning: No RASSP data found for {smiles}")
            continue
        else:
            hits +=1
        
        if count > 1:
            print(f"Warning: Multiple RASSP data found for {smiles}")
            multihits += 1
    print(f"Found RASSP data for {hits} SMILES")
    print(f"Found multiple RASSP data for {multihits} SMILES")

def process_massformer_predictions(input_folder: Union[str, Path] , aggregator: TableAggregator):
    massformer_spectra = []
    spectra_list = []
    hit_map = {}
    input_folder = Path(input_folder)
    
    if 'massformer_spectra' in aggregator.df.columns:
        return

    for input_file in input_folder.glob("GNPSnew_chunk*.csv"):
        print (f"Massformer - Processing {input_file}")
        massformer_spectra.extend(msp_utils.read_massformer_file(input_file, aggregator))

    extractor = msp_utils.SpectraLookup(massformer_spectra)
    
    for index, row in aggregator.df.iterrows():
        if not row['is_primary']:
            spectra_list.append([0.0 for i in range(2800)])
            continue
        smiles = str(row['canonical_smiles']).strip()
        hit_map[smiles] = 0
        if not smiles in extractor.smiles_set:
            print(f"Smiles not found in dictionary {smiles}")
            spectra_list.append([0.0 for i in range(2800)])
        else:
            hit_map[smiles] = extractor.smiles_set.count(smiles)
            spectra_list.append(extractor.get_merged_normalized_cutoff_binned_spectra(smiles))
    aggregator.df['massformer_spectra'] = spectra_list
    
    hits = 0
    multihits = 0
    for smiles, count in hit_map.items():
        if count == 0:
            print(f"Warning: No MASSFORMER data found for {smiles}")
            continue
        else:
            hits +=1
        
        if count > 1:
            print(f"Warning: Multiple MASSFORMER data found for {smiles}")
            multihits += 1
    print(f"Found MASSFORMER data for {hits} SMILES")
    print(f"Found multiple MASSFORMER data for {multihits} SMILES")



def write_descriptors_to_csv(output_csv: Union[str, Path], aggregator: TableAggregator):
    df = None
    if 'MolWt' in aggregator.df.columns:
        return

    for index, row in aggregator.df.iterrows():
        #if not row['is_primary']:
        #    continue
        mol = MolFromSmiles(row['canonical_smiles'])
        if mol is None:
            raise ValueError("Failed to parse SMILES string")

        # Generate 3D conformation for 3D descriptors
        #AllChem.EmbedMolecule(mol, randomSeed=42)
        #AllChem.MMFFOptimizeMolecule(mol)

        # Calculate all available descriptors
        descripts = Descriptors.CalcMolDescriptors(mol)
        # Get descriptor names (the keys of the descriptors dictionary)
        descriptor_names = [desc[0] for desc in Descriptors._descList]

        if df is None:
            # Create a DataFrame with descriptor values as a new row
            df = pd.DataFrame([descripts], columns=descriptor_names)
        else:
            # Append to existing DataFrame
            df = pd.concat([df, pd.DataFrame([descripts], columns=descriptor_names)], ignore_index=True)
        if (index + 1) % 100 == 0:
            print(f"Processed {index + 1} rows.")
    #print(df.corr())
    correlation_matrix = df.corr()
    # Initialize a dictionary to count high correlations for each column
    high_correlation_counts = {col: 0 for col in correlation_matrix.columns}
    # Iterate through the correlation matrix
    for row in correlation_matrix.index:
        for col in correlation_matrix.columns:
            if row != col and correlation_matrix.loc[row, col] > float('0.7'):
                #print(f"High correlation between {row} and {col}: {correlation_matrix.loc[row, col]}")
                high_correlation_counts[row] += 1
                high_correlation_counts[col] += 1
    
    # Sort columns by the count of high correlations
    sorted_columns = sorted(high_correlation_counts.items(), key=lambda item: item[1], reverse=True)
    
    # Print the sorted column names
    print("Columns sorted by the number of high correlations:")
    for col, count in sorted_columns:
        print(f"{col}: {count} high correlations")

    aggregator.df = pd.concat([aggregator.df, df], axis=1)

#
# Main - comment previous stages when re-running
# 
aggregator = TableAggregator(csv_path)
#write_descriptors_to_csv(output_path, aggregator)
#process_cfmid_predictions(f"{rootdir}/GNPSnew/GNPSnew_cfmid_chunked_files", aggregator)
#process_rassp_predictions(f"{rootdir}/GNPSnew/GNPSnew_rassp_chunked_files", aggregator)
#process_scarf_predictions(f"{rootdir}/GNPSnew/GNPSnew_scarf1_chunked_files", aggregator)
#process_massformer_predictions(f"{rootdir}/GNPSnew/GNPSnew_massformer_chunked_files", aggregator)
#process_iceberg_predictions(f"{rootdir}/GNPSnew/GNPSnew_iceberg1_chunked_files", aggregator)
#process_referense_spectra(f"{rootdir}/GNPSnew/GNPSnew_reference_spectra/ALL_GNPS_NO_PROPOGATED_filtered_since_20240101_noMissing_noAmbig_lkup2.json", aggregator)
#aggregator.save_processed(output_path)
#aggregator.check_binned_spect_size(2800,['iceberg_spectra', 'iceberg_probability', 'scarf_spectra', 'cfmid_spectra', 'rassp_spectra', 'massformer_spectra', 'reference_spectra'])
#aggregator.split_columns_to_rows(['iceberg_spectra', 'iceberg_probability', 'scarf_spectra', 'cfmid_spectra', 'rassp_spectra', 'massformer_spectra', 'reference_spectra'], output_path_split)
#aggregator.save_processed(output_path)
aggregator.mod_ms_intensity(1.0,['iceberg_spectra', 'iceberg_probability', 'scarf_spectra', 'cfmid_spectra', 'rassp_spectra', 'massformer_spectra', 'reference_spectra'])
aggregator.save_processed(output_path_mod1)
#aggregator.mod_ms_intensity(0.33,['iceberg_spectra', 'iceberg_probability', 'scarf_spectra', 'cfmid_spectra', 'rassp_spectra', 'massformer_spectra', 'reference_spectra'])
#aggregator.save_processed(output_path_mod2)
