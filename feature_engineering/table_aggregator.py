import pandas as pd
from typing import List, Dict, Optional, Union, Callable
from pathlib import Path
from datetime import datetime
import csv
import json
import ast
from rdkit.Chem import MolFromSmiles, MolToSmiles

class TableAggregator:
    def __init__(self, csv_path: Union[str, Path]):
        """
        Initialize the aggregator with a CSV file path.
        
        Args:
            csv_path: Path to CSV file containing molecule data
        """
        self.csv_path = Path(csv_path)
        self.df = None
        self.primary_mol_ids = {}  # Maps canonical smiles to primary mol_id number
        self.synonym_groups = {}     # Maps primary mol_id to all related mol_ids
        self.load_and_process()

    def canonize_smiles(self, smiles: str) -> str:
        """
        Convert chemical smiles to canonical form.
        
        Args:
            smiles: Chemical smiles string
        
        Returns:
            Canonized smiles string
        """
        return MolToSmiles(MolFromSmiles(smiles.strip()))
    
    def load_and_process(self) -> None:
        """
        Load CSV file and process the data:
        1. Keep all original rows
        2. Mark primary mol_ids for each canonical smiles
        3. Create mappings for synonym groups
        """
        # Load original data
        self.df = pd.read_csv(self.csv_path)
        
        # Ensure required columns exist
        required_cols = {'mol_id', 'smiles'}
        if not required_cols.issubset(self.df.columns):
            raise ValueError(f"CSV must contain columns: {required_cols}")
        
        # Check if processing is needed
        process_cols = {'canonical_smiles', 'is_primary', 'primary_mol_id'}
        if not process_cols.isdisjoint(self.df.columns):
            print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
            print("Loaded CSV already contains processed data")
            return

        # Create canonical forms
        self.df['canonical_smiles'] = self.df['smiles'].apply(
            self.canonize_smiles
        )
        
        # Identify primary mol_ids (first occurrence of each canonical smiles)
        primary_rows = self.df.groupby('canonical_smiles')['mol_id'].first()
        self.primary_mol_ids = primary_rows.to_dict()
        self.df['is_primary'] = False  # Initialize
        self.df['primary_mol_id'] = self.df['canonical_smiles'].map(self.primary_mol_ids)
        self.df['is_primary'] = self.df['mol_id'] == self.df['primary_mol_id']
        
        # Create synonym groups mapping
        for canon_smiles, group in self.df.groupby('canonical_smiles'):
            primary_mol_id = self.primary_mol_ids[canon_smiles]
            self.synonym_groups[primary_mol_id] = sorted(group['mol_id'].tolist())
           
    def add_properties_row(self, canon_smiles: str, properties_dict: Dict[str, any]) -> None:
        """
        Add a set of properties to a specific row in the DataFrame.
        
        Args:
            canon_smiles (str): The smiles as row key.
            properties_dict (Dict[str, any]): A dictionary where keys are column headers and values are the data to be added.
        """
        row_index = self.df[self.df['canonical_smiles'] == canon_smiles].index[0]
        for column, value in properties_dict.items():
            self.df.loc[row_index, column] = value

    def split_columns_to_rows(self, columns: List[str], output_csv_path: Union[str, Path]) -> None:
        """
        Split columns into rows, creating a new row for each value.
        
        Args:
            columns: List of column names to split
        """
        size = 2800
        output_csv_path = Path(output_csv_path)

        # Open the output CSV file in write mode and write the header
        with open(output_csv_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=self.df.columns.tolist() + ['mz', 'bin_index'])
            writer.writeheader()

            # Iterate over each row in the DataFrame
            for index, row in self.df.iterrows():
                print(f"Splitting row {index}")
                #Ensure all the lists of values in spectra columns have the same length/size
                evaluated_columns = {col: ast.literal_eval(row[col]) for col in columns}
                for col in columns:
                    column_size = len(evaluated_columns[col])
                    assert (column_size == size), f"The row {index} in column {col} must have length {size}, but has length {column_size}"
                
                # Iterate over the indices of the lists
                new_rows = []
                for i in range(size):
                    # Construct a new row with text columns and the current index of list columns
                    new_row = {col: row[col] for col in self.df.columns if col not in columns}
                    try:
                        new_row.update({col: evaluated_columns[col][i] for col in columns})
                    except IndexError:
                        print(f"Error! Index: {index} Bin: {i}")
                        print (row)
                        exit(1)
                    #new_row.update({col: row[col][i] for col in columns})

                    # Add the bin_index column
                    new_row['mz'] = 100 + (i / 2)
                    new_row['bin_index'] = i
                    # Check if at least one value in the lists is non-zero
                    if any(evaluated_columns[col][i] > 0 for col in columns):
                        new_rows.append(new_row)
                # Write accumulated rows for the current original row to the CSV file
                writer.writerows(new_rows)
        #result_df = pd.DataFrame(new_rows)
        #self.df = result_df

    def mod_ms_intensity(self, bin_width: float, columns: List[str]) -> None:
        """
        Modify the intensity values of the given columns.
        
        Args:
            columns: List of column names to modify
        """
        for col in columns:
            # Assert that all values in the column are floats
            assert all(isinstance(value, float) for value in self.df[col]), f"All values in column {col} must be floats"
            self.df[col] = self.df[col].apply(lambda x: (int(x // bin_width)+1 if x > 0.0 else 0))

    def check_binned_spect_size(self, size: int, columns: List[str]) -> None:
        """
        Check if the binned spectra have the correct size.
        
        Args:
            size: The expected size of the binned spectra
            columns: List of column names to check
        """
        for col in columns:
            count = 0
            for value in self.df[col]:
                count += 1
                print(value)
                assert len(ast.literal_eval(value)) == size, f"The row {count} in column {col} must have length {size}, but has length {len(value)}"

    def add_property(self, name: str, values_dict: Dict[int, any]) -> None:
        """
        Add a new property column, mapping values only to primary rows.
        
        Args:
            name: Name of the property column
            values_dict: Dictionary mapping mol_id numbers to values
        """
        # Initialize column with None
        self.df[name] = None
        
        # Map values only to primary rows
        for canon_smiles, primary_mol_id in self.primary_mol_ids.items():
            if primary_seq in values_dict:
                mask = self.df['mol_id'] == primary_seq
                self.df.loc[mask, name] = values_dict[primary_seq]
    
    def get_primary_mol_id(self, mol_id: int) -> int:
        """
        Get primary mol_id number for given mol_id.
        
        Args:
            mol_id: Original mol_id number
            
        Returns:
            Primary mol_id number for this molecule
        """
        row = self.df[self.df['mol_id'] == mol_id].iloc[0]
        return self.primary_mol_ids[row['canonical_smiles']]
    
    def get_synonyms(self, mol_id: int) -> List[int]:
        """
        Get all synonym mol_ids for given mol_id number.
        
        Args:
            mol_id: Original mol_id number
            
        Returns:
            List of all related mol_id numbers
        """
        primary_seq = self.get_primary_mol_id(mol_id)
        return self.synonym_groups[primary_seq]
    
    def print_synonym_summary(self) -> None:
        """Print a summary of synonym groups"""
        # Count synonyms per group
        synonym_counts = {
            primary: len(seqs) 
            for primary, seqs in self.synonym_groups.items()
        }
        
        print(f"Total smiless: {len(self.df)}")
        print(f"Unique molecules: {len(self.primary_mol_ids)}")
        print(f"Groups with synonyms: {len([c for c in synonym_counts.values() if c > 1])}")
        
        print("\nGroups with synonyms:")
        for primary_seq, seqs in self.synonym_groups.items():
            if len(seqs) > 1:
                smiles = self.df[self.df['mol_id'] == primary_seq]['smiles'].iloc[0]
                print(f"\nPrimary mol_id: {primary_seq} (smiles: {smiles})")
                print(f"Synonyms: {seqs}")
    
    def save_processed(self, output_path: Union[str, Path]) -> None:
        """
        Save the processed data to a new CSV file.
        """
        self.df.to_csv(Path(output_path), index=False)

    def get_primary_rows(self) -> pd.DataFrame:
        """
        Get only the primary rows with their data.
        
        Returns:
            DataFrame containing only primary rows
        """
        return self.df[self.df['is_primary']].copy()

if __name__ == '__main__':
    # Example usage
    csv_path = '/mnt/c/MDU/GNPSnew/GNPSnew_canonical.csv'
    output_path = '/mnt/c/MDU/GNPSnew/GNPSnew_test.csv'
    
    # Initialize the aggregator
    aggregator = TableAggregator(csv_path)
    
    # Add a new property
    #values = [1, 2, 3, 4]
    #aggregator.add_property('property_name', values)
    
    # Import property values from files
    #def extract_value(content: str) -> int:
    #    return int(content.strip())
    #aggregator.import_property_from_files('property_file', 'property_{}.out', extract_value)
    
    # Import property values from a dictionary
    #data_dict = {1: 10, 2: 20, 3: 30, 4: 40}
    #aggregator.import_property_from_dict('property_dict', data_dict, key_type='original')
    
    # Save the processed data
    aggregator.save_processed(output_path)
    
    # Print summary of redundant smiles groups
    #aggregator.print_synonym_summary()