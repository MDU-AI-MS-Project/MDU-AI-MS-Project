# msp_utils.py
from dataclasses import dataclass
from typing import List, Dict, Optional, Set, Any
from collections import defaultdict
import json
import csv
from ast import literal_eval
from rdkit.Chem import MolFromSmiles, MolToSmiles, MolFromInchi, Descriptors
from table_aggregator import TableAggregator

@dataclass
class MassSpectrum:
    metadata: Dict[str, str]
    peaks: List[tuple[float, float]]
    
    @property
    def smiles(self) -> Optional[str]:
        return self.metadata.get('SMILES', None)
    
    @property
    def name(self) -> str:
        return self.metadata.get('Name', 'Unknown')
    
    @property
    def ion_mode(self) -> str:
        return self.metadata.get('Ion_mode', 'Unknown')
    
    @property
    def collision_energy(self) -> str:
        return self.metadata.get('Collision_energy', 'Unknown')
    
    def get_measurement_conditions(self) -> str:
        """Get a string describing the measurement conditions"""
        conditions = []
        if self.ion_mode != 'Unknown':
            conditions.append(f"Ion mode: {self.ion_mode}")
        if self.collision_energy != 'Unknown':
            conditions.append(f"CE: {self.collision_energy}")
        if 'Instrument_type' in self.metadata:
            conditions.append(f"Instrument: {self.metadata['Instrument_type']}")
        return '; '.join(conditions)
    
    def get_peaks_fingerprint(self) -> str:
        """Create a fingerprint of the peaks list for comparison"""
        sorted_peaks = sorted(
            [(round(mz, 4), round(intensity, 4)) 
             for mz, intensity in self.peaks]
        )
        return ';'.join(f"{mz},{intensity}" for mz, intensity in sorted_peaks)

class SpectraLookup:
    def __init__(self, spectra: List[MassSpectrum]):
        """Initialize with path to MSP file."""
        self.spectra = spectra
        self.smiles_set = list({spectrum.smiles for spectrum in self.spectra if spectrum.smiles is not None})
        print(f"Loaded {len(self.spectra)} total spectra")

    def get_matching_spectra(self, smiles: str = None, name: str = None) -> List[MassSpectrum]:
        """Get spectra matching either SMILES or name."""
        if smiles is None and name is None:
            raise ValueError("Must provide either SMILES or name")
        
        # Get matching spectra
        matching_spectra = []
        for spectrum in self.spectra:
            #print (f'@{spectrum.smiles}@-\n@{smiles}@')
            if smiles and spectrum.smiles == smiles:
                matching_spectra.append(spectrum)
            elif name and spectrum.name == name:
                matching_spectra.append(spectrum)
        
        return matching_spectra
    
    def get_merged_normalized_cutoff_binned_spectra(self, smiles: str = None, name: str = None) -> List[tuple[int, float]]:
        """Get merged, normalized, cutoff, binned spectra matching either SMILES or name.
        Returns:
        List of (bin_index, normalized_intensity) tuples
        """

        if smiles is None and name is None:
            raise ValueError("Must provide either SMILES or name")
        
        matching_spectra = self.get_matching_spectra(smiles=smiles, name=name)
        if not matching_spectra:
            return []
        #print(matching_spectra)
        # Calculate total number of bins
        min_mz = 100
        max_mz = 1500
        bin_width = 0.5
        num_bins = int((max_mz - min_mz) / bin_width)
    
        # Initialize all bins with zeros
        bins = {i: 0.00 for i in range(num_bins)}
    
        # Process each spectrum
        for spectrum in matching_spectra:
            # Find max intensity for the spectrum
            max_intensity = max((intensity for mz, intensity in spectrum.peaks if mz >= min_mz and mz < max_mz), default=0.00)
            # Skip if max intensity is negligible or zero
            if abs(max_intensity) < 0.01:
                continue
            for mz, intensity in spectrum.peaks:
                # Skip if mz is out of range
                if mz < min_mz or mz >= max_mz:
                    continue
                # Calculate bin index
                bin_idx = int((mz - min_mz) / bin_width)
                # Update bin with normalized maximum intensity
                bins[bin_idx] = max(bins[bin_idx], round(intensity / max_intensity, 2))
                #bins[bin_idx] = 1.00 if (intensity / max_intensity) > 0.01
                ### bins[bin_idx] = max(bins[bin_idx], round(round(intensity / max_intensity, 2)) / 0.33) + 1

        # Normalize and filter if we have any non-zero intensities
        #max_intensity = max(bins.values())
        #if max_intensity > 0:
            # Normalize all intensities
            #for idx in bins:
                #bins[idx] = round(bins[idx] / max_intensity, 2)

        # Convert to sorted list of tuples
        return [(bins[idx]) for idx in range(num_bins)]



    def get_unique_spectra(self, smiles: str = None, name: str = None) -> List[MassSpectrum]:
        """Get unique spectra matching either SMILES or name."""
        if smiles is None and name is None:
            raise ValueError("Must provide either SMILES or name")
        
        # Get matching spectra
        matching_spectra = self.get_matching_spectra(smiles=smiles, name=name)        
        if not matching_spectra:
            return []
        
        # Group spectra by their peak fingerprints
        unique_spectra = {}
        for spectrum in matching_spectra:
            fingerprint = spectrum.get_peaks_fingerprint()
            if fingerprint not in unique_spectra:
                unique_spectra[fingerprint] = spectrum
        
        result = list(unique_spectra.values())
        
        # Print summary if there were matches
        total = len(matching_spectra)
        unique = len(result)
        print(f"\nFound {total} total spectra")
        print(f"Unique peak patterns: {unique}")
        
        # Group by conditions to show distribution
        conditions_summary = defaultdict(int)
        for spectrum in matching_spectra:
            key = (spectrum.ion_mode, spectrum.collision_energy)
            conditions_summary[key] += 1
        
        print("\nMeasurement conditions distribution:")
        for (ion_mode, ce), count in sorted(conditions_summary.items()):
            print(f"  Ion mode: {ion_mode}, Collision Energy: {ce} - {count} spectra")
        
        return result

def read_msp_file(file_path: str) -> List[MassSpectrum]:
    """Read a Massbank MSP file and return a list of MassSpectrum objects."""
    spectra = []
    current_metadata = {}
    current_peaks = []
    in_peaks_section = False
    
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('Name:') and current_metadata:
                if current_metadata and current_peaks:
                    print(current_metadata)
                    spectra.append(MassSpectrum(current_metadata.copy(), current_peaks.copy()))
                current_metadata = {}
                current_peaks = []
                in_peaks_section = False
            
            if in_peaks_section:
                try:
                    mz, intensity = map(float, line.split())
                    current_peaks.append((mz, intensity))
                except ValueError:
                    in_peaks_section = False
            
            if not in_peaks_section:
                if ': ' in line:
                    key, value = line.split(': ', 1)
                    current_metadata[key] = value
                    if key == 'Num Peaks':
                        in_peaks_section = True
        
        if current_metadata and current_peaks:
            spectra.append(MassSpectrum(current_metadata.copy(), current_peaks.copy()))
    
    return spectra

def read_cfmid_file(file_path: str, smiles: str) -> List[MassSpectrum]:
    """
    Read a CFM-ID predicted spectrum file and return a list of MassSpectrum objects.
    
    The method extracts spectra for different energy levels and creates a MassSpectrum 
    for each energy level.
    """
    spectra = []
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
        if not content:
            print(f"Warning: File {file_path} is empty. Skipping SMILES {smiles}")
            return spectra
        if smiles not in content:
            print(f"Error: SMILES {smiles} not found in {file_path}")
            exit(5)

    # Extract metadata
    metadata = {}
    metadata_patterns = {
        'SMILES': r'#SMILES=(.*)',
        'InChiKey': r'#InChiKey=(.*)', 
        'Formula': r'#Formula=(.*)',
        'Precursor_mass': r'#PMass=(.*)'
    }
    
    for key, pattern in metadata_patterns.items():
        import re
        match = re.search(pattern, content)
        if match:
            metadata[key] = match.group(1)
            if key == 'SMILES':
                metadata[key] = MolToSmiles(MolFromSmiles(metadata[key]))
    
    # Find energy level sections
    import re
    energy_sections = re.findall(r'(energy\d+)(.*?)(?=energy\d+|$)', content, re.DOTALL)
    
    for energy_level, section in energy_sections:
        # Create a copy of metadata for this spectrum
        spectrum_metadata = metadata.copy()
        spectrum_metadata['Collision_energy'] = energy_level.replace('energy', '') + '%'
        spectrum_metadata['Name'] = f"CFMID_{metadata.get('Formula', 'Unknown')}_{energy_level}"
        
        # Parse peaks
        peaks = []
        lines = section.strip().split('\n')
        for line in lines:
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        peaks.append((mz, intensity))
                    except ValueError:
                        continue
            else:
                # Skip the rest - seems end of the peak section and start of the fragment section
                break
        # Create and append spectrum if peaks exist
        if peaks:
            spectra.append(MassSpectrum(
                metadata=spectrum_metadata, 
                peaks=peaks
            ))
    
    return spectra

def read_reference_file(file_path: str, aggregator: TableAggregator) -> List[MassSpectrum]:

    """
    Read reference spectra from JSON file and convert to List[MassSpectrum]
    
    :param file_path: Path to the JSON file
    :return: List of MassSpectrum objects
    """
    # Read the JSON file
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            ref_data = json.load(f)
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return []
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in file {file_path}.")
        return []
    
    spectra = []
    count = 0
    skip_count = 0
    # Iterate through reference data
    for item in ref_data:
        count += 1
        # Skip not relevant adducts
        if not (item.get('Adduct') == 'M+H' or item.get('Adduct') == '[M+H]+'):
            print (f"Skipping adduct {item.get('Adduct')}")
            skip_count += 1
            continue
        # Get SMILES
        ref_smiles = item.get('Smiles', '')
        if not ref_smiles:
            print(f"Error: Missing SMILES in reference data")
            continue
        else:
            print (f"Processing reference spectrum number {count}: for SMILES {ref_smiles}")
        smiles = MolToSmiles(MolFromSmiles(ref_smiles))
        if smiles not in aggregator.df['canonical_smiles'].values:
            #print(f"Warning: SMILES {smiles} not found in the table aggregator")
            continue
        if 'C(=O)O[C@@H]1C(CO)=C[C@@H]2C(=O)[C@]3(C=C(C)[C@H](O)[C@@]13O)[C@H](C)CC1C2C1(C)C' in smiles:
            print (f"Found missing {smiles}: {ref_smiles}")
            exit(0)
        # Prepare metadata
        metadata = {
            'SMILES': smiles,
            'Name': 'ref_spectra'  # Default name
        }
        
        # Convert spectral data to list of (mz, intensity) tuples
        peaks_json = item.get('peaks_json', '')
        peaks = []
        try:
            peaks = json.loads(peaks_json)
            #print (type(peaks))
            #print (f"Peaks: {peaks}")
            #exit(0)
        except (ValueError, json.JSONDecodeError) as e:
            print(f"Error: Could not parse peaks for SMILES {smiles}: {e}")
            continue        
        if len(peaks) == 0:
            print(f"Error: Missing peaks in reference data")
            continue
        # Create MassSpectrum object
        spectrum = MassSpectrum(metadata, peaks)
        spectra.append(spectrum)
    print(f"Processed {count} reference spectra, skipped {skip_count}")
    return spectra

def read_rassp_file(file_path: str) -> List[MassSpectrum]:

    """
    Read RASSP format JSON file and convert to List[MassSpectrum]
    
    :param file_path: Path to the RASSP JSON file
    :return: List of MassSpectrum objects
    """
    # Read the JSON file
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            rassp_data = json.load(f)
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return []
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in file {file_path}.")
        return []
    
    spectra = []
    
    # Iterate through predictions
    for prediction in rassp_data.get('predictions', []):
        # Skip unsuccessful predictions
        if not prediction.get('success', False):
            continue
        # Get SMILES
        rassp_smiles = prediction.get('smiles', '')
        smiles = MolToSmiles(MolFromSmiles(rassp_smiles))
        # Prepare metadata
        metadata = {
            'SMILES': smiles,
            'InChI': prediction.get('inchi', ''),
            'Name': 'RASSP'  # Default name
        }
        
        # Convert spectral data to list of (mz, intensity) tuples
        peaks = prediction.get('spect', [])

        # Create MassSpectrum object
        spectrum = MassSpectrum(metadata, peaks)
        spectra.append(spectrum)
    
    return spectra

def read_massformer_file(file_path: str, aggregator: TableAggregator) -> List[MassSpectrum]:
    """
    Read Massformer CSV file and convert to List[MassSpectrum]
    
    :param file_path: Path to the Massformer CSV file
    :return: List of MassSpectrum objects
    """
    spectra = []
    print (file_path)
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            # Use csv reader to handle potential complexities
            csv_reader = csv.DictReader(f)
            
            for row in csv_reader:
                # Prepare metadata
                metadata = {
                    'Name': f"Spec_{row['spec_id']}",
                    'Mol_ID': row['mol_id'],
                    'Group_ID': row['group_id'],
                    'Precursor_Type': row['prec_type'],
                    'Precursor_MZ': row['prec_mz'],
                    'Normalized_Collision_Energy': row['nce'],
                    'Instrument_type': row['inst_type'],
                    'Fragmentation_Mode': row['frag_mode'],
                    'Spectrum_Type': row['spec_type'],
                    'Ion_mode': row['ion_mode']
                }
                # massformer mol_id starts from 1 while indexes in df start from 0
                index=int(row['mol_id'])-1
                if index < 0 or index >= len(aggregator.df):
                    print(f"Warning: Invalid mol_id {row['mol_id']} for spectrum {row['spec_id']}")
                    continue
                smiles = aggregator.df.loc[index, 'canonical_smiles']
                metadata['SMILES'] = smiles
                molecular_weight = Descriptors.MolWt(MolFromSmiles(smiles))
                # Check if molecular weight matches precursor MZ minus H+
                if not 0.0 < (float(row['prec_mz']) - float(molecular_weight)) <= 1.0:
                    print(f"Warning: Mismatch in molecular weight for spectrum {row['spec_id']}, smiles {smiles}: {molecular_weight} vs {row['prec_mz']}")
                    continue
                # Parse peaks using json.loads()
                try:
                    peaks = literal_eval(row['peaks'])
                except (ValueError, json.JSONDecodeError) as e:
                    print(f"Warning: Could not parse peaks for spectrum {row['spec_id']}: {e}")
                    continue
                
                # Create MassSpectrum object
                spectrum = MassSpectrum(metadata, peaks)
                spectra.append(spectrum)
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
    except csv.Error as e:
        print(f"CSV Error: {e}")
    
    return spectra

def read_scarf_file(file_path: str) -> List[MassSpectrum]:
    """
    Read SCARF format json file and convert to List[MassSpectrum]
    
    :param file_path: Path to the SCARF json file
    :return: List of MassSpectrum objects
    """
    # Read the JSON file
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            spectrum_data = json.load(f)
            #print (spectrum_data['output_tbl'])
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return []
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in file {file_path}.")
        return []
    
    # Extract pairs (mz_charge, inten) and (mz_charge, probability)
    mz_inten_pairs = []
    for mz, inten in zip(spectrum_data['output_tbl']['mono_mass'], spectrum_data['output_tbl']['rel_inten']):
        if inten > 0.001:
            mz_inten_pairs.append((mz, inten))
    
    # Prepare metadata
    metadata = { 
    'SMILES': MolToSmiles(MolFromSmiles(spectrum_data['smiles']))
    }
    
    # Create MassSpectrum object
    spectrum = MassSpectrum(metadata, sorted(mz_inten_pairs, reverse=True))
    return [spectrum]

def read_iceberg_file(file_path: str) -> List[MassSpectrum]:
    """
    Read Iceberg format json file and convert to List[MassSpectrum]
    
    :param file_path: Path to the Iceberg CSV file
    :return: List of MassSpectrum objects
    """
    # Read the JSON file
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            iceberg_data = json.load(f)
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return []
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in file {file_path}.")
        return []
    
    # Extract pairs (mz_charge, inten) and (mz_charge, probability)
    mz_inten_pairs = []
    mz_prob_pairs = []

    for frag_hash, frag_data in iceberg_data['frags'].items():
        mz_values = frag_data.get('mz_charge', [])
        intens_values = frag_data.get('intens', [])
        prob_gen = frag_data.get('prob_gen', 0)
        
        for mz, inten in zip(mz_values, intens_values):
            if inten > 0.01 and prob_gen > 0.2:
                mz_inten_pairs.append((mz, inten))
                mz_prob_pairs.append((mz, prob_gen))

    # Prepare metadata
    metadata = { 
    'SMILES': MolToSmiles(MolFromInchi(iceberg_data['root_inchi']))
    }
    
    # Create MassSpectrum object
    spectrum = MassSpectrum(metadata, sorted(mz_inten_pairs, reverse=True))
    probability_spectrum = MassSpectrum(metadata, sorted(mz_prob_pairs, reverse=True))
    return [spectrum], [probability_spectrum]

# Example usage
if __name__ == "__main__":
    # Initialize the lookup
    #msp_spectra = read_msp_file("/mnt/c/MDU/filtered_massbank_2024.msp")
    #extractor = SpectraLookup(msp_spectra)
    rassp_spectra = read_rassp_file("/mnt/c/MDU/GNPSnew/GNPSnew_rassp_chunked_files/GNPSnew_rassp.chunk-0000.json")
    extractor = SpectraLookup(rassp_spectra)
    # Example lookups
    #smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
    
    # Look up by SMILES
    for smiles in extractor.smiles_set:
        unique_spectra = extractor.get_unique_spectra(smiles=smiles)
    
        # Process unique spectra
        for i, spectrum in enumerate(unique_spectra, 1):
            print(f"\nSpectrum {i}:")
            print(f"Name: {spectrum.name}")
            print(f"Conditions: {spectrum.get_measurement_conditions()}")
            print(f"Number of peaks: {len(spectrum.peaks)}")
