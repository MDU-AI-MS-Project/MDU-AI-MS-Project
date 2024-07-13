
# MDU AI MS Project

There are different algorithms for MS/MS predicition but each one is installed, trained and run differently, accepting the input in a different way and producing the output in a different format which makes running different models and comparing the results difficult.

the MDU AI MS Project allows to easily run 5 different models to be run on a specific molecule in different energy level and compare the results.

## Description

The tool is a UI wrapper which runs 5 different trained mass spectra predcition models on a specific molecule,  at different energy levels, collects the resuts, filters, normalizes them and then compares them.

The UI allows entering a molecule in smiles format,  setting a threshold for the specta peaks and choosing which models to run and on whice energy level. 

After running the models, converting ,cleaning and normalizing the results they are presented in a table and a histogram for comparsion.

The tool currently supports 5 algorithms: SCARF, ICEBERG, MassFromer, RASSP, CFM-ID


## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)

## Prerequisites

1) a **Linux** operating system (we tested on Ubuntu 22.04) 

2) The system must have **Anaconda** installed 

3) The system must have **git** installed

4) The system must have **Docker** installed 


## Installation
### Project root directory

The project needs to be install under a specific *project root directory* which has to exist before the installation.
All the repos clones and installations needs to be done from that directory.
Currently, the *project root directory* has to be */home/greg*


### Download and install the algorithms 


#### 1. Download & install MassFormer:
    
In the *project root directory*, clone the [MassFormer](https://github.com/Roestlab/massformer) repo from github:

```bash
git clone https://github.com/Roestlab/massformer
```

Follow the *CPU Environment* section in the  [MassFormer installation instructions](https://github.com/Roestlab/massformer?tab=readme-ov-file#massformer) 


#### 2. Download & Install SCARF & ICEBERG
In the *project root directory*, clone the [SCARF & ICEBERG](https://github.com/samgoldman97/ms-pred) repo from github:

```bash
git clone https://github.com/samgoldman97/ms-pred md-pred-new
```

Follow the [SCARF & ICEBERG installation instructions](https://github.com/samgoldman97/ms-pred?tab=readme-ov-file#install--setup-) 
    

#### 3. Download & Install RASSP

In the *project root directory*, clone the [RASSP](https://github.com/thejonaslab/rassp-public) repo from github:

```bash
git clone https://github.com/thejonaslab/rassp-public
```

follow [RASPP Installation instructions](https://github.com/thejonaslab/rassp-public?tab=readme-ov-file#option-1-local-installation) 


#### 4. Download CFM-ID 

Download the CFM-ID Docker image:

```bash
docker pull wishartlab/cfmid
```

### Download and install MDU AI MS Project

Under the project root directory create a directory called MDU_outputs

Clone the current repository into a subdirectory of the *project root directory* named MDU_outputs.

```bash
git clone https://github.com/MDU-AI-MS-Project/MDU-AI-MS-Project MDU_outputs
```

Create an conda environemnt and install the dependencies
```bash
conda create -n mduaims numpy pandas matchms matplotlib rdkit
```


## Usage

### Run the tool 
Go into the *project root directory* ,  activate the conda environment and run the python code 

```bash
conda activate mduaims
cd /home/greg/MDU_outputs
python msms_gui.py
```
### Fill the details in the UI
After the UI pops up,  fill in the following details:
- Moleculde name :  in *smiles* format, for example: *CCCO*
- Threshold: as float, for example *0.001*
- Choose which of the 10 possible predictions you want to run by checking in the empty checkboxes 

![GUI](screenshots/gui.jpg)

After filling it the details , press *Plot Graph* to show the histograms

### Restuls:

![Histogram](screenshots/figure.png)
*Histogram*
 
![Histogram3D](screenshots/figure3d.png)
*3D Histogram*


## License

[Specify the license under which your project is distributed.]

This project is licensed under the [License Name] License - see the [LICENSE](LICENSE) file for details.


