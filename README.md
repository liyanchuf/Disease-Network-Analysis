# Disease-Network-Analysis
This repository provides code for computing the Salton Cosine Index (SCI) and constructing disease comorbidity networks, as used in the study:
**"Disease network analysis to reveal comorbidity patterns in hospitalized patients with COPD using large-scale administrative health data."**
In addition, interactive HTML versions are provided for both male and female comorbidity networks. 


## Contents
- `SCI_construct.py`: Python script for calculating cosine similarity and building the network.
- `Figure 2a_male comorbidity network.html`: Interactive visualization of the male comorbidity network (Figure 2a).
- `Figure 2a_female comorbidity network.html`: Interactive visualization of the female comorbidity network (Figure 2a).
- No external files required; a demo dataset is included in the `__main__` section.

## Usage
```bash
python SCI_construct.py

## Requirements
- Python ≥ 3.6  
- Libraries:
  - `numpy`
  - `pandas`
  - `math`
