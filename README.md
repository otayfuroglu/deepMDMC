# MLP Besed Flexible Framework Adsorption Simulation with MD/GCMC

This repository contains code for performing adsorption simulations in metal-organic frameworks (MOFs) using both MLP-based rigid GCMC and hybrid MD/GCMC simulations that incorporate framework flexibility. The simulations leverage a machine-learned potential (MLP), specifically a equivariant interatomic potentials (NeqIP)

The code is actively under development. We warmly welcome contributions from the community — whether it's fixing issues, improving documentation, or extending functionality.

## MD/GCMC 

(a)The overall workflow of the simulation algorithm is illustrated in the following scheme, summarizing the key steps of the MLP-based GCMC and hybrid MD/GCMC approach, (b) isotermal and (c) isobaric adsoption results.

<br/><br/>
<p align="center">
  <img src="workflow_ads_2.png" alt="Hybrid MD/GCMC Workflow" width="900"/>
</p>
<br/><br/>

### Releted paper
Modeling CO2 Adsorption in Flexible MOFs with Open Metal Sites via Fragment-Based Neural Network Potentials, Omer Tayfuroglu andSeda Keskin.
https://doi.org/10.26434/chemrxiv-2025-c85xt

## Installation

Before running the simulation code, ensure the following Python packages are installed:

- `nequip`
- `pair_nequip`
- `ase`
- `molmod`
- `CoolProp`

You can install them using pip:

```bash
pip install numpy torch ase molmod CoolProp
```

Clone this repository using:
```bash
git clone https://github.com/yourusername/your-repo-name.git
cd your-repo-name
```


## How To Use

## Cite us
Modeling CO2 Adsorption in Flexible MOFs with Open Metal Sites via Fragment-Based Neural Network Potentials, Omer Tayfuroglu andSeda Keskin.
https://doi.org/10.26434/chemrxiv-2025-c85xt
