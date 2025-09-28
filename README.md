# IsoFrag for Phytol
# Phytol Orbitrap Isotopologue Processing Pipeline

This repository contains a minimal reproducible pipeline for processing Orbitrap MS2 isotopologue data of phytol fragments.  
It implements the workflow described in *DOI*.

## Contents
- `phytol_processing.Rmd`: R Markdown script reproducing the key steps:
  - Import IsoX `.isox` files
  - Filter isotopologues
  - Calculate ion counts (Nio) and fractional abundances
  - Apply theoretical corrections for missing ²H isotopologues
  - Summarize corrected isotopologue ratios

- `data/`: Example input files (replace with your own IsoX files).
- `scripts/`: Custom helper functions (optional).
- `results/`: Example output files.

## Requirements
- R ≥ 4.2
- R packages:
  - `dplyr`, `ggplot2`, `tibble`, `isoorbi`

## Usage
1. Clone this repository:
   ```bash
   git clone https://github.com/YOURNAME/phytol-orbitrap.git
