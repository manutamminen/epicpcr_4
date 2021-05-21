# epicPCR experiment May 14 2021

The aim of this experiment was to find an optimal concentration gradient for synthetic and biological standards

The following samples are included:

| Sample name | Sample description |
| :---: | :---: |
| RhodoM100XDilBC10e5 | Only biological controls; 100X dilution of standards; 1e5 barcode molecules |
| WWRhodoM100XDilBC10e5 | Wastewater + biological controls; 100X dilution of standards; 1e5 barcode molecules |
| WWM100XDilBC10e5 | Only wastewater; 100X dilution of standards; 1e5 barcode molecules |
| RhodoM100XDilBC10e6 | Only biological controls; 100X dilution of standards; 1e6 barcode molecules |
| WWRhodoM100XDilBC10e6 | Wastewater + biological controls; 100X dilution of standards; 1e6 barcode molecules |
| WWM100XDilBC10e6 | Only wastewater; 100X dilution of standards; 1e6 barcode molecules |
| RhodoM10XDilBC10e5 | Only biological controls; 10X dilution of standards; 1e5 barcode molecules |
| WWRhodoM10XDilBC10e5 | Wastewater + biological controls; 10X dilution of standards; 1e5 barcode molecules |
| WWM10XDilBC10e5 | Only wastewater; 10X dilution of standards; 1e5 barcode molecules |
| RhodoM10XDilBC10e6 | Only biological controls; 100X dilution of standards; 1e6 barcode molecules |
| WWRhodoM10XDilBC10e6 | Wastewater + biological controls; 100X dilution of standards; 1e6 barcode molecules |
| WWM10XDilBC10e6 | Only wastewater; 100X dilution of standards; 1e6 barcode molecules |

Raw data available at https://zenodo.org/record/4766882

Lab protocols available at https://github.com/manutamminen/epicpcr_4/blob/main/docs/protocols.md

Summary of the results available at https://github.com/manutamminen/epicpcr_4/blob/main/docs/index.md


# Building

## Dependencies

- Snakemake
- VSEARCH
- SINA
- FastTree
- Tidyverse

Download the raw data into data/raw.

Start the processing pipeline by invoking `snakemake --cores all report`.



