#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 09:34:01 2021

@author: Manu Tamminen
"""

import logging
import pandas as pd


def process_uc(input_file):
    with open(input_file) as uc:
        uc_acc = []
        for line in uc:
            line = line.split()
            if line[0] == "S":
                uc_acc.append([line[8], line[8]])
            elif line[0] == "H":
                uc_acc.append([line[9], line[8]])
    return pd.DataFrame(uc_acc, columns=('OTU', 'Read'))


def process_tax(input_file):
    with open(input_file) as tax:
        tax_acc = []
        for line in tax:
            try:
                otu, taxonomy, _ = line.split("\t")
            except ValueError:
                continue
            otu, _ = otu.split(";")
            tax_acc.append([otu, taxonomy])
    return pd.DataFrame(tax_acc, columns=('OTU', 'Taxonomy'))


def merge_bcs_to_otus(input_uc, input_tax, input_bc, output_file):
    uc_tbl = process_uc(input_uc)
    tax_tbl = process_tax(input_tax)
    bc_tbl = pd.read_csv(input_bc,
                         names=('Sample', 'Read', 'BC'))
    uc_bc_tbl = pd.merge(uc_tbl, bc_tbl, how="left", on="Read")
    uc_bc_tax_tbl = pd.merge(uc_bc_tbl, tax_tbl, how="left", on="OTU")
    count_tbl = uc_bc_tax_tbl[['Sample', 'BC', 'Taxonomy']].value_counts()
    count_tbl = count_tbl.reset_index(name="Count")
    count_tbl.to_csv(output_file, index=False, sep="\t")


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    logger = logging.getLogger(__name__)
    logger.info('Merging the bacterial barcodes with taxonomic information...')
    merge_bcs_to_otus(
            input_uc=snakemake.input[0],
            input_tax=snakemake.input[1],
            input_bc=snakemake.input[2],
            output_file=snakemake.output[0])
    logger.info('Merging the eukaryotic barcodes with taxonomic information...')
    merge_bcs_to_otus(
            input_uc=snakemake.input[3],
            input_tax=snakemake.input[4],
            input_bc=snakemake.input[5],
            output_file=snakemake.output[1])
