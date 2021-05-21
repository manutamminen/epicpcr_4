#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 04 2021

@author: Manu Tamminen
"""

import logging


def extract_fasta(otu_seq_tbl, fasta_output):
    with open(otu_seq_tbl) as tbl, open(fasta_output, 'w') as fas:
        for line in tbl:
            if "d:Mock" not in line:
                tax_id, count, _, seq = line.split("\t")
                tax_id = tax_id.replace(":", "_").replace("(", "__")
                tax_id = tax_id.replace(")", "_").replace(",", "_")
                tax_id = tax_id.replace("[", "").replace("]", "")
                tax_id = tax_id.replace(" ", "_")
                fas.write(">{}___{}\n{}".format(tax_id, count, seq))


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    logger = logging.getLogger(__name__)

    logger.info('Extracting the bacterial taxonomic OTUs...')
    extract_fasta( snakemake.input[2], snakemake.output[0])

    logger.info('Extracting the eukaryotic taxonomic OTUs...')
    extract_fasta( snakemake.input[3], snakemake.output[1])

