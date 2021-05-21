#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 20:50:31 2021

@author: Manu Tamminen
"""

import logging
import pandas as pd
from collections import defaultdict
from itertools import groupby


def read_fasta(fasta, remove_commas=False):
    """ Read a fasta file and yield pairs of fasta ids and DNA sequences.
    """
    with open(fasta) as f:
        grouped = groupby(f, lambda x: x[0] == ">")
        for cond, entry in grouped:
            if cond:
                fasta_id = next(entry)
                fasta_id = fasta_id.replace(">", "").strip()
                if remove_commas:
                    fasta_id = fasta_id.replace(",", "")
                _, seq_iter = next(grouped)
                seq = ''.join([line.strip() for line in seq_iter]).upper()
                yield ([fasta_id, seq])


def create_mock_seq(ribo, number):
    mock_otu = """
    d:Mock{0}_{1}(1.00),p:Mock{0}_{1}(1.00),c:Mock{0}_{1}(1.00),
    o:Mock{0}_{1}(1.00),f:Mock{0}_{1}(1.00),g:Mock{0}_{1}(1.00),
    s:Mock{0}_{1}(1.00)""".format(number, ribo)
    return mock_otu.replace("\n", "").replace(" ", "")


def prepare_tax_dict(tax_file):
    tax_dict = defaultdict(list)
    with open(tax_file) as tax:
        for ix, line in enumerate(tax):
            try:
                otu, tax, _ = line.split("\t")
                if "Mock" in tax:
                    mock_seq = tax.split(",")[0]
                    mock_n, ribo_type = mock_seq.split("_", 1)
                    mock_n = mock_n.split("ock")[1]
                    ribo_type = ribo_type.split("(")[0]
                    tax = create_mock_seq(ribo_type, mock_n)
                tax_dict[tax].append(otu)
            except ValueError:
                pass
    return(tax_dict)


def sort_and_sum_otus(otus):
    otu_acc = []
    otu_size_acc = 0
    for otu in otus:
        otu_id, size = otu.split(";")
        _, size = size.split("=")
        otu_acc.append([otu, int(size)])
        otu_size_acc += int(size)
    sorted_otu_acc = sorted(otu_acc, reverse=True, key=lambda x: x[1])
    return([otu_size_acc, sorted_otu_acc[0][0]])


def prepare_tax_abund(taxonomy_dictionary):
    tax_list = []
    for tax, otus in taxonomy_dictionary.items():
        size_sum, abund_otu = sort_and_sum_otus(otus)
        tax_list.append([tax, size_sum, abund_otu])
    return tax_list


def collapse_taxonomies(input_taxonomy, input_fasta, output_table):
    tax_dict = prepare_tax_dict(input_taxonomy)
    tax_lst = pd.DataFrame(prepare_tax_abund(tax_dict),
                                columns=['Taxonomy', 'Size', 'OTU'])
    seqs = pd.DataFrame(read_fasta(input_fasta),
                             columns=['OTU', 'Sequence'])
    otu_seqs = pd.merge(left=tax_lst, right=seqs, how="left")
    otu_seqs.to_csv(output_table, header=False, index=False, sep="\t")


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    logger = logging.getLogger(__name__)
    logger.info('Collapsing the bacterial taxonomies...')
    collapse_taxonomies(
            input_taxonomy=snakemake.input[0],
            input_fasta=snakemake.input[1],
            output_table=snakemake.output[0])
    logger.info('Collapsing the eukaryotic taxonomies...')
    collapse_taxonomies(
            input_taxonomy=snakemake.input[2],
            input_fasta=snakemake.input[3],
            output_table=snakemake.output[1])
