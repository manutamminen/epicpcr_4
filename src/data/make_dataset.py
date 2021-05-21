#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 11:15:03 2021

@author: Manu Tamminen
"""

import logging
from itertools import groupby


def read_fasta(fasta, remove_commas=False):
    """ Read a fasta file and yield pairs of fasta ids and DNA sequences.
    """
    with open(fasta) as f:
        grouped = groupby(f, lambda x: x[0] == ">")
        for cond, entry in grouped:
            if cond:
                fasta_id = next(entry)
                if remove_commas:
                    fasta_id = fasta_id.replace(",", "")
                _, seq_iter = next(grouped)
                seq = ''.join([line.strip() for line in seq_iter]).upper()
                yield ([fasta_id, seq])


def process_fasta(fasta, ribo_type):
    """ Split a sequence at a barcode and only keep those where the total
    sequence is > 100 bp and the barcode 20 bp.
    """
    seq_acc = []
    for ix, (seq_id, seq) in enumerate(read_fasta(fasta)):
        if len(seq) > 100 and len(seq.split("GATCATGACCCATTTGGAGAAGATG")) == 2:
            line_acc = []
            line_acc.append(fasta.rsplit("/", 1)[1].split("_filtered")[0])
            line_acc.append(fasta.rsplit("/", 1)[1].split("_filtered")[0] + "_" + str(ix))
            bc, ribo = seq.split("GATCATGACCCATTTGGAGAAGATG")
            if ribo_type == "18S":
                try:
                    _, ribo = ribo.split("TGGTGGTGCATGGCCGTTCTT")
                    ribo, _ = ribo.split("ATAACAGGTCTGTGATGCC")
                except ValueError:
                    continue
            elif ribo_type == "16S":
                try:
                    _, ribo = ribo.split("GCCGCGGTAAT")
                    ribo = ribo[2:]
                    other_end = ribo[-20:]
                    if "GTAGTCC" in other_end:
                        ribo, _ = ribo.rsplit("GTAGTCC", 1)
                        ribo = ribo[:-13]
                    else:
                        continue
                except ValueError:
                    continue
            if len(bc) == 20:
                line_acc.append(bc)
                line_acc.append(ribo)
                seq_acc.append(line_acc)
    return seq_acc


def prepare_ribosomal_output(output_fasta, output_bc, ribo_type, glob_pattern):
    """ Pool the sequences and barcodes and write them into separate files.
    """
    with open(output_fasta, "w") as sf, open(output_bc, "w") as bc:
        fasta_files = glob_pattern
        for file_name in fasta_files:
            for sample, seq_id, barcode, ribo in process_fasta(file_name, ribo_type):
                sf.write(">" + seq_id + "\n" + ribo + "\n")
                bc.write(sample + "," + seq_id + "," + barcode + "\n")


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    logger = logging.getLogger(__name__)

    bact_seqs = [seqname for seqname in snakemake.input if "16S" in seqname]
    euk_seqs = [seqname for seqname in snakemake.input if "18S" in seqname]

    logger.info('Merging bacterial sequences and extracting barcodes...')
    prepare_ribosomal_output(
            output_fasta=snakemake.output[0],
            output_bc=snakemake.output[1],
            ribo_type="16S",
            glob_pattern=bact_seqs)

    logger.info('Merging eukaryotic sequences and extracting barcodes...')
    prepare_ribosomal_output(
            output_fasta=snakemake.output[2],
            output_bc=snakemake.output[3],
            ribo_type="18S",
            glob_pattern=euk_seqs)
