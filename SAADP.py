#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     :  SAADP.py
@Time     :  2022/10/28
@Author   :  weylz
@Version  :  1.0.2
@Desc     :  Sequence Alignment Algorithm Dynamic Programming
'''

from pywebio import start_server
from pywebio.input import input_group, input, TEXT, FLOAT, textarea
from pywebio.output import put_markdown, put_text, put_collapse, put_grid, put_code
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def saa_dp():
    
    def check_seq(seq):
        for i in range(len(seq)):
            if seq[i] not in ["A", "T", "C", "G", "a", "t", "c", "g", "U", "u", "N", "n"]:
                return "Please input correct DNA sequence!"

    def check_mismatch_score(num):
        if float(num) > 0:
            return "Mismatch score should be non-positive."

    def check_open_gap_penalty(num):
        if float(num) < 0:
            return "Open gap penalty should be non-negative."

    def check_extend_gap_penalty(num):
        if float(num) < 0:
            return "Extend gap penalty should be non-negative."

    data = input_group("Sequence Alignment - Dynamic Programming",[textarea("Please input the first sequence :", name = "seq_1", type = TEXT, validate = check_seq, required = True, value = 'TCTGCCTTCCCCCAAAGACCGCACTTCGCTG'), textarea("Please input the second sequence :", name = "seq_2", type = TEXT, validate = check_seq,  required = True, value = 'TCGTCCTTCCGCCTACTTCCCGCCTCGACTG'), input("Match score :", name = "num_match", type = FLOAT, value = 1), input("Mismatch score :", name = "num_mismatch", type = FLOAT, validate = check_mismatch_score, value = 0), input("Open gap penalty :", name = "num_gap_open", type = FLOAT, validate = check_open_gap_penalty, value = 0), input("Extend gap penalty :", name = "num_gap_extend", type = FLOAT, validate = check_extend_gap_penalty, value = 0, help_text = "Notice : open gap penalty should be higher than extend gap penalty (or equal).")])

    alignments = pairwise2.align.globalms(data["seq_1"], data["seq_2"], data["num_match"], data["num_mismatch"], -data["num_gap_open"], -data["num_gap_extend"])
    put_markdown("# Sequence Alignment - Dynamic Programming Results")
    put_grid([[put_markdown("**Original Sequence 1**"), put_text(data["seq_1"])], [put_markdown("**Original Sequence 2**"), put_text(data["seq_2"])]])
    put_text(' ')

    put_grid([[put_markdown("**Parameters:**"), put_text("Match score"), put_text("Mismatch score"), put_text("Open GAP penalty"), put_text("Extend GAP penalty")], [put_text(" "), put_text(data["num_match"]), put_text(data["num_mismatch"]), put_text(data["num_gap_open"]), put_text(data["num_gap_extend"])]])
    
    put_markdown("## Alignment Type: *Global*")
    with put_collapse('Global Alignment Results', open = True):
        put_code(format_alignment(*alignments[0]))

    put_markdown("## Alignment Type: *Local*")
    with put_collapse('Local Alignment Results', open = False):
        for a in pairwise2.align.localms(data["seq_1"], data["seq_2"], data["num_match"], data["num_mismatch"], -data["num_gap_open"], -data["num_gap_extend"]):
            put_code(format_alignment(*a))


if __name__ == '__main__':
    start_server(saa_dp, port = 8080, host = '127.0.0.1')
