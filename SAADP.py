#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     :  SAADP.py
@Time     :  2022/10/28
@Author   :  weylz
@Version  :  1.0.1
@Desc     :  Sequence Alignment Algorithm Dynamic Programming
'''

from pywebio import start_server
from pywebio.input import input_group, input, TEXT, FLOAT, textarea
from pywebio.output import put_markdown, put_text, put_collapse, put_grid, put_code
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def sadp():
    
    def check_seq(seq):
        for i in range(len(seq)):
            if seq[i] not in ["A", "T", "C", "G", "a", "t", "c", "g", "U", "u", "N", "n"]:
                return "Please input correct DNA sequence!"

    def check_penalty(num):
        if float(num) > 0:
            return "Gap penalties should be non-positive."

    data = input_group("Sequence Alignment - Dynamic Programming",[textarea("Please input the first sequence :", name = "seq_1", type = TEXT, validate = check_seq, required = True, value = 'TCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTG'), textarea("Please input the second sequence :", name = "seq_2", type = TEXT, validate = check_seq,  required = True, value = 'TTGCTAGCTTCCCCACCACCAAGACCTACTTTCCTCACTTTGATGTA'), input("Match score :", name = "num_match", type = FLOAT, value = 1), input("Mismatch score :", name = "num_mismatch", type = FLOAT, value = 0), input("Open GAP penalty :", name = "num_gap_open", type = FLOAT, validate = check_penalty, value = 0), input("Extend GAP penalty :", name = "num_gap_extend", type = FLOAT, validate = check_penalty, value = 0, help_text = "Gap opening penalty should be higher than gap extension penalty (or equal).")])

    alignments = pairwise2.align.globalms(data["seq_1"], data["seq_2"], data["num_match"], data["num_mismatch"], data["num_gap_open"], data["num_gap_extend"])
    put_markdown("# Sequence Alignment - Dynamic Programming Results")
    put_grid([[put_markdown("**Original Sequence 1**"), put_text(data["seq_1"])], [put_markdown("**Original Sequence 2**"), put_text(data["seq_2"])]])
    put_text(' ')

    put_grid([[put_markdown("**Parameters:**"), put_text("Match score"), put_text("Mismatch score"), put_text("Open GAP penalty"), put_text("Extend GAP penalty")], [put_text(" "), put_text(data["num_match"]), put_text(data["num_mismatch"]), put_text(data["num_gap_open"]), put_text(data["num_gap_extend"])]])
    
    put_markdown("## Alignment Type: *Global*")
    with put_collapse('Global Alignment Results', open = True):
        put_code(format_alignment(*alignments[0]))

    put_markdown("## Alignment Type: *Local*")
    with put_collapse('Local Alignment Results', open = False):
        for a in pairwise2.align.localms(data["seq_1"], data["seq_2"], data["num_match"], data["num_mismatch"], data["num_gap_open"], data["num_gap_extend"]):
            put_code(format_alignment(*a))

if __name__ == '__main__':
    start_server(sadp, port = 7698)
