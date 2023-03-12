#!/usr/bin/env python
import pandas as pd

def _encode_fasta_header(assay_num, target, type_str):
    if "Forward" in type_str:
        dir = "fwd"
    else:
        dir = "rev"
    return F"{assay_num}|{target}|{dir}"

def _decode_fasta_header(header_str):
    spl = header_str.split("|")
    atrrib_list = ["assay_number", "target", "direction"]
    return {x : y for x,y in zip(spl, attrib_list)}

def idt_to_fasta(file):
    '''
    Coverts to the output of IDT Primer Quest export to a BLASTable FASTA.

    file: file object of inputted primer file
    '''
    idt_df = pd.read_excel(file, sheet_name=0) # read first sheet in file
    expected_cols = "Type", "AssaySet", "Sequence"
    if len([x for x in expected_cols if x in idt_df.columns]) < len(expected_cols):
        raise ValueError(F"Please verify that the primers file: {file} is formatted as expected from IDT PrimerQuest. See test file at ")
    buffer = ""
    primer_dict = {}
    for i in range(len(idt_df)):
        type_str = idt_df.Type[i] # fwd, rev, or product
        if type_str != "Product": # if a primer
            assay_num, target = idt_df.AssaySet[i].strip(")").split(" (") # get assay set and target name
            header = _encode_fasta_header(assay_num, target, type_str)
            buffer += F">{header}\n{idt_df.Sequence[i]}\n"
            primer_dict[header] = idt_df.Sequence[i]
    return buffer, primer_dict
