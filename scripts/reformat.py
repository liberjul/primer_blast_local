#!/usr/bin/env python
import pandas as pd
import os

def _encode_fasta_header(assay_num, target, type_str):
    if "Forward" in type_str:
        dir = "fwd"
    else:
        dir = "rev"
    return F"{assay_num}|{target}|{dir}"

def _decode_fasta_header(header_str, line_num=0):
    spl = header_str.strip().split("|")
    attrib_list = ["assay_number", "target", "direction"]
    if len(spl) != len(attrib_list) or spl[2] not in ("fwd", "rev") or " " in header_str:
        if line_num > 0:
            raise ValueError(F"Please reformat the primer file as specified in the help.\nHeader on line {line_num}: {header_str} is misformatted.")
        else:
            raise ValueError(F"Please reformat the primer file as specified in the help.\nHeader: {header_str} is misformatted.")
    else:
        return {x : y for x,y in zip(attrib_list, spl)}

def _check_and_read_valid_FASTA(file, primers=False):
    if primers:
        primer_dict = {}
    with open(file, "r") as ifile:
        lines = ifile.readlines()
        headers, seq_count, seq = 0, 0, ""
        for i in lines:
            if i[0] == ">":
                if seq != "":
                    if primers:
                        primer_dict[header] = seq
                    seq_count += 1
                if primers:
                    header = i[1:].strip()
                    _ = _decode_fasta_header(i)
                headers += 1

                seq = ""
            else:
                seq += i.strip()
        if seq != "":
            if primers:
                primer_dict[header] = seq
            seq_count += 1
    if headers == seq_count and headers > 0 and not primers:
        return True
    elif headers == seq_count and headers > 0 and primers:
        return True, primer_dict
    elif primers:
        return False, primer_dict
    else:
        return False

        if lines[0][0] == ">" and len(lines) > 1:
            return "FASTA"
        else:
            return None

def _determine_primerfile_type(file):
    try:
        df = pd.read_excel(file, sheet_name=0)
        return "EXCEL"
    except ValueError:
        qual, _ = _check_and_read_valid_FASTA(file, primers=True)
        if qual:
            return "FASTA"
        else:
            return None

def _idt_to_fasta(file):
    '''
    Coverts to the output of IDT Primer Quest export to a BLASTable FASTA.

    file: file object of inputted primer file
    '''
    idt_df = pd.read_excel(file, sheet_name=0) # read first sheet in file
    expected_cols = "Type", "AssaySet", "Sequence"
    if len([x for x in expected_cols if x in idt_df.columns]) < len(expected_cols):
        raise ValueError(F"Please verify that the primers file: {file} is formatted as expected from IDT PrimerQuest.\nSee test file at https://github.com/liberjul/primer_blast_local/blob/main/test/IDT_PrimerQuest_Export.xls")
    buffer = ""
    primer_dict = {}
    for i in range(len(idt_df)):
        type_str = idt_df.Type[i] # fwd, rev, or product
        if type_str != "Product": # if a primer
            assay_num, target = idt_df.AssaySet[i].strip(")").split(" (") # get assay set and target name
            header = _encode_fasta_header(assay_num, target, type_str).replace(" ", "_")
            buffer += F">{header}\n{idt_df.Sequence[i]}\n"
            primer_dict[header] = idt_df.Sequence[i]
    return buffer.strip(), primer_dict
