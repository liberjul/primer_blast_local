#!/usr/bin/env python
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from reformat import _decode_fasta_header
import numpy as np
import subprocess

def _other_dir(dir):
    if dir == "fwd":
        return "rev"
    else:
        return "fwd"
def _count_matches(seq1, seq2, shift = 0):
    count = 0
    for i in range(min([len(seq1), len(seq2)])):
        if seq1[i+shift] == seq2[i]:
            count += 1
    return count


def _find_3prime_mms(pseq, aseq): # find 3' end mismatches between primer and aligned sequence
    match_ct_arr = np.zeros(len(pseq))
    for shift in range(len(pseq)):
        match_ct_arr[shift] = _count_matches(pseq, aseq, shift = shift)
    best_shift = np.argmax(match_ct_arr)
    return len(pseq)+shift - len(aseq)

def _check_primer_quals(hit1, hit2, fwd_seq, rev_seq, tm_thresh = 45., size_max=9999, size_min=20, Na=50, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5):
    if hit1["sseqid"] == hit2["sseqid"] and hit1["sstrand"] != hit2["sstrand"]: # Check opposite strand annealing
        end_diff = int(hit2["send"]) - int(hit1["send"]) # amplicon size
        if (hit1["sstrand"] == "+" and end_diff > 0) or (hit1["sstrand"] == "-" and end_diff < 0): # primers are convergent
            try:
                tm_fwd = mt.Tm_NN(hit1["qseq"], c_seq = Seq(hit1["sseq"]).reverse_complement(), Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=saltcorr) # Calculate Tm
            except ValueError:
                tm_fwd = 0
            try:
                tm_rev = mt.Tm_NN(hit2["qseq"], c_seq = Seq(hit2["sseq"]).reverse_complement(), Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=saltcorr)
            except ValueError:
                tm_rev = 0
            threeprime_end_mm_fwd = _find_3prime_mms(fwd_seq, hit1["sseq"]) # find 3' end mismatch for the forward primer
            threeprime_end_mm_rev = _find_3prime_mms(rev_seq, hit2["sseq"]) # find 3' end mismatch for the reverse primer
            if tm_fwd >= tm_thresh and tm_rev >= tm_thresh and size_min <= abs(end_diff) <= size_max and threeprime_end_mm_fwd == 0 and threeprime_end_mm_rev == 0:
                return True, tm_fwd, tm_rev, abs(end_diff)
            else:
                return False, 0., 0. , 0
        else:
            return False, 0., 0. , 0
    else:
        return False, 0., 0. , 0

def _call_makeblastdb(fasta, log_file):
    with open(log_file, "a") as log:
        db_basename = os.path.splitext(fasta)[0]
        subprocess.run(F" makeblastdb -in {fasta} -dbtype nucl -out {db_basename}__BLAST", check=True, shell=True, stderr=log")
    return F"{db_basename}__BLAST"

def _call_blastn(query, db, nt, log_file, out_file):
    with open(log_file, "a") as log:
        subprocess.run(F"blastn -query {query} -db {db} -num_threads {nt} -outfmt \"6 qseqid sseqid qstart qend sstart send evalue pident qcovs qseq sseq strand\" -max_target_seqs 30000 > {out_file}", check=True, shell=True, stderr=log")


def _blast_to_dict(file):
    '''
    Coverts the output BLASTN with FMT=6 to a dictionary of hits by query key.

    blastn -query <query> -db <db> -num_threads <nt> -outfmt "6 qseqid sseqid qstart qend sstart send evalue pident qcovs qseq sseq strand" -max_target_seqs 30000 > <output>

    file: file object of inputted primer file
    '''
    hit_keys = ["sseqid", "qstart", "qend", "sstart", "send", "evalue", "pident", "qcovs", "qseq", "sseq", "sstrand"]
    hit_dict = {}
    with open(file):
        line = ifile.readline()
        while line != "":
            spl = line.strip().split("\t")
            p_dict = _decode_fasta_header(spl[0])
            assay_num_target, dir = p_dict["assay_num"] + "|" + p_dict["target"], p_dict["direction"]
            if assay_num_target not in hit_dict:
                hit_dict[assay_num_target] = {dir : [{x : y for x,y in zip(hit_keys, spl[1:])}], _other_dir(dir) : []} # make dict for each hit line
            else:
                hit_dict[assay_num_target][dir].append({x : y for x,y in zip(hit_keys, spl[1:])})
            line = ifile.readline()
    return hit_dict

def _evaluate_hit_loc(hit_dict, primer_dict, tm_thresh = 45., size_max=9999, size_min=20, Na=50, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5):
    buffer = "Assay_name_and_target,Forward_primer_seq,Reverse_primer_seq,Tm_forward,Tm_reverse,amplicon_size\n"
    for assay_num_target in hit_dict:
        for x in hit_dict[assay_num_target]["fwd"]:
            for y in hit_dict[assay_num_target]["rev"]:
                fwd_seq, rev_seq = primer_dict[F"{assay_num_target}|fwd"], primer_dict[F"{assay_num_target}|rev"] # get primer sequences
                passing, tm_fwd, tm_rev, amp_size = _check_primer_quals(x, y, fwd_seq, rev_seq, tm_thresh=tm_thresh, size_max=size_max, size_min=size_min, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=saltcorr)
                if passing:
                    buffer += F"{assay_num_target},{fwd_seq},{rev_seq},{tm_fwd},{tm_rev},{amp_size}\n"
