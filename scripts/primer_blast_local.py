#!/usr/bin/env python

import argparse, time
from run_parse_blastn import _call_blastn, _blast_to_dict, _evaluate_hit_loc

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genomes", type=str, help="FASTA formatted file containing concatenated genome records. It may be helpful to include the genome name in each record header.")
parser.add_argument("-p", "--primers", type=str, help="A FASTA containing primer sequences OR an excel file output by IDT PrimerQuest. If a FASTA, records should be formatted >Assay_Name(unique)|Target_Name|Direction(fwd|rev)")
parser.add_argument("-m", "--tm_thresh", type=float, help="Minimum melting temperature at which primers should be included.")
parser.add_argument("--min_size", type=int, help="Minimum amplicon size to include in results.")
parser.add_argument("--max_size", type=int, help="Maximum amplicon size to include in results.")
parser.add_argument("--na", type=float, help="Sodium concentration, in millimolar.")
parser.add_argument("-k", type=float, help="Potassium concentration, in millimolar.")
parser.add_argument("--tris", type=float, help="Tris concentration, in millimolar.")
parser.add_argument("--mg", type=float, help="Magnesium concentration, in millimolar.")
parser.add_argument("--dNTPs", type=float, help="dNTP concentration, in millimolar.")
parser.add_argument("--saltcorr", type=int, help="Salt correction method. See https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction")
args = parser.parse_args()

log_file = F"primer_blast_local_{time.strftime('%Y-%m-%d_%H-%M-%S')}.stderr.log"

#Check the inputs
