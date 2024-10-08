#!/usr/bin/env python

import argparse, time, os
from reformat import _determine_primerfile_type, _check_and_read_valid_FASTA, _idt_to_fasta
from run_parse_blastn import _call_makeblastdb, _call_blastn, _blast_to_dict, _evaluate_hit_loc, _pull_amp_seqs

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genomes", type=str, help="FASTA formatted file containing concatenated genome records. It may be helpful to include the genome name in each record header.")
parser.add_argument("-p", "--primers", type=str, help="A FASTA containing primer sequences OR an excel file output by IDT PrimerQuest. If a FASTA, records should be formatted without spaces like >Assay_Name(unique)|Target_Name|Direction(fwd|rev)")
parser.add_argument("-o", "--out", type=str, help="Path and prefix of outputs.")
parser.add_argument("-m", "--tm_thresh", type=float, default=45., help="Minimum melting temperature at which primers should be included. (default:45.)")
parser.add_argument("-e", "--evalue", type=float, default=10., help="Maximum e-value of BLAST hits to evaluate. (default:10)")
parser.add_argument("-t", "--n_threads", type=int, default=1, help="Number of threads to use for BLASTN. (default:1)")
parser.add_argument("--min_size", type=int, default=20, help="Minimum amplicon size to include in results. (default:20)")
parser.add_argument("--max_size", type=int, default=9999, help="Maximum amplicon size to include in results. (default:9999)")
parser.add_argument("--na", type=float, default=50., help="Sodium concentration, in millimolar. (default:50)")
parser.add_argument("-k", "--pot", type=float, default=0., help="Potassium concentration, in millimolar. (default:0)")
parser.add_argument("--tris", type=float, default=0., help="Tris concentration, in millimolar. (default:0)")
parser.add_argument("--mg", type=float, default=0., help="Magnesium concentration, in millimolar. (default:0)")
parser.add_argument("--dntps", type=float, default=0, help="dNTP concentration, in millimolar. (default:0)")
parser.add_argument("--saltcorr", type=int, default=5, help="Salt correction method. See https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction  (default: 5)")
parser.add_argument("--no_blast", action="store_true", help="If specified, don't rerun blast but just change parameters for matches.")
parser.add_argument("--use_existing_db", action="store_true", help="If specified, don't rebuild the database.")
parser.add_argument("--amp_seq", action="store_true", help="If specified, include the sequence of the amplicon.")
args = parser.parse_args()

log_file = F"primer_blast_local_{time.strftime('%Y-%m-%d_%H-%M-%S')}.stderr.log"

#Check the inputs
if not args.use_existing_db:
    print("Verifying genome file...")
    if not _check_and_read_valid_FASTA(args.genomes):
        raise ValueError("Please verify that the genomes file is in FASTA format.")
    print("Genome file verfied...")

print("Verifying primer file...")
primerfile_type = _determine_primerfile_type(args.primers)
if primerfile_type == None:
    raise ValueError("Please verify that the primer file type is a valid XLS or FASTA file.")
if primerfile_type == "EXCEL":
    buffer, primer_dict = _idt_to_fasta(args.primers)
    primer_fasta = os.path.splitext(args.primers)[0] + "__formatted.fasta"
    with open(primer_fasta, "w") as ofile:
        ofile.write(buffer)
else:
    primer_fasta = args.primers
    qual, primer_dict = _check_and_read_valid_FASTA(primer_fasta, primers = True)
    if not qual:
        raise ValueError("Please reformat the primer file as specified in the help.")
print("Primer file verified...")

if not os.path.isdir(os.path.dirname(args.out)) and os.path.dirname(args.out) != "":
    os.mkdir(os.path.dirname(args.out))

# call commands
print("Running BLASTn...")
blast_out = args.out + "__blastn.out"
if not args.no_blast:
    if not args.use_existing_db:
        blast_db = _call_makeblastdb(args.genomes, log_file)
    else:
        blast_db = os.path.splitext(args.genomes)[0] + "__BLAST"
    _call_blastn(primer_fasta, blast_db, args.n_threads, args.evalue, log_file, blast_out)
print("BLASTn finished...")
print("Parsing results...")
blast_d = _blast_to_dict(blast_out)
buffer_passing, buffer_all = _evaluate_hit_loc(blast_d, primer_dict, tm_thresh = args.tm_thresh, size_max=args.max_size, size_min=args.min_size, Na=args.na, K=args.pot, Tris=args.tris, Mg=args.mg, dNTPs=args.dntps, saltcorr=args.saltcorr)
if args.amp_seq:
    buffer_passing = _pull_amp_seqs(buffer_passing, args.genomes, log_file, Na=args.na, K=args.pot, Tris=args.tris, Mg=args.mg, dNTPs=args.dntps, saltcorr=args.saltcorr)
with open(args.out + "__results.pass.csv", "w") as ofile:
    ofile.write(buffer_passing)
with open(args.out + "__results.all.csv", "w") as ofile:
    ofile.write(buffer_all)
print("Finished parsing results...")
