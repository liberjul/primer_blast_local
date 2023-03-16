# primer_blast_local
A program for screening primers against a local database of genomes/sequences.

```
usage: primer_blast_local.py [-h] [-g GENOMES] [-p PRIMERS] [-o OUT] [-m TM_THRESH] [-t N_THREADS] [--min_size MIN_SIZE] [--max_size MAX_SIZE] [--na NA] [-k POT] [--tris TRIS] [--mg MG] [--dntps DNTPS] [--saltcorr SALTCORR] [--no_blast] [--use_existing_db]

optional arguments:
  -h, --help            show this help message and exit
  -g GENOMES, --genomes GENOMES
                        (str) FASTA formatted file containing concatenated genome records. It may be helpful to include the genome name in each record header.
  -p PRIMERS, --primers PRIMERS
                        (str) A FASTA containing primer sequences OR an excel file output by IDT PrimerQuest. If a FASTA, records should be formatted without spaces like >Assay_Name(unique)|Target_Name|Direction(fwd|rev)
  -o OUT, --out OUT     (str) Path and prefix of outputs.
  -m TM_THRESH, --tm_thresh TM_THRESH
                        (float) Minimum melting temperature at which primers should be included. Default 45.0.
  -t N_THREADS, --n_threads N_THREADS
                        (int) Number of threads to use for BLASTN. Default 1.
  --min_size MIN_SIZE   (int) Minimum amplicon size to include in results. Default 20.
  --max_size MAX_SIZE   (int) Maximum amplicon size to include in results. Default 9999.
  --na NA               (float) Sodium concentration, in millimolar. Default 50.0.
  -k POT, --pot POT     (float) Potassium concentration, in millimolar. Default 0.0.
  --tris TRIS           (float) Tris concentration, in millimolar. Default 0.0.
  --mg MG               (float) Magnesium concentration, in millimolar. Default 0.0.
  --dntps DNTPS         (float) dNTP concentration, in millimolar. Default 0.0.
  --saltcorr SALTCORR   (int) Salt correction method. See https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction Default 5.
  --no_blast            If specified, don't rerun blast but just change parameters for matches.
  --use_existing_db     If specified, don't rebuild the database
  ```
