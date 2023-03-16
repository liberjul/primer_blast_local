# primer_blast_local
A program for screening primers against a local database of genomes/sequences.

```
usage: primer_blast_local.py [-h] [-g GENOMES] [-p PRIMERS] [-o OUT] [-m TM_THRESH] [-t N_THREADS] [--min_size MIN_SIZE] [--max_size MAX_SIZE] [--na NA] [-k POT] [--tris TRIS] [--mg MG] [--dntps DNTPS] [--saltcorr SALTCORR] [--no_blast] [--use_existing_db]

optional arguments:
  -h, --help            show this help message and exit
  -g GENOMES, --genomes GENOMES
                        FASTA formatted file containing concatenated genome records. It may be helpful to include the genome name in each record header.
  -p PRIMERS, --primers PRIMERS
                        A FASTA containing primer sequences OR an excel file output by IDT PrimerQuest. If a FASTA, records should be formatted without spaces like
                        >Assay_Name(unique)|Target_Name|Direction(fwd|rev)
  -o OUT, --out OUT     Path and prefix of outputs.
  -m TM_THRESH, --tm_thresh TM_THRESH
                        Minimum melting temperature at which primers should be included.
  -t N_THREADS, --n_threads N_THREADS
                        Nuber of threads to use for BLASTN.
  --min_size MIN_SIZE   Minimum amplicon size to include in results.
  --max_size MAX_SIZE   Maximum amplicon size to include in results.
  --na NA               Sodium concentration, in millimolar.
  -k POT, --pot POT     Potassium concentration, in millimolar.
  --tris TRIS           Tris concentration, in millimolar.
  --mg MG               Magnesium concentration, in millimolar.
  --dntps DNTPS         dNTP concentration, in millimolar.
  --saltcorr SALTCORR   Salt correction method. See https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
  --no_blast            If specified, don't rerun blast but just change parameters for matches.
  --use_existing_db     If specified, don't rebuild the database
  ```
