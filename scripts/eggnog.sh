#!/bin/bash -l
#SBATCH -A snic2022-22-263
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -J eggnogmapper

module load bioinfo-tools
module load bioinfo-tools eggNOG-mapper/1.0.3

#python /sw/bioinfo/eggNOG-mapper/1.0.3/rackham/emapper.py -i spiro.faa --output spiro_NOG -m diamond --data_dir /sw/data/eggNOG/4.5.1
#python /sw/bioinfo/eggNOG-mapper/1.0.3/rackham/emapper.py -i wb.faa --output wb_NOG -m diamond --data_dir /sw/data/eggNOG/4.5.1
#python /sw/bioinfo/eggNOG-mapper/1.0.3/rackham/emapper.py -i muris.faa --output muris_NOG -m diamond --data_dir /sw/data/eggNOG/4.5.1
python /sw/bioinfo/eggNOG-mapper/1.0.3/rackham/emapper.py -i carpe_aa.fasta --output carpe_NOG -m diamond --data_dir /sw/data/eggNOG/4.5.1
python /sw/bioinfo/eggNOG-mapper/1.0.3/rackham/emapper.py -i trepo_aa.fasta --output trepo_NOG -m diamond --data_dir /sw/data/eggNOG/4.5.1
python /sw/bioinfo/eggNOG-mapper/1.0.3/rackham/emapper.py -i kbiala_aa.fasta --output kbiala_NOG -m diamond --data_dir /sw/data/eggNOG/4.5.1