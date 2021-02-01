#!/usr/bin/bash

#SBATCH --cpus-per-task=32
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --partition=m
#SBATCH --mem=500G

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o genome_index.out
#SBATCH -e genome_index.err
#SBATCH --job-name bowtie2index

ml load bowtie2/2.3.5.1-foss-2018b
genome='/groups/tanaka/People/current/jiwang/Genomes/axolotl/AmexG_v6.DD.corrected.round2.chr_mtDNA.NCBI.AJ584639.1.fa'
out="/groups/tanaka/People/current/jiwang/Genomes/axolotl/Bowtie2/am6.DD.mtDNA/am6"

bowtie2-build -f --threads 32 $genome $out
