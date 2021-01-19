#!/usr/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=480
#SBATCH --mem=80G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o genome_index.out
#SBATCH -e genome_index.err
#SBATCH --job-name index_am6_bowtie2

module load bowtie2/2.3.4.2-foss-2018b
genome='/groups/tanaka/People/current/jiwang/Genomes/axolotl/AmexG_v6.DD.corrected.round2.chr.fa'
out="/groups/tanaka/People/current/jiwang/Genomes/axolotl/Bowtie2/am6.DD"

bowtie2-build $genome $out
