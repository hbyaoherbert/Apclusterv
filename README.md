# Apclusterv: Refinement of Viral Genome Clustering with Affinity Propagation
Apclusterv is a novel clustering software for viral genomes. The input genomes can be either complete genomes or contigs from metagenomic assembly. The program is based on protein-protein alignment and written in python<br>
The current stable version is 1.2.5
## Dependencies:<br>
### Operating System:  Unix Platforms
   python>=3.8<br>
   pandas<br>
   numpy<br>
   networkx >= 2.8.4 <br>
   scipy >=1.8.1<br>
   scikit-learn >= 1.1.2<br>
   MCL <br>
   diamond >= 0.9.14 <br>
   prodigal >= 2.6.3<br>
   R>=3.6.1<br>
   
## Installation: <br>
   Suppose you are in a conda environment, you need to install MCL, prodigal (for ORF prediction),diamond (for alignment) and R(if not already installed, we just need stats library in r-base)
   ```bash
   conda install diamond -c bioconda 
   conda install mcl -c bioconda
   conda install prodigal -c bioconda
   conda install r-base 
   
   pip install apclusterv==1.2.5
   ```
## Getting Started:<br>
### option 1. start with contigs 
   step1. preduct ORFs from the DNA file with the following command:  <br>
   ```bash
   prepare  contig_dna_fasta 
   ```
   (contig_dna_fasta is the path to the dna sequences for clustering)
   
   step2. execute clustering with the following command:
   ```bash
   apclusterv -contig contig_dna_fasta 
   ```
### option 2. if you already have protein sequences from the contigs, you can run apclusterv by the proteins and a protein-contig map file. 
   An example of protein file and mapfile are data/experiment1.faa and data/experiment1.csv
   ```bash
   apclusterv -protein experiment1.faa -csv experiment1.csv
   ```
### Help message and parameter setting
   apclusterv -h 
## Results <br>
   The program will create tmp/ directory. The clustering result is tmp/cluster_result.i.r.csv
   (cluster_result.3.4.csv by default)
   Simulation profile used in the manuscript is in data/profile.csv
   RI and ARI for evalation script is data/eval.py

