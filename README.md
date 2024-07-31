# Apclusterv: Clustering viral genomes with Affinity Propagation
This software works in Python3.<br>
Dependencies:<br>
   Python>=3.7<br>
   Pandas<br>
   Numpy<br>
   Networkx<br>
   Scipy<br>

   Diamond <br>

   R<br>

Installation:
   conda install diamond -c bioconda <br>
   conda install r-base (if R is not already installed,we just need to run Rscript with stats package)<br>
   pip install apclusterv <br>
   
Usage:<br>

   step1. preduct ORFs from the DNA file with the following command:  <br>
   prepare contig_dna_fasta <br>
   (contig_dna_fasta is the path to the dna sequences for clustering)
   
   step2. execute clustering with the following command:
   apclusterv contig_dna_fasta <br

Results <br>
   The program will create tmp/ directory. The clustering result is tmp/result.csv

