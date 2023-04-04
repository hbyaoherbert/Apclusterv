# Apclusterv: Clustering viral genomes with Affinity Propagation
This software works in Python3.<br>
Dependencies:<br>
   Pandas<br>
   Numpy<br>
   Networkx<br>
   Scipy<br>
   vConTACT2 (bioconda channel)<br>
   cluster_one (http://www.paccanarolab.org/cluster-one/)<br>
   
Usage:<br>
   apcluster.sh input_gene_fasta gene_to_contig_map out_dir<br>
   input_gene_fasta is the protein sequence file output from Prodigal gene prediciton for viral contigs. An example gene_contig_map can be found in data/<br>
Parameter configuration:<br>
   configure parameters in config file.<br>
