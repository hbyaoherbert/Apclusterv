import subprocess
import os



#util code adopted from vcontact2
def make_diamond_db(aa_fp, db_dir, cpu: int):

    diamond_db_bp = os.path.join(db_dir, os.path.basename(aa_fp).rsplit('.', 1)[0])
    make_diamond_cmd = ['diamond', 'makedb', '--threads', str(cpu), '--in', aa_fp, '-d', diamond_db_bp]


    res = subprocess.run(make_diamond_cmd, check=True, stdout=subprocess.PIPE)

    diamond_db_fp = diamond_db_bp + '.dmnd'

    return diamond_db_fp


def run_diamond(aa_fp, db_fp, cpu: int, evalue: float, alignments: int, diamond_out_fn):

    # More sensitive as an option?
    diamond_cmd = ['diamond', 'blastp', '--threads', str(cpu),
                   '--sensitive', '--evalue', str(evalue), '--max-target-seqs', str(alignments),
                   '-d', db_fp, '-q', aa_fp, '-o', diamond_out_fn]


    res = subprocess.run(diamond_cmd, check=True, stdout=subprocess.PIPE)


    return diamond_out_fn


def make_protein_clusters_mcl(blast_fp, out_p, inflation=2):
    """
    Args: 
        blast_fp (str): Path to blast results file
        inflation (float): MCL inflation value
        out_p (str): Output directory path
    Returns:
        str: fp for MCL clustering file
    """


    blast_fn = os.path.basename(blast_fp)
    abc_fn = '{}.abc'.format(blast_fn)
    abc_fp = os.path.join(out_p, abc_fn)
    subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(blast_fp, abc_fp), shell=True)



    mci_fn = '{}.mci'.format(blast_fn)
    mci_fp = os.path.join(out_p, mci_fn)
    mcxload_fn = '{}_mcxload.tab'.format(blast_fn)
    mcxload_fp = os.path.join(out_p, mcxload_fn)
    subprocess.check_call("mcxload -abc {0} --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o {1}"
                          " -write-tab {2}".format(abc_fp, mci_fp, mcxload_fp), shell=True)

    mcl_clstr_fn = "{0}_mcl{1}.clusters".format(blast_fn, int(inflation*10))
    mcl_clstr_fp = os.path.join(out_p, mcl_clstr_fn)

    subprocess.check_call("mcl {0} -I {1} -use-tab {2} -o {3}".format(
        mci_fp, inflation, mcxload_fp, mcl_clstr_fp), shell=True)

    return mcl_clstr_fp
