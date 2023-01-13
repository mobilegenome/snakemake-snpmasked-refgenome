
import tempfile
import numpy as np

from pybedtools import BedTool
from Bio import SeqIO, Seq
from Bio import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser


test_BED_string = """
chr1   2 3 C       T
chr1   9 10 G       A
chr1   11 12 G       A
chr2   2 3 G       T
chr2   3 4 A       G
chr2   6 7 A       G
"""

test_FASTA_string = """>chr1
GCATCAGGCGAG
>chr2
AAGAATAACAA
"""


tmpfn = tempfile.NamedTemporaryFile()
tmpfn = tmpfn.name

with open(tmpfn, "w") as fout:
    fout.write(test_FASTA_string)

def filter_indels(df_snv):
    """
    For the following data format used for SNV data from
    https://www.ebi.ac.uk/research/flicek/publications/fog09/
    this method filters out any InDels.

    InDels are defined by coordinate inervals larger than 1 and
    by non-existing bases in the reference genome indicated by "*" for
    an reference allele.

    chr10   7337993 7337994 *       -T/-T
    chr10   7337991 7337992 T       T
    chr10   7337991 7337992 *       -G/-G
    chr10   7337984 7337985 A       G

    :param df_snv:
    :return: df_snv
    """
    df_snv = df_snv[(df_snv.end - df_snv.start) == 1]
    df_snv = df_snv[df_snv.ref.isin(list("ACTG"))]  # ignore indels

    return df_snv



snvs = BedTool(test_BED_string, from_string=True)
df_snv = snvs.to_dataframe().rename(columns=
                                    {"name": "ref",
                                     "score": "alt"})
df_snv = filter_indels(df_snv)

fasta = SeqIO.index(tmpfn, format="fasta", )


fasta_in = tmpfn
fasta_out = "sequence.fa"

with open(fasta_in) as fh_in, \
        open(fasta_out, "w") as fh_out:
    fasta = SimpleFastaParser(fh_in)

    for seqname, seq in fasta:
        seq = np.array(Seq(seq))
        snvs = df_snv[df_snv.chrom == seqname]
        coords, ref_alleles, alt_alleles = snvs.start.to_list(),snvs.ref.to_list(),snvs.alt.to_list()
        for coord, ref, alt in zip(coords, ref_alleles, alt_alleles):
            if seq[coord] == ref:
                seq[coord] = alt
                print(f"Successfully replaced {ref} with {alt} at index {seqname}:{coord}")
            else:
                print(f"Warning: {ref} does not match {seq[elem]} at index  {seqname}:{coord}")

        seqr = SeqRecord.SeqRecord(Seq("".join(seq.tolist())),
                                   id=seqname,
                                   description="")
        SeqIO.write(seqr, fh_out, "fasta")
