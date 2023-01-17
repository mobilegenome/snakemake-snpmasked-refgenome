
import tempfile
import numpy as np
import pandas as pd
import textwrap

from pathlib import Path

from pybedtools import BedTool
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser

input_fasta = Path(snakemake.input.fasta)
input_bed = Path(snakemake.input.bed[0])

output_fasta = Path(snakemake.output.fasta)

logfile = open(snakemake.log[0], "w")

print(f"Input files: {input_fasta}, {input_bed}", file=logfile)
print(f"Output file: {output_fasta}", file=logfile)

def test_inject_snvs():
    test_BED_string = \
        textwrap.dedent("""\chr1   2 3 C       T
        chr1   9 10 G       A
        chr1   11 12 G       A
        chr2   2 3 G       T
        chr2   3 4 A       G
        chr2   6 7 A       G
        """)

    test_FASTA_string = \
        textwrap.dedent("""\
        >chr1
        GCATCAGGCGAG
        >chr2
        AAGAATAACAA
        """)


    tmpfn = tempfile.NamedTemporaryFile()
    tmpfn = tmpfn.name

    with open(tmpfn, "w") as fout:
        fout.write(test_FASTA_string)

    snvs = BedTool(test_BED_string, from_string=True)

    df_snv = snvs.to_dataframe().rename(columns=
                                        {"name": "ref",
                                         "score": "alt"})
    df_snv = filter_indels(df_snv)

    fasta_in = tmpfn
    fasta_out = "sequence.fa"

    inject_snvs(fasta_in, fasta_out, df_snv)


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


def inject_snvs(fasta_in: Path,
                fasta_out: Path,
                df_snv: pd.DataFrame):

    seqindex = SeqIO.index(fasta_in, "fasta")
    modified_sequences = []

    for seqname in seqindex:
        if not seqname in df_snv.chrom.drop_duplicates().to_list():
            print(f"Chrom {seqname} not found in Bedfile skipping...", file=logfile)
            continue

        seq = seqindex[seqname].seq
        seq = np.array(seq)

        snvs = df_snv[df_snv.chrom == seqname]
        coords, ref_alleles, alt_alleles = snvs.start.to_list(), snvs.ref.to_list(),snvs.alt.to_list()
        for coord, ref, alt in zip(coords, ref_alleles, alt_alleles):
            if seq[coord] == ref:
                seq[coord] = alt
                print(f"Successfully replaced {ref} with {alt} at index {seqname}:{coord}", file=logfile)
            else:
                print(f"Warning: {ref} does not match {seq[coord]} at index  {seqname}:{coord}", file=logfile)

        seqr = SeqRecord.SeqRecord(Seq.Seq("".join(seq.tolist())),
                                   id=seqname,
                                   description="")
        modified_sequences.append(seqr)

        SeqIO.write(modified_sequences, fasta_out, "fasta")


snvs = BedTool(input_bed)
df_snv = snvs.to_dataframe().rename(columns=
                                    {"name": "ref",
                                     "score": "alt"})

df_snv["chrom"] = df_snv.chrom.astype(str)  # convert chromosome column to strings

df_snv = filter_indels(df_snv)

print(df_snv.chrom.drop_duplicates(), file=logfile)
inject_snvs(input_fasta, output_fasta, df_snv)

logfile.close()