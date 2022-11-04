import os
import tempfile
from pathlib import Path
from snakemake.script import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


input_vcf = Path(snakemake.input.vcf_file).absolute()
output_dir = Path(snakemake.output.dir_snp).parent.absolute()
fasta_dir = Path(snakemake.input.ref_genome).parent.absolute()

shell(
    f"mkdir -p {output_dir}; "
    f"cd {output_dir}; "
    "SNPsplit_genome_preparation "
    "--vcf_file {input_vcf} "
    "--strain {snakemake.params.strain} "
    f"--reference_genome {fasta_dir} {log}"
)
