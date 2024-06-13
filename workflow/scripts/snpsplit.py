import os
import tempfile
from pathlib import Path
from snakemake.script import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


input_vcf = Path(snakemake.input.vcf_file).absolute()
output_dir = Path(snakemake.output.dir_snp).parent.absolute()
fasta_dir = Path(snakemake.input.ref_genome).parent.absolute()

print(input_vcf)
print(output_dir)
print(fasta_dir)

extra = ""
if snakemake.params.incorporate_snvs:
    extra = "--full_sequence"

shell(
    f"mkdir -p {output_dir}; "
    f"cd {output_dir}; "
    "SNPsplit_genome_preparation "
    "--vcf_file {input_vcf} "
    "--strain {snakemake.params.strain} "
    f"{extra} "
    f"--reference_genome {fasta_dir} {log}"
)
