import os
import tempfile
from pathlib import Path
from snakemake.script import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

output_dir = Path(snakemake.output.dir_snp).parent

shell(
    f"cd {output_dir}; ",
    "SNPsplit_genome_preparation ",
    "--vcf_file ../{snakemake.input.vcf_file} ",
    "--strain {snakemake.params.strain} ",
    "--reference_genome ../{snakemake.input.ref_genome} {log}"
)
