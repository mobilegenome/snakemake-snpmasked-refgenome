import os
import tempfile
from pathlib import Path
from snakemake.script import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

output_dir = Path(snakemake.output.dir_snp).parent
fasta_dir = Path(f"../{snakemake.input.ref_genome}").parent

shell(
    f"mkdir -p {output_dir}; "
    f"cd {output_dir}; "
    "SNPsplit_genome_preparation "
    "--vcf_file ../../../../{snakemake.input.vcf_file} "
    "--strain {snakemake.params.strain} "
    f"--reference_genome {fasta_dir} {log}"
)
