
from pathlib import Path
from snakemake.script import shell


log = snakemake.log_fmt_shell(stdout=True, stderr=True)

input_fasta = Path(snakemake.input.fasta).absolute()
output_dir = Path(snakemake.output).parent.absolute()

shell(
    f"cd {output_dir} && "
    "cellranger mkref "
    "--genome={snakemake.params.genome} "
    f"--fasta={input_fasta} "
    "--genes={snakemake.input.annotation} "
    "--memgb={snakemake.params.mem} {log}"
)
