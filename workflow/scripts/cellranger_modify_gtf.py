import re

input_gtf = snakemake.input.gtf[0]
output_gtf = snakemake.output[0]

regex_pattern = snakemake.params.regex_pattern
repl_pattern = snakemake.params.repl_pattern


def _extract_version_number(annoline):
    if re.search(regex_pattern, annoline):
        return re.sub(pattern=regex_pattern, repl=repl_pattern, string=annoline)
    else:
        return annoline


with open(input_gtf) as fh_gtf_in, \
        open(output_gtf, "w") as fh_gtf_out:
    modified = map(_extract_version_number, fh_gtf_in.readlines())
    fh_gtf_out.writelines(modified)
