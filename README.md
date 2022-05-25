# Prerequisites

Make sure `nextflow`, `conda` and `singularity` are in the PATH.

# Quick start (lsf)

1. `git pull` this repo
2. Go to interactive session (you will need access to imaging queue):
    bsub -q imaging -n1 -M6000 -R"select[mem>6000] rusage[mem=6000]" -Is $SHELL
3. Run `bash [path-to-this-repo]/lsf_run [path-to-yaml-config]`, in a folder where you have `write` permission
