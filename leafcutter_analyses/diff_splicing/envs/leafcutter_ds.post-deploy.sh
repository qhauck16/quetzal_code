#!/usr/bin/env bash
# Snakemake runs this once after building envs/leafcutter_ds.yml: installs leafcutter's R
# package (github.com/davidaknowles/leafcutter, the commit used for the manuscript).
set -euo pipefail
Rscript -e 'remotes::install_github("davidaknowles/leafcutter/leafcutter", ref = "2c9907ef66adf0bfb3092f0ceb6886ee5c046fbb", upgrade = "never")'
