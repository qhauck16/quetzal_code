#!/bin/bash
# Build figure 2's three expensive caches as a SLURM job.
#
# Each cache loads the 3.5 GB flashier object and peaks around 7.5 GB, so it runs
# as one SLURM job (the three caches in sequence). Add -A/-p for your cluster.
#
#   bash scripts/submit_caches.sh          # all three: b, c, g
#   bash scripts/submit_caches.sh c        # just one
#
# Each writes into cache/; compose.R reads only those. Run once (or
# whenever the upstream flashier RDS changes).

set -euo pipefail
cd "$(dirname "$0")/.."
FIG2="$PWD"
RSCRIPT=Rscript
mkdir -p "$FIG2/logs"

PANELS="${*:-c g b}"   # c first: it is the one panel-c needs soonest

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=fig2caches
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=05:45:00
#SBATCH --output=$FIG2/logs/caches_%j.out
#SBATCH --error=$FIG2/logs/caches_%j.out

cd "$FIG2/scripts"
for p in $PANELS; do
  echo "=========== panel \$p cache START \$(date) on \$(hostname) ==========="
  $RSCRIPT --vanilla build_cache_panel_\${p}.R
  echo "=========== panel \$p cache EXIT=\$? \$(date) ==========="
done
echo "=== all requested caches done \$(date) ==="
ls -la "$FIG2/cache/"
EOF

squeue -u "$USER" -o "%.10i %.14j %.9P %.8T %.10M %R"
