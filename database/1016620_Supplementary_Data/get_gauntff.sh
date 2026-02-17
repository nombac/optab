#!/bin/bash
set -euo pipefail

BASE_URL="https://data.nublado.org/gauntff"
DIR="$(cd "$(dirname "$0")" && pwd)"

for f in gauntff.dat gauntff_freqint.dat gauntff_nonav.dat; do
  echo "Downloading ${f}..."
  curl -f -o "${DIR}/${f}" "${BASE_URL}/${f}"
done

echo "Done."
