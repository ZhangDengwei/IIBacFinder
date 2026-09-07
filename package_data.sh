#!/usr/bin/env bash
set -euo pipefail

# Run this script on the author's machine to build the data-only archive.
# Upload the resulting IIBacFinder_data.tar.gz to Zenodo and publish its SHA-256.
OUT="IIBacFinder_data.tar.gz"

tar -czf "$OUT" \
  domains \
  models \
  AMP_database \
  scripts/cleavage_pred/training_data \
  test_fasta

sha256sum "$OUT" > "$OUT.sha256"
echo "Created $OUT"
echo "SHA-256: $(awk '{print $1}' "$OUT.sha256")"
