#!/usr/bin/env bash
set -euo pipefail

# Set these variables before running the script.
# IIBACFINDER_DATA_URL should point to the Zenodo data-only archive.
# IIBACFINDER_DATA_SHA256 should be the SHA-256 checksum published with the archive.
IIBACFINDER_DATA_URL="${IIBACFINDER_DATA_URL:?https://zenodo.org/records/22211103/files/IIBacFinder_data.tar.gz?download=1}"
IIBACFINDER_DATA_SHA256="${IIBACFINDER_DATA_SHA256:?https://zenodo.org/records/22211103/files/IIBacFinder_data.tar.gz.sha256?download=1}"

OUT="IIBacFinder_data.tar.gz"

if command -v wget >/dev/null 2>&1; then
  wget -O "$OUT" "$IIBACFINDER_DATA_URL"
else
  curl -fL "$IIBACFINDER_DATA_URL" -o "$OUT"
fi

echo "$IIBACFINDER_DATA_SHA256  $OUT" | sha256sum -c -
tar -xzf "$OUT"
echo "IIBacFinder data resources downloaded, verified, and unpacked."
