#!/usr/bin/env bash
set -euo pipefail

# Default data-only Zenodo archive for IIBacFinder.
IIBACFINDER_DATA_URL="${IIBACFINDER_DATA_URL:-https://zenodo.org/records/22211103/files/IIBacFinder_data.tar.gz?download=1}"
IIBACFINDER_DATA_SHA256_URL="${IIBACFINDER_DATA_SHA256_URL:-https://zenodo.org/records/22211103/files/IIBacFinder_data.tar.gz.sha256?download=1}"

OUT="IIBacFinder_data.tar.gz"
SHA256_FILE="IIBacFinder_data.tar.gz.sha256"

download_file() {
  local url="$1"
  local output="$2"
  if command -v wget >/dev/null 2>&1; then
    wget -O "$output" "$url"
  else
    curl -fL "$url" -o "$output"
  fi
}

download_file "$IIBACFINDER_DATA_URL" "$OUT"
download_file "$IIBACFINDER_DATA_SHA256_URL" "$SHA256_FILE"

sha256sum -c "$SHA256_FILE"
tar -xzf "$OUT"
echo "IIBacFinder data resources downloaded, verified, and unpacked."
