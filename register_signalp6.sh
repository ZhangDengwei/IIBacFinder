#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "ERROR: Activate the IIBacFinder conda environment before running this script." >&2
  exit 1
fi

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 /path/to/signalp-6.0X.fast.tar.gz" >&2
  exit 1
fi

LIBSTDCXX="$CONDA_PREFIX/lib/libstdc++.so.6"
if [[ ! -f "$LIBSTDCXX" ]]; then
  echo "ERROR: $LIBSTDCXX was not found." >&2
  echo "Install the conda C++ runtime with:" >&2
  echo "  mamba install -y -p \"$CONDA_PREFIX\" -c conda-forge 'libstdcxx-ng>=13' 'libgcc-ng>=13'" >&2
  exit 1
fi

if command -v strings >/dev/null 2>&1 && ! strings "$LIBSTDCXX" | grep -q 'GLIBCXX_3.4.29'; then
  echo "ERROR: $LIBSTDCXX does not contain GLIBCXX_3.4.29." >&2
  echo "Update the conda C++ runtime with:" >&2
  echo "  mamba install -y -p \"$CONDA_PREFIX\" -c conda-forge 'libstdcxx-ng>=13' 'libgcc-ng>=13'" >&2
  exit 1
fi

# Force the loader to use the conda C++ runtime before the system /lib64 copy.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export LD_PRELOAD="$LIBSTDCXX${LD_PRELOAD:+:$LD_PRELOAD}"

exec signalp6-register "$1"
