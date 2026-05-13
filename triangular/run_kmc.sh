#!/usr/bin/env bash
# Build and run the Fortran KMC. Run on your laptop where gfortran is available.

set -e
cd "$(dirname "$0")"

echo "Building..."
gfortran -O3 -fno-range-check -ffree-line-length-none kmc_triangular_jmvr.f90 -o kmc_triangular

echo "Running..."
./kmc_triangular | tee kmc_run.log

echo ""
echo "Now compare with corrected exact theory:"
echo "    python3 fig11_final_hex.py"
echo ""
echo "For the PI-friendly gnuplot version:"
echo "    python3 export_fig11_final_hex_gnuplot_data.py"
echo "    gnuplot fig11_final_hex.gnu"
