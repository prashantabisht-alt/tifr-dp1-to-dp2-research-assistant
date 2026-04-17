#!/usr/bin/env bash
# lapack_check.sh — detect whether LAPACK+BLAS is callable from gfortran.
# Usage:  bash lapack_check.sh
#
# It tries (in order):
#   1. macOS Accelerate framework  (-framework Accelerate)
#   2. Netlib LAPACK + BLAS         (-llapack -lblas)
#   3. OpenBLAS (Homebrew/Linux)    (-lopenblas)
# and prints the first one that works so we know which link flag to use.
set -euo pipefail

echo "=== checking gfortran ==="
if ! command -v gfortran >/dev/null 2>&1; then
  echo "✗ gfortran not found."
  echo "  macOS:  brew install gcc    (gfortran ships inside gcc)"
  echo "  Linux:  sudo apt install gfortran"
  exit 1
fi
gfortran --version | head -1
echo

# Write a minimal LAPACK test: eigenvalues of diag(1,2) via dgeev.
tmp=$(mktemp -d)
cat > "$tmp/lapack_check.f90" << 'EOF'
program lapack_check
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
  real(dp) :: a(2,2), wr(2), wi(2), vl(2,2), vr(2,2), work(20)
  integer :: info
  a = reshape([1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [2,2])
  call dgeev("N", "N", 2, a, 2, wr, wi, vl, 2, vr, 2, work, 20, info)
  if (info == 0) then
     print '(A,2F8.4)', " eigenvalues = ", wr
     stop 0
  else
     stop 1
  end if
end program
EOF

try_link() {
  local label="$1"; shift
  if gfortran "$tmp/lapack_check.f90" "$@" -o "$tmp/lc" 2>"$tmp/err" && "$tmp/lc" 2>&1; then
    echo "LINKED_WITH: $*"
    echo "LABEL: $label"
    exit 0
  fi
}

echo "=== trying link options ==="
echo "-> attempt: macOS Accelerate"
try_link "macOS Accelerate"            -framework Accelerate
echo "   (not available on this system)"
echo "-> attempt: netlib -llapack -lblas"
try_link "netlib LAPACK+BLAS"          -llapack -lblas
echo "   (not available on this system)"
echo "-> attempt: -lopenblas"
try_link "OpenBLAS"                    -lopenblas
echo "   (not available on this system)"

echo
echo "✗ No LAPACK backend found. Install suggestions:"
echo "  macOS:  xcode-select --install         (Accelerate ships with the CLT)"
echo "          brew install openblas lapack   (alternative)"
echo "  Linux:  sudo apt install liblapack-dev libblas-dev"
echo "  HPC:    module load lapack   OR   module load intel"
exit 1
