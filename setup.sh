#!/bin/sh
# Toolchain for pic/ (em2) and ccp1d/ in a fresh Ubuntu container
# (Claude cloud session, 2 cores, 7 GB): gfortran + Open MPI + numpy/matplotlib.
# Takes about one minute. Safe to run again.
set -e
SUDO=""; [ "$(id -u)" -ne 0 ] && command -v sudo >/dev/null 2>&1 && SUDO=sudo
if ! command -v mpif90 >/dev/null 2>&1 || ! command -v gfortran >/dev/null 2>&1; then
  export DEBIAN_FRONTEND=noninteractive
  $SUDO apt-get update -qq 2>/dev/null || true      # a 403 from an unrelated mirror is harmless
  $SUDO apt-get install -y -qq gfortran libopenmpi-dev openmpi-bin
fi
python3 -c "import numpy, matplotlib" 2>/dev/null || \
  pip install --break-system-packages -q numpy matplotlib
echo "---- toolchain ----"
gfortran --version | head -1
mpirun --version 2>/dev/null | head -1
python3 -c "import numpy, matplotlib; print('numpy', numpy.__version__, 'matplotlib', matplotlib.__version__)"
echo "cores: $(nproc)   memory: $(free -g | awk '/Mem/{print $2}') GB"
echo "mpirun in this container needs: mpirun --allow-run-as-root -np N ./prog   (add --oversubscribe for N > cores)"
