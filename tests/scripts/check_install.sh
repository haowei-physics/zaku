#!/usr/bin/env bash
# Verify an installed Zaku tree from a fresh environment:
#   bash tests/scripts/check_install.sh <install-prefix>
# Sources <prefix>/setup.sh, then checks ROOT autoload, PyROOT import and the
# zaku command. Exit status is non-zero on any failure.
set -u
prefix="${1:?usage: check_install.sh <install-prefix>}"
if [ ! -f "$prefix/setup.sh" ]; then
  echo "FAIL: $prefix/setup.sh not found (run cmake --install first)"
  exit 1
fi
# shellcheck disable=SC1090
source "$prefix/setup.sh"
cd /
status=0

echo "--- ROOT (header from ROOT_INCLUDE_PATH, library autoloaded through the rootmap) ---"
out=$(root -l -b -q -e '#include "zaku/Zaku.h"' -e 'std::cout << "ROOT:" << zaku::buildInfo() << " nucleons=" << zaku::Detector::water(188.0).nNucleons() << std::endl;' 2>&1)
echo "$out" | tail -3
echo "$out" | grep -q "^ROOT:Zaku" || { echo "FAIL: ROOT autoload"; status=1; }

echo "--- PyROOT ---"
out=$(python3 -c 'import zaku; print("PY:%s %.6e %s" % (zaku.version(), zaku.Detector.water(188).nNucleons(), zaku.DATA_DIR))' 2>&1)
echo "$out" | tail -3
echo "$out" | grep -q "^PY:" || { echo "FAIL: python import"; status=1; }
echo "$out" | grep -q "1.132" || { echo "FAIL: nNucleons"; status=1; }

echo "--- zaku command ---"
if ! command -v zaku > /dev/null; then echo "FAIL: zaku not on PATH"; status=1; else zaku --version; fi

echo "--- layout ---"
for f in lib/libZaku.so lib/libZaku_rdict.pcm lib/libZaku.rootmap include/zaku/Zaku.h python/zaku/__init__.py bin/zaku; do
  if [ -e "$prefix/$f" ]; then echo "ok   $f"; else echo "FAIL: missing $f"; status=1; fi
done

[ $status -eq 0 ] && echo "check_install: OK" || echo "check_install: FAILED"
exit $status
