#!/bin/bash
# Gather the pieces conda cannot supply into the build context.
#
# These are third-party sources installed on this machine; they are vendored
# into the build context rather than committed, so the image can be rebuilt
# without redistributing them.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
V="$HERE/vendor"

SISSIZ_SRC="${SISSIZ_SRC:-$HOME/C/SISSIz}"
ALIFOLDZ="${ALIFOLDZ:-$HOME/.local/share/RNAz/perl/alifoldz.pl}"

rm -rf "$V" && mkdir -p "$V/perl"

for f in "$ALIFOLDZ"; do
  [ -f "$f" ] || {
    echo "missing: $f" >&2
    echo "  alifoldz.pl ships in the RNAz source tarball (perl/)" >&2
    echo "  it is not in the conda package" >&2
    exit 1
  }
  cp "$f" "$V/perl/"
done
[ -d "$SISSIZ_SRC" ] || {
  echo "missing SISSIz source: $SISSIZ_SRC" >&2
  echo "clone it with: git clone https://github.com/mtw/SISSIz \"$SISSIZ_SRC\"" >&2
  exit 1
}
cp -R "$SISSIZ_SRC" "$V/SISSIz"
# Drop host build artefacts so the container configures cleanly.
find "$V/SISSIz" \( -name '*.o' -o -name 'config.status' -o -name 'config.log' \
     -o -name 'config.cache' -o -name 'autom4te.cache' \) -exec rm -rf {} + 2>/dev/null || true
# Autotools aux files are symlinks into the host's automake and would dangle in
# the image; the Dockerfile regenerates them with autoreconf.
find "$V/SISSIz" -type l ! -exec test -e {} \; -delete 2>/dev/null || true

# The package itself, without runs, venvs or git history.
rm -rf "$HERE/rnaconsnake" && mkdir -p "$HERE/rnaconsnake"
( cd "$REPO" && tar --exclude='.git' --exclude='.venv*' --exclude='*_benchmark_review' \
    --exclude='container' --exclude='dist' --exclude='__pycache__' \
    -cf - src snakefile config.yaml pyproject.toml setup.py RNAcs README.md LICENSE resources docs \
) | tar -xf - -C "$HERE/rnaconsnake"

echo "context ready:"
du -sh "$V/SISSIz" "$V/perl" "$HERE/rnaconsnake" | sed 's/^/  /'
