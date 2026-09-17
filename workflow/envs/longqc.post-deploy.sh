#!/usr/bin/env bash
# Snakemake runs this once, right after it creates the longqc conda env, with
# CONDA_PREFIX pointing at the new env.
#
# Why this is needed
# ------------------
# LongQC resolves its two helper binaries relative to its own source file:
#
#     longQC.py:100   path_minimap2 = <dir of longQC.py>/minimap2-coverage
#     longQC.py:283   LqMask(os.path.join(path_minimap2, "sdust"), ...)
#     longQC.py:441   LqExec(os.path.join(path_minimap2, "minimap2-coverage"))
#
# In a conda env, dirname(longQC.py) is $CONDA_PREFIX/bin, so LongQC looks for
# $CONDA_PREFIX/bin/minimap2-coverage/sdust. But bioconda's minimap2-coverage
# recipe installs the *binary* as the file $CONDA_PREFIX/bin/minimap2-coverage,
# so that path can never resolve, and neither recipe builds `sdust` at all.
#
# The failure is silent: lq_mask.py:110 submits sdust through
# pool.apply_async() and never retrieves the result, so the NotADirectoryError
# is swallowed. _sdust opens its output file before running the binary, so the
# per-chunk tmp_N.txt files are created empty, concatenated into an empty
# longqc_sdust.txt, and longQC.py:372 then dies with
# pandas.errors.EmptyDataError after ~30 min of work.
#
# Fix: unpack the pinned upstream source into the env and build both binaries in
# place, which is the layout LongQC hard-codes. Rule longqc runs this copy of
# longQC.py rather than the one bioconda puts in bin/.
#
# LongQC's sdust is NOT minimap2's sdust. Upstream patched it (sdust.c, guarded
# by -D_SDUST_MAIN) to print one line per read -- name, masked bases, length,
# masked fraction, mean Q, Q7 count -- which is exactly what longQC.py:373-374
# and :455-460 index by column. The stock sdust shipped by the minimap2 package
# prints 3-column BED intervals and must not be substituted for it.
set -euo pipefail

version="1.2.0c"
# sha256 of the release tarball, as pinned by the bioconda longqc and
# minimap2-coverage recipes for this same version.
sha256="522837f655379881102233c69f8881866ef3d65116a5be61453428e7e989b01e"
url="https://github.com/yfukasawa/LongQC/archive/refs/tags/${version}.tar.gz"
dest="${CONDA_PREFIX}/share/LongQC"

tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT

curl -fsSL -o "${tmp}/longqc.tar.gz" "${url}"
python -c "
import hashlib, sys
digest = hashlib.sha256(open(sys.argv[1], 'rb').read()).hexdigest()
if digest != sys.argv[2]:
    sys.exit(f'LongQC tarball sha256 mismatch: got {digest}, expected {sys.argv[2]}')
" "${tmp}/longqc.tar.gz" "${sha256}"

rm -rf "${dest}"
mkdir -p "${dest}"
tar -xzf "${tmp}/longqc.tar.gz" -C "${dest}" --strip-components=1

# Same make invocation as the bioconda minimap2-coverage recipe, so the build
# picks up the env's zlib. The default target is `all`, which is
# `minimap2-coverage sdust`.
make -C "${dest}/minimap2-coverage" -j"$(nproc)" \
    INCLUDES="-I${CONDA_PREFIX}/include" \
    CFLAGS="-g -Wall -O2 -Wc++-compat -L${CONDA_PREFIX}/lib"

test -x "${dest}/minimap2-coverage/sdust"
test -x "${dest}/minimap2-coverage/minimap2-coverage"
test -f "${dest}/longQC.py"
