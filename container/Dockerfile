# RNAcs container — linux/arm64 (Apple Silicon native).
#
# Most of the toolchain comes from conda-forge/bioconda. Three pieces do not and
# are handled explicitly: SISSIz is absent from bioconda and is compiled from
# source; alifoldz.pl ships in the RNAz source tarball but not in its conda
# package; ps2eps and epstopdf come from Debian. Refolding needs no binary at
# all -- it runs through the ViennaRNA Python module.
#
# Build context must contain vendor/ — run container/prepare-context.sh first.

# Platform is set by the build command (--platform linux/arm64); pinning it
# here as a constant makes buildx warn.
FROM mambaorg/micromamba:1.5.8

USER root

# ps2eps and epstopdf are needed for the per-candidate EPS/PDF plots, which the
# workflow produces regardless of do_png. build-essential and autotools are for
# SISSIz.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential autoconf automake \
        ps2eps texlive-font-utils \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY environment.container.yaml /tmp/environment.yaml
RUN micromamba install -y -n base -f /tmp/environment.yaml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH

# --- SISSIz: compiled from source, since bioconda has no package -------------
COPY vendor/SISSIz /tmp/SISSIz
# The source tree carries autotools auxiliary files as symlinks into whatever
# automake the host had; they dangle here. Drop them and regenerate with the
# container's own toolchain, so the build does not depend on the host at all.
RUN cd /tmp/SISSIz \
    && find . -name '*.o' -delete \
    && rm -f config.status config.cache config.log Makefile \
    && find . -type l ! -exec test -e {} \; -delete \
    && autoreconf -fi \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd / && rm -rf /tmp/SISSIz \
    && SISSIz --version

# --- Perl helper the conda package omits -------------------------------------
# It carries the shebang of whatever perl the host had (e.g. Homebrew's),
# which does not exist here. Repoint it at the container's own perl.
COPY vendor/perl/alifoldz.pl /usr/local/bin/
RUN sed -i "1s|^#!.*perl|#!/opt/conda/bin/perl|" /usr/local/bin/alifoldz.pl \
    && chmod +x /usr/local/bin/alifoldz.pl \
    && head -1 /usr/local/bin/alifoldz.pl \
    && perl -c /usr/local/bin/alifoldz.pl

# --- RNAConSnake itself ------------------------------------------------------
COPY rnaconsnake /opt/rnaconsnake
RUN cd /opt/rnaconsnake && pip install --no-deps -e . && RNAcs --version

WORKDIR /work
ENTRYPOINT ["RNAcs"]
CMD ["--help"]
