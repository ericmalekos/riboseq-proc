# riboseq-proc image.
#   /opt/env       STAR, cutadapt, bowtie2, samtools, pigz, bedtools, rclone (locks/core.linux-64.lock)
#   /opt/ribotish  Ribo-TISH (locks/ribotish.linux-64.lock)
#   /opt/umitools  UMI-tools (locks/umitools.linux-64.lock)
# Three environments because their Python pins differ. Every package is pinned by URL and md5 in the lockfiles.
# The tag is the content of VERSION; a published tag is never rebuilt.
FROM mambaorg/micromamba:1.5.8-jammy
LABEL org.opencontainers.image.source="https://github.com/ericmalekos/riboseq-proc" \
      org.opencontainers.image.description="Ribo-seq reads to uniquely mapped, ncRNA-filtered BAMs with a 3-nt periodicity check"
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates gawk procps \
 && rm -rf /var/lib/apt/lists/*
COPY locks/ /opt/riboseq-proc/locks/
RUN micromamba create -y -p /opt/env --file /opt/riboseq-proc/locks/core.linux-64.lock \
 && micromamba create -y -p /opt/ribotish --file /opt/riboseq-proc/locks/ribotish.linux-64.lock \
 && micromamba create -y -p /opt/umitools --file /opt/riboseq-proc/locks/umitools.linux-64.lock \
 && micromamba clean -a -y
COPY VERSION /opt/riboseq-proc/VERSION
COPY bin/ /opt/riboseq-proc/bin/
COPY refs/ /opt/riboseq-proc/refs/
RUN chmod 0755 /opt/riboseq-proc/bin/riboseq-proc /opt/riboseq-proc/bin/periodicity_gate.py /opt/riboseq-proc/refs/build_refs.sh \
 && mkdir -m 1777 /work
ENV PATH=/opt/riboseq-proc/bin:/opt/riboseq-proc/refs:/opt/env/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    PYTHONNOUSERSITE=1 \
    LC_ALL=C.UTF-8 \
    WORK=/work
# fails the build when a tool is missing or does not start
RUN riboseq-proc version && bedtools --version && python3 -I /opt/riboseq-proc/bin/periodicity_gate.py --help > /dev/null \
 && command -v flock setsid > /dev/null
USER mambauser
WORKDIR /work
ENTRYPOINT []
CMD ["riboseq-proc", "version"]
