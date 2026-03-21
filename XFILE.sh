#!/bin/bash
PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}" \
PYTHONWARNINGS='ignore:invalid escape sequence:SyntaxWarning' \
snakemake --cores all --rerun-incomplete --rerun-triggers mtime -p --latency-wait 20
