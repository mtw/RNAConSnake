#!/bin/bash
PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}" \
PYTHONWARNINGS='ignore:invalid escape sequence:SyntaxWarning' \
snakemake --cores all --rerun-incomplete -p --latency-wait 20
