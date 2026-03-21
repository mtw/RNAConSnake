#!/bin/bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: sh XFILE.sh /path/to/input.stk [additional rnaconsnake-run args...]"
  exit 2
fi

input_alignment="$1"
shift

PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}" \
PYTHONWARNINGS='ignore:invalid escape sequence:SyntaxWarning' \
rnaconsnake-run --input-alignment "$input_alignment" --conservative --cores all "$@"
