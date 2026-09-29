#!/bin/bash
# Mirrors the CI lint and test jobs. Run: ./scripts/check-before-push.sh
set -e

echo "1. ruff check"
python -m ruff check .
echo "2. ruff format --check"
python -m ruff format --check .
echo "3. pytest"
python -m pytest --tb=short -q

echo "All checks passed."
