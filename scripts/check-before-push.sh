#!/bin/bash
# Local CI checks before pushing to master/dev
# Run: ./scripts/check-before-push.sh

set -e

echo "Running pre-push checks..."
echo

# Lint
echo "1. Linting with Ruff..."
python -m ruff check src tests
echo "✓ Lint passed"
echo

# Tests
echo "2. Running pytest (full suite)..."
python -m pytest --tb=short -q
echo "✓ Tests passed"
echo

echo "All checks passed! Safe to push."
