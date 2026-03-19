PYTHON ?= python3
VENV ?= .venv
ACTIVATE = . $(VENV)/bin/activate

.PHONY: venv install test run check-deps dist clean

venv:
	$(PYTHON) -m venv $(VENV)
	$(ACTIVATE) && python -m pip install --upgrade pip

install:
	$(ACTIVATE) && pip install -e .

test:
	$(ACTIVATE) && pytest -v

dist:
	$(ACTIVATE) && python -m build

run:
	$(ACTIVATE) && sh XFILE.sh

check-deps:
	$(ACTIVATE) && rnaconsnake-run --check-deps

clean:
	rm -rf $(VENV) build dist *.egg-info
