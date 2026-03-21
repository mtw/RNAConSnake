PYTHON ?= python3
VENV ?= .venv
ACTIVATE = . $(VENV)/bin/activate

.PHONY: venv install test run check-deps dist clean clean-data

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
	@if [ -z "$(INPUT)" ]; then echo "Usage: make run INPUT=/path/to/input.stk"; exit 2; fi
	$(ACTIVATE) && sh XFILE.sh "$(INPUT)" $(if $(OUTPUT),--output-dir "$(OUTPUT)")

check-deps:
	$(ACTIVATE) && rnaconsnake-run --check-deps

clean-data:
	$(ACTIVATE) && snakemake --cores 1 clean

clean:
	rm -rf $(VENV) build dist *.egg-info
