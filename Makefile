PYTHON ?= python3
VENV := .venv
PY := $(VENV)/bin/python

.PHONY: setup pipeline dashboard test lint
setup:
	$(PYTHON) -m venv $(VENV)
	$(PY) -m pip install -r requirements.txt

pipeline:
	$(PY) load_data.py
	$(PY) analysis.py

dashboard:
	$(PY) -m streamlit run app.py --server.address 0.0.0.0 --server.port 8501 --server.headless true

test:
	$(PY) -m pytest -q

lint:
	$(PY) -m ruff check .
	$(PY) -m ruff format --check .
