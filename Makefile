MAIN ?= main
SRC := $(wildcard *.tex)
BIB := biblio.bib ../optimal-gartner/cryptobib/abbrev3.bib ../optimal-gartner/cryptobib/crypto.bib
PDFLATEX := $(shell command -v pdflatex 2>/dev/null)
# Prefer TeX Live over an older MiKTeX executable in /usr/local/bin.
TEXLIVE_BIB := $(firstword $(wildcard /usr/local/texlive/*/bin/*/bibtex8))
BIBTEX8 ?= $(if $(TEXLIVE_BIB),$(TEXLIVE_BIB),bibtex8)

SAGE ?= sage
ESTIMATOR ?= ../../lattice-estimator

.PHONY: all verify security
all: $(MAIN).pdf

ifneq ($(PDFLATEX),)
$(MAIN).pdf: $(SRC) $(BIB) iacrcc.cls after-hyperref.sty
	$(PDFLATEX) -recorder -halt-on-error -interaction=nonstopmode $(MAIN).tex
	$(BIBTEX8) -7 --huge $(MAIN)
	$(PDFLATEX) -recorder -halt-on-error -interaction=nonstopmode $(MAIN).tex
	$(PDFLATEX) -recorder -halt-on-error -interaction=nonstopmode $(MAIN).tex
else
$(MAIN).pdf: $(SRC) $(BIB) iacrcc.cls after-hyperref.sty
	tectonic --pass tex --keep-intermediates $(MAIN).tex
	$(BIBTEX8) -7 --huge $(MAIN)
	tectonic --pass tex --keep-intermediates $(MAIN).tex
	tectonic --pass tex --keep-intermediates $(MAIN).tex
	xdvipdfmx -q $(MAIN).xdv
endif

verify:
	python3 verify.py
	python3 compare.py
	python3 compare_four.py
	python3 check_protocol.py
	python3 check_multimodal.py
	python3 check_four_compression.py

security:
	$(SAGE) -python security_estimates.py --estimator $(ESTIMATOR)
