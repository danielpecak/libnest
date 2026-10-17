IN=plots/
OUT=../_static/
GEN=../_generated/
PYTHON?=python3

# Figures and generated text are rebuilt when their script OR the library code changes.
LIBSRC = $(wildcard ../../libnest/*.py ../../libnest/*/*.py)

# Scripts always plot the code of this checkout, never another installed copy.
export PYTHONPATH := $(abspath ../..)$(if $(PYTHONPATH),:$(PYTHONPATH))

# Figures used in the pages; all of them are also listed in the gallery.
FILES = energySM.png energySMZoom.png rho_vs_kf.png
FILES+= tau_vs_kf.png tau_vs_rho.png
FILES+= vlandau_vs_kf.png vlandau_vs_rho.png
FILES+= pairing_vs_rho.png pairing_vs_kf.png eminigap_vs_kf.png eminigap_vs_rho.png
FILES+= Cs_vs_rho.png pairing_schemes.png functionals.png
# Images that are not figures (not in the gallery).
EXTRA = favicon.png
# Generated rst fragments and files included by the pages.
TEXTS = functionals.inc gallery.inc citation.bib

all: $(addprefix $(OUT), $(FILES) $(EXTRA)) $(addprefix $(GEN), $(TEXTS))

$(OUT) $(GEN):
	mkdir -p $@

$(OUT)%.png: $(IN)plot_%.py $(IN)_common.py $(LIBSRC) | $(OUT)
	$(PYTHON) $< $@

$(GEN)functionals.inc: text/functionals.py $(LIBSRC) | $(GEN)
	$(PYTHON) $< $@

$(GEN)gallery.inc: text/gallery.py $(addprefix $(IN)plot_, $(FILES:.png=.py)) Makefile | $(GEN)
	$(PYTHON) $< $@ $(FILES)

$(GEN)citation.bib: text/citation.py ../../CITATION.cff | $(GEN)
	$(PYTHON) $< ../../CITATION.cff $@

clean:
	rm -f $(addprefix $(OUT), $(FILES) $(EXTRA)) $(addprefix $(GEN), $(TEXTS))

.PHONY: all clean
