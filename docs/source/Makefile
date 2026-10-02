IN=plots/
OUT=../_static/
PYTHON?=python3

# Figures are regenerated when their script OR the library code changes.
LIBSRC = $(wildcard ../../libnest/*.py)

# Scripts always plot the code of this checkout, never another installed copy.
export PYTHONPATH := $(abspath ../..)$(if $(PYTHONPATH),:$(PYTHONPATH))

FILES = energySM.png energySMZoom.png rho_vs_kf.png
FILES+= tau_vs_kf.png tau_vs_rho.png
FILES+= vlandau_vs_kf.png vlandau_vs_rho.png
FILES+= pairing_vs_rho.png pairing_vs_kf.png eminigap_vs_kf.png eminigap_vs_rho.png
FILES+= Cs_vs_rho.png

all: $(addprefix $(OUT), $(FILES))

$(OUT):
	mkdir -p $(OUT)

$(OUT)%.png: $(IN)plot_%.py $(IN)_common.py $(LIBSRC) | $(OUT)
	$(PYTHON) $< $@

clean:
	rm -f $(addprefix $(OUT), $(FILES))

.PHONY: all clean
