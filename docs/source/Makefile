IN=_static/
OUT=plots/

FILES = energySM.png energySMZoom.png rho_vs_kf.png
FILES+= tau_vs_kf.png tau_vs_rho.png
FILES+= vlandau_vs_kf.png vlandau_vs_rho.png
FILES+= pairing_vs_rho.png pairing_vs_kf.png eminigap_vs_kf.png eminigap_vs_rho.png
all: $(addprefix $(IN), $(FILES))


$(IN)%.png: $(OUT)plot_%.py
	$(PYTHON) $^ $@
