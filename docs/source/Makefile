IN=plots/
OUT=../_static/



FILES = energySM.png energySMZoom.png rho_vs_kf.png
FILES+= tau_vs_kf.png tau_vs_rho.png
FILES+= vlandau_vs_kf.png vlandau_vs_rho.png
FILES+= pairing_vs_rho.png pairing_vs_kf.png eminigap_vs_kf.png eminigap_vs_rho.png
FILES+= Cs_vs_rho.png
all: $(addprefix $(OUT), $(FILES))
	mkdir -p $(OUT)


$(OUT)%.png: $(IN)plot_%.py
	$(PYTHON) $^ $@

clean:
	rm $(addprefix $(OUT), $(FILES))
