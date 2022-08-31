IN=_static/
OUT=plots/

FILES= energySM.png energySMZoom.png rho_vs_kf.png eminigap_vs_kf.png eminigap_vs_rho.png tau_vs_kf.png tau_vs_rho.png vlandau_vs_kf.png vlandau_vs_rho.png pairing_vs_rho.png pairing_vs_kf.png

all: $(addprefix $(IN), $(FILES))

$(IN)%.png: $(OUT)plot_%.py
	./$^ $@
