IN=_static/
OUT=plots/

FILES= energySM.png energySMZoom.png

all: $(addprefix $(IN), $(FILES))

$(IN)%.png: $(OUT)plot_%.py
	./$^ $@
