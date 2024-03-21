OUTPUT = feynman_diagrams.pdf
N_ROWS_IN_MERGED_OUTPUT = 3
N_COLUMNS_IN_MERGED_OUTPUT = 2
COMPILE_LOG = compilation.log

ALLGRAPHS=all_graphs.pdf
GRAPHS_SRC=$(wildcard *.tex)
GRAPHS_PDF=$(patsubst %.tex,%.pdf,$(wildcard *.tex))

# targets: ./0_lbl_box_0.pdf
all: $(OUTPUT)

%.pdf: %.tex
	@echo "\nCompiling "$*.pdf" ...\n"
	latex $*.tex 2>&1 >> $(COMPILE_LOG)
	mpost diagram_$* 2>&1 >> $(COMPILE_LOG)
	latex $*.tex 2>&1 >> $(COMPILE_LOG)
	dvips $*.dvi &> /dev/null
	ps2pdf $*.ps 2>&1 >> $(COMPILE_LOG)
	pdfcrop $*.pdf $*.pdf 2>&1 >> $(COMPILE_LOG)

$(ALLGRAPHS): $(sort $(GRAPHS_PDF))
	pdfunite $^ $@

$(OUTPUT): $(ALLGRAPHS)
	python3 ./combine_pages.py $(ALLGRAPHS) $(OUTPUT) --margins 0.15 0.15 --grid_shape $(N_ROWS_IN_MERGED_OUTPUT) $(N_COLUMNS_IN_MERGED_OUTPUT)
	@echo "\n> Feynman diagrams all written to file '"$(OUTPUT)"'\n"

clean:
	rm -f *.t1 *.log *.dvi *.mp *.log *.1 *.ps *.aux

clean_pdf:
	rm *.pdf