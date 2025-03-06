OUTPUT = feynman_diagrams.pdf
N_ROWS_IN_MERGED_OUTPUT = 3
N_COLUMNS_IN_MERGED_OUTPUT = 2
COMPILE_LOG = compilation.log

ALLGRAPHS=all_graphs.pdf
GRAPHS_SRC=$(wildcard *.dot)
# GRAPHS_PDF=$(sort $(patsubst %.tex,%.pdf,$(wildcard *.tex)))
GRAPHS_PDF= ./0_0_GL0.pdf ./1_0_GL4.pdf ./2_0_GL8.pdf ./3_0_GL12.pdf ./4_0_GL16.pdf ./5_0_GL20.pdf ./6_0_GL24.pdf ./7_0_GL28.pdf ./8_0_GL32.pdf ./9_0_GL34.pdf ./10_0_GL40.pdf ./11_0_GL42.pdf ./12_0_GL44.pdf ./13_0_GL50.pdf ./14_0_GL56.pdf ./15_0_GL60.pdf ./16_0_GL69.pdf ./17_0_GL72.pdf ./18_0_GL76.pdf ./19_0_GL82.pdf ./20_0_GL88.pdf ./21_0_GL92.pdf ./22_0_GL120.pdf ./23_0_GL124.pdf ./24_0_GL128.pdf ./25_0_GL130.pdf ./26_0_GL136.pdf ./27_0_GL139.pdf ./28_0_GL140.pdf ./29_0_GL146.pdf ./30_0_GL156.pdf ./31_0_GL158.pdf ./32_0_GL166.pdf ./33_0_GL174.pdf ./34_0_GL178.pdf ./35_0_GL209.pdf ./36_0_GL210.pdf ./37_0_GL218.pdf ./38_0_GL224.pdf ./39_0_GL235.pdf ./40_0_GL236.pdf ./41_0_GL244.pdf ./42_0_GL264.pdf ./43_0_GL268.pdf ./44_0_GL272.pdf ./45_0_GL279.pdf ./46_0_GL280.pdf ./47_0_GL284.pdf ./48_0_GL290.pdf ./49_0_GL302.pdf ./50_0_GL306.pdf ./51_0_GL313.pdf ./52_0_GL316.pdf ./53_0_GL322.pdf ./54_0_GL336.pdf ./55_0_GL340.pdf ./56_0_GL344.pdf ./57_0_GL348.pdf ./58_0_GL352.pdf ./59_0_GL356.pdf ./60_0_GL360.pdf ./61_0_GL364.pdf ./62_0_GL368.pdf ./63_0_GL370.pdf ./64_0_GL372.pdf ./65_0_GL374.pdf ./66_0_GL376.pdf ./67_0_GL383.pdf ./68_0_GL384.pdf ./69_0_GL386.pdf ./70_0_GL388.pdf ./71_0_GL390.pdf ./72_0_GL404.pdf ./73_0_GL410.pdf ./74_0_GL428.pdf ./75_0_GL434.pdf ./76_0_GL440.pdf ./77_0_GL444.pdf ./78_0_GL448.pdf ./79_0_GL453.pdf ./80_0_GL458.pdf ./81_0_GL470.pdf ./82_0_GL488.pdf ./83_0_GL491.pdf ./84_0_GL493.pdf ./85_0_GL499.pdf ./86_0_GL506.pdf ./87_0_GL508.pdf ./88_0_GL511.pdf ./89_0_GL514.pdf ./90_0_GL533.pdf ./91_0_GL551.pdf ./92_0_GL553.pdf ./93_0_GL582.pdf ./94_0_GL584.pdf ./95_0_GL588.pdf ./96_0_GL594.pdf ./97_0_GL606.pdf ./98_0_GL610.pdf ./99_0_GL614.pdf ./100_0_GL622.pdf ./101_0_GL627.pdf ./102_0_GL635.pdf ./103_0_GL636.pdf ./104_0_GL640.pdf ./105_0_GL644.pdf ./106_0_GL648.pdf ./107_0_GL653.pdf ./108_0_GL656.pdf ./109_0_GL662.pdf ./110_0_GL664.pdf ./111_0_GL668.pdf ./112_0_GL672.pdf ./113_0_GL679.pdf ./114_0_GL697.pdf ./115_0_GL712.pdf ./116_0_GL723.pdf ./117_0_GL732.pdf ./118_0_GL767.pdf ./119_0_GL768.pdf ./120_0_GL776.pdf ./121_0_GL828.pdf ./122_0_GL830.pdf ./123_0_GL832.pdf ./124_0_GL834.pdf ./125_0_GL837.pdf ./126_0_GL844.pdf ./127_0_GL852.pdf ./128_0_GL854.pdf ./129_0_GL856.pdf ./130_0_GL858.pdf ./131_0_GL860.pdf ./132_0_GL862.pdf ./133_0_GL864.pdf ./134_0_GL866.pdf ./135_0_GL872.pdf ./136_0_GL874.pdf ./137_0_GL880.pdf ./138_0_GL883.pdf ./139_0_GL888.pdf ./140_0_GL890.pdf ./141_0_GL897.pdf ./142_0_GL904.pdf ./143_0_GL906.pdf ./144_0_GL912.pdf ./145_0_GL914.pdf ./146_0_GL916.pdf ./147_0_GL923.pdf ./148_0_GL944.pdf ./149_0_GL946.pdf ./150_0_GL952.pdf ./151_0_GL955.pdf ./152_0_GL958.pdf ./153_0_GL964.pdf ./154_0_GL988.pdf ./155_0_GL1012.pdf ./156_0_GL1054.pdf ./157_0_GL1104.pdf ./158_0_GL1108.pdf

all: $(OUTPUT)

%.pdf: %.dot
	@echo "\nCompiling "$*.pdf" ...\n"
	dot -Tpdf $*.dot -o $*.pdf 2>&1 >> $(COMPILE_LOG)
	sleep 0.1
	pdfcrop $*.pdf $*.pdf 2>&1 >> $(COMPILE_LOG)

$(ALLGRAPHS): $(GRAPHS_PDF)
	pdfunite $^ $@

$(OUTPUT): $(ALLGRAPHS)
	python3 ./combine_pages.py $(ALLGRAPHS) $(OUTPUT) --margins 0.15 0.15 --grid_shape $(N_ROWS_IN_MERGED_OUTPUT) $(N_COLUMNS_IN_MERGED_OUTPUT)
	@echo "\n> Feynman diagrams all written to file '"$(OUTPUT)"'\n"

clean:
	rm -f *.t1 *.log *.dvi *.mp *.log *.1 *.ps *.aux

clean_pdf:
	rm *.pdf
