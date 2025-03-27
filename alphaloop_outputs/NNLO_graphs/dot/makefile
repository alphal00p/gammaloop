OUTPUT = feynman_diagrams.pdf
N_ROWS_IN_MERGED_OUTPUT = 3
N_COLUMNS_IN_MERGED_OUTPUT = 2
COMPILE_LOG = compilation.log

ALLGRAPHS=all_graphs.pdf
GRAPHS_SRC=$(wildcard *.dot)
# GRAPHS_PDF=$(sort $(patsubst %.tex,%.pdf,$(wildcard *.tex)))
GRAPHS_PDF= ./0_0_GL0.pdf ./1_0_GL1.pdf ./2_0_GL4.pdf ./3_0_GL5.pdf ./4_0_GL8.pdf ./5_0_GL9.pdf ./6_0_GL12.pdf ./7_0_GL13.pdf ./8_0_GL17.pdf ./9_0_GL18.pdf ./10_0_GL24.pdf ./11_0_GL26.pdf ./12_0_GL28.pdf ./13_0_GL30.pdf ./14_0_GL32.pdf ./15_0_GL35.pdf ./16_0_GL36.pdf ./17_0_GL38.pdf ./18_0_GL40.pdf ./19_0_GL42.pdf ./20_0_GL44.pdf ./21_0_GL47.pdf ./22_0_GL60.pdf ./23_0_GL62.pdf ./24_0_GL64.pdf ./25_0_GL66.pdf ./26_0_GL72.pdf ./27_0_GL74.pdf ./28_0_GL76.pdf ./29_0_GL78.pdf ./30_0_GL80.pdf ./31_0_GL82.pdf ./32_0_GL84.pdf ./33_0_GL86.pdf ./34_0_GL88.pdf ./35_0_GL94.pdf ./36_0_GL96.pdf ./37_0_GL102.pdf ./38_0_GL104.pdf ./39_0_GL106.pdf ./40_0_GL112.pdf ./41_0_GL114.pdf ./42_0_GL120.pdf ./43_0_GL122.pdf ./44_0_GL124.pdf ./45_0_GL126.pdf ./46_0_GL128.pdf ./47_0_GL131.pdf ./48_0_GL132.pdf ./49_0_GL138.pdf ./50_0_GL140.pdf ./51_0_GL146.pdf ./52_0_GL148.pdf ./53_0_GL154.pdf ./54_0_GL160.pdf ./55_0_GL162.pdf ./56_0_GL164.pdf ./57_0_GL166.pdf ./58_0_GL168.pdf ./59_0_GL170.pdf ./60_0_GL172.pdf ./61_0_GL174.pdf ./62_0_GL176.pdf ./63_0_GL178.pdf ./64_0_GL180.pdf ./65_0_GL182.pdf ./66_0_GL187.pdf ./67_0_GL189.pdf ./68_0_GL192.pdf ./69_0_GL194.pdf ./70_0_GL196.pdf ./71_0_GL198.pdf ./72_0_GL205.pdf ./73_0_GL207.pdf ./74_0_GL212.pdf ./75_0_GL214.pdf ./76_0_GL216.pdf ./77_0_GL218.pdf ./78_0_GL220.pdf ./79_0_GL222.pdf ./80_0_GL227.pdf ./81_0_GL228.pdf ./82_0_GL234.pdf ./83_0_GL236.pdf ./84_0_GL238.pdf ./85_0_GL240.pdf ./86_0_GL242.pdf ./87_0_GL245.pdf ./88_0_GL246.pdf ./89_0_GL248.pdf ./90_0_GL254.pdf ./91_0_GL260.pdf ./92_0_GL263.pdf ./93_0_GL272.pdf ./94_0_GL274.pdf ./95_0_GL276.pdf ./96_0_GL278.pdf ./97_0_GL284.pdf ./98_0_GL286.pdf ./99_0_GL292.pdf ./100_0_GL294.pdf ./101_0_GL300.pdf ./102_0_GL302.pdf ./103_0_GL304.pdf ./104_0_GL307.pdf ./105_0_GL308.pdf ./106_0_GL310.pdf ./107_0_GL312.pdf ./108_0_GL314.pdf ./109_0_GL316.pdf ./110_0_GL318.pdf ./111_0_GL320.pdf ./112_0_GL322.pdf ./113_0_GL324.pdf ./114_0_GL333.pdf ./115_0_GL338.pdf ./116_0_GL342.pdf ./117_0_GL346.pdf ./118_0_GL358.pdf ./119_0_GL362.pdf ./120_0_GL368.pdf ./121_0_GL386.pdf ./122_0_GL388.pdf ./123_0_GL390.pdf ./124_0_GL392.pdf ./125_0_GL394.pdf ./126_0_GL396.pdf ./127_0_GL398.pdf ./128_0_GL400.pdf ./129_0_GL403.pdf ./130_0_GL404.pdf ./131_0_GL406.pdf ./132_0_GL408.pdf ./133_0_GL410.pdf ./134_0_GL412.pdf ./135_0_GL418.pdf ./136_0_GL421.pdf ./137_0_GL426.pdf ./138_0_GL428.pdf ./139_0_GL430.pdf ./140_0_GL432.pdf ./141_0_GL434.pdf ./142_0_GL440.pdf ./143_0_GL442.pdf ./144_0_GL448.pdf ./145_0_GL451.pdf ./146_0_GL453.pdf ./147_0_GL454.pdf ./148_0_GL464.pdf ./149_0_GL466.pdf ./150_0_GL473.pdf ./151_0_GL474.pdf ./152_0_GL476.pdf ./153_0_GL478.pdf ./154_0_GL491.pdf ./155_0_GL496.pdf ./156_0_GL508.pdf ./157_0_GL552.pdf ./158_0_GL556.pdf

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
