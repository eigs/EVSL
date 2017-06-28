set term post enhanced color landscape 		  # enhanced PostScript, essentially PostScript
  		 		  # with bounding boxes
set out 'OUT/testerG.eps'              # output file

set title "LanDos Test -- comparison with exact histogram"
set xlabel '{/Symbol l}'
set ylabel 'Spectral Density (DOS)'
set parametric

plot 'OUT/LanDosG_Approx_DOS.txt' lt rgb "red" with line title 'Lanczos DOS ~{/Symbol  f}{.4\~}({/Symbol l})', \
     'OUT/LanDosG_Exact_DOS.txt' with line title 'Exact Histogram {/Symbol f}({/Symbol l})', \
