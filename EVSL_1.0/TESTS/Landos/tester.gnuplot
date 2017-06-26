set term post enhanced color landscape 		  # enhanced PostScript, essentially PostScript
  		 		  # with bounding boxes
set out 'OUT/tester.eps'              # output file

set title "LanDos Test -- comparison with exact histogram"
set xlabel '{/Symbol l}'
set ylabel 'Spectral Density (DOS)'
set parametric

plot 'OUT/LanDos_Approx_DOS.txt' lt rgb "blue" with line title 'Lanczos DOS ~{/Symbol  f}{.4\~}({/Symbol l})', \
     'OUT/LanDos_Exact_DOS.txt' with line title 'Exact Histogram {/Symbol f}({/Symbol l})' 
