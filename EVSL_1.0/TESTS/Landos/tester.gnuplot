set term post enhanced color landscape 		  # enhanced PostScript, essentially PostScript
  		 		  # with bounding boxes
set out 'tester.eps'              # output file

set title "LanDos Test -- comparison with exact histogram"
set xlabel '{/Symbol l}'
set ylabel 'Spectral Density (DOS)'
set parametric

plot 'OUT/myydos.txt' lt rgb "red" with line title 'Lanczos DOS (G function)~{/Symbol  f}{.4\~}({/Symbol l})', \
     'OUT/myydos2.txt' lt rgb "blue" with line title 'Lanczos DOS (S function) ~{/Symbol  f}{.4\~}({/Symbol l})', \
     'OUT/Exydos.txt' with line title 'Exact Histogram {/Symbol f}({/Symbol l})' 
