set term post enh 		  # enhanced PostScript, essentially PostScript
  		 		  # with bounding boxes
set out 'tester.eps'              # output file

set title "LanDos Test -- comparison with exact histogram"
set xlabel '{/Symbol l}'
set ylabel 'Spectral Density (DOS)'
set xrange[-2.5:12]
set yrange[0:0.20]
set parametric

plot 'OUT/myydos.txt' lt rgb "red" with line title 'Lanczos DOS ~{/Symbol  f}{.4\~}({/Symbol l})', \
     'OUT/Exydos.txt' with line title 'Exact Histogram {/Symbol f}({/Symbol l})', \
