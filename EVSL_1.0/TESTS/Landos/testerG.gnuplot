set term post enh 		  # enhanced PostScript, essentially PostScript
  		 		  # with bounding boxes
set out 'testerG.eps'              # output file

set title "LanDos Test -- comparison with exact histogram"
set xlabel '{/Symbol l}'
set ylabel 'Spectral Density (DOS)'
set xrange[0:.05]
set yrange[0:300]
set parametric

plot 'OUT/myydosG.txt' lt rgb "red" with line title 'Lanczos DOS ~{/Symbol  f}{.4\~}({/Symbol l})', \
     'OUT/ExydosG.txt' with line title 'Exact Histogram {/Symbol f}({/Symbol l})', \
