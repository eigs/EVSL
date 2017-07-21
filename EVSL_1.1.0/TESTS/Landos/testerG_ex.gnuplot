set term post enhanced color landscape 		  # enhanced PostScript, essentially PostScript
  		 		  # with bounding boxes
set out filename.'.eps'              # output file

set title "LanDos Test -- comparison with exact histogram"
set xlabel '{/Symbol l}'
set ylabel 'Spectral Density (DOS)'
set parametric

plot filename lt rgb "red" with line title 'Lanczos DOS ~{/Symbol  f}{.4\~}({/Symbol l})', \
     exact_dos with line title 'Exact Histogram {/Symbol f}({/Symbol l})'
