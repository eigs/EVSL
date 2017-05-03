set term post enh 		  # enhanced PostScript, essentially PostScript
 		 		  # with bounding boxes
set out 'tester.eps'              # output file

set title "test LanDos"
set xlabel '{/Symbol l}'
set ylabel 'Spectral Density (DOS)'
set xrange[-2.5:12]
set yrange[0:0.20]
set parametric

plot 'OUT/myydos.txt' lt rgb "red" with line title 'Lan DOS {/Symbol f}({/Symbol l})', \
     'OUT/Exydos.txt' with line title 'Ex Dos {/Symbol f}({/Symbol l})', \
