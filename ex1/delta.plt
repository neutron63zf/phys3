set title "Behavior of Error of Numerical Differencial of sin(x) by h(step value)"
set xlabel "h"
set ylabel "error"
set logscale y
set logscale x
set key opaque
set key title "procedure name"
plot \
    "point2.txt" title "2 point differential" with lines linewidth 1, \
    "point3.txt" title "3 point differential" with lines linewidth 1, \
    "point5.txt" title "5 point differential" with lines linewidth 1
set terminal postscript eps
set output "delta.eps"
replot
set terminal x11