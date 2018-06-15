set title "Newton behavior by start point"
set ylabel ""
set xlabel ""
#set logscale y
#set logscale x
set key title "start point"
set key opaque
plot \
    "first_1.txt" title "1" with lines linewidth 1, \
    "first_1p5.txt" title "1.5" with lines linewidth 1, \
    "first_1p25.txt" title "1.25" with lines linewidth 1, \
    "first_m1.txt" title "-1" with lines linewidth 1, \
    "first_m100.txt" title "-100" with lines linewidth 1
set terminal postscript eps
set output "newton.eps"
replot