set terminal x11
set title "Euler Simulation behavior by h(step value)"
set ylabel "error"
set xlabel "time[s]"
#A=10;kappa=0.2;k=2;m=1;
#gamma=kappa/(2*m);w=sqrt(k/m-gamma**2.0);
set key title "value of h"
set key opaque
plot \
    "dh4.txt" title "pow(2, -4)" with lines linewidth 1, \
    "dh5.txt" title "pow(2, -5)" with lines linewidth 1, \
    "dh6.txt" title "pow(2, -6)" with lines linewidth 1, \
    "dh7.txt" title "pow(2, -7)" with lines linewidth 1
set terminal postscript eps
set output "d_euler.eps"
replot