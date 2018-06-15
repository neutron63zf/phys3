set terminal x11
set title "Euler Simulation behavior by h(step value)"
set ylabel "x[m]"
set xlabel "time[s]"
A=10;kappa=0.2;k=2;m=1;
gamma=kappa/(2*m);w=sqrt(k/m-gamma**2.0);
set key title "value of h"
set key opaque
plot \
    "h3.txt" title "pow(2, -3)" with lines linewidth 1, \
    "h4.txt" title "pow(2, -4)" with lines linewidth 1, \
    "h5.txt" title "pow(2, -5)" with lines linewidth 1, \
    "h6.txt" title "pow(2, -6)" with lines linewidth 1, \
    "h7.txt" title "pow(2, -7)" with lines linewidth 1, \
    A*exp(-gamma*x)*cos(w*x) title "zero (real value)"
set terminal postscript eps
set output "euler.eps"
replot