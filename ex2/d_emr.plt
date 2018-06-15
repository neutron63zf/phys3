set terminal x11
set title "Simulation behavior by procedure type"
set ylabel "error"
set xlabel "time[s]"
#A=10;kappa=0.2;k=2;m=1;
#gamma=kappa/(2*m);w=sqrt(k/m-gamma**2.0);
set key title "type"
set key opaque
plot \
    "dh4.txt" title "Euler" with lines linewidth 1, \
    "dh4m.txt" title "Middle Point" with lines linewidth 1, \
    "dh4r.txt" title "Runge Kutta" with lines linewidth 1
set terminal postscript eps
set output "d_emr.eps"
replot