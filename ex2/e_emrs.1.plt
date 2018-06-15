set title "Simulation behavior by procedure type"
set ylabel "Energy"
set xlabel "time[s]"
#A=10;kappa=0.2;k=2;m=1;
#gamma=kappa/(2*m);w=sqrt(k/m-gamma**2.0);
set yrange[99.9965:100.0005]
set key title "type"
set key opaque
plot \
    "eh7r.txt" title "Runge Kutta" with lines linewidth 1, \
    "eh7s.txt" title "Simplectic" with lines linewidth 1
set terminal postscript eps
set output "e_emrs1.eps"
replot
set terminal x11