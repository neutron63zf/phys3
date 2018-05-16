set title "Euler and Middle-Point and RK(deg4) behavior by h(step value)"
set ylabel "x[m]"
set xlabel "time[s]"
set key title "value of h"
set key opaque
plot \
    "h4.txt" title "Euler, h = pow(2, -4)" with lines linewidth 1, \
    "h4m.txt" title "Middle Point, h = pow(2, -4)" with lines linewidth 1, \
    "h4r.txt" title "RK(deg4), h = pow(2, -4)" with lines linewidth 1
#    "h7.txt" title "Euler, h = pow(2, -7)" with lines linewidth 1, \
