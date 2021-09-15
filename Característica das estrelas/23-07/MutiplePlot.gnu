set terminal postscript eps enhanced color "Helvetica" 24 dashed dashlength 0.5
set datafile fortran
#multiplot com cada plot por interaccao

set autoscale
set output 'MR-new-07.eps'
set size 1.,2.
set origin 0.,0

set multiplot


unset logscale
set size 1,1.5
set origin 0.,0.2
set xrange [9:15]
set yrange [0:1]
set xtics 9, 1, 18
#set ytics 0.4, 0.4, 2.
set mytics 5
set mxtics 10
set xlabel 'R    [km]'
set ylabel 'rh_c [fm^-4]'
unset key
set key below
set label 'eos-bbp-nl3-1star' at 14.1,2.8



p 'tov0.out' u 1:4 w l lt 1 lw 2 lc 0 dt 1 t'1,0sM',\
'tov2.out' u 1:4 w l lt 1 lw 2 lc 1 dt 1 t'1,2sM',\
'tov4.out' u 1:4 w l lt 1 lw 2 lc 2 dt 1 t'1,4sM',\
'tov6.out' u 1:4 w l lt 1 lw 2 lc 3 dt 1 t'1,6sM',\
'tov8.out' u 1:4 w l lt 1 lw 2 lc 4 dt 1 t'1,8sM',\
'tovmax.out' u 1:4 w l lt 1 lw 2 lc 6 dt 1 t'Max. mass',\




#'FamMagCMFFreq0Mu2e31.d' u 8:11 w l lt 1 lw 4 lc 5 dt 4 t'',\


unset multiplot
