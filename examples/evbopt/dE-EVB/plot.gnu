set terminal svg lw 1.5 size 500,400
set output "evbopt_de.svg"
set xlabel "Frame"
set ylabel "Energie (Hartree)"
set xrange [0:70]
set yrange [-272.23:-272.11]
plot "ref.dat" title "Reference" w l, "energies.qmdff" title "EVB-QMDFF" w l, "single_qmdff.dat" u 0:1 title "QMDFF1" w l, "single_qmdff.dat" u 0:2 title "QMDFF2" w l
pause -1
