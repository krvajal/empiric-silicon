# set terminal pdfcairo
# set output "fig1.pdf"
unset key
unset xtics
set xtics ('Γ' 0, "X" 0.92102998110955914, "W" 1.3815449716643387,\
            "L" 2.0328115169830259,\
            "Γ" 2.8304468782710055,\
            "K" 3.8073466962490361)
set ylabel "Energia [Ry]"
set yrange [0.367:0.4116    ]
set xrange [0:3.8073466962490361]
set arrow from 0.92102998110955914,-0.05 to  0.92102998110955914, 0.7 nohead dt 2
set arrow from 1.3815449716643387,-0.05 to  1.3815449716643387, 0.7 nohead dt 2
set arrow from 2.0328115169830259,-0.05 to  2.0328115169830259, 0.7 nohead dt 2
set arrow from 2.8304468782710055,-0.05 to  2.8304468782710055, 0.7 nohead dt 2
plot '../bands.txt' u 1:2 w l , '../bands.txt' u 1:4 w l, '../bands.txt' u 1:5 w l, '../bands.txt' u 1:3 w l, '../bands.txt' u 1:6 w l, '../bands.txt' u 1:7 w l , '../bands.txt' u 1:8 w l , '../bands.txt' u 1:9 w l, '../bands.txt' u 1:10 w l

set output
