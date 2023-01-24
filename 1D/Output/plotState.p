reset
set key font ",12"

nx=512
nt=2000
tskip=50


#set ytics font "Verdana,6" 


do for [j=0:nt/tskip -1] {

  set terminal png size 2560,1440
  
  i = j*tskip
  
  filename="Frames/fields_t=".i.".png"
  set output filename
  
  set title "Time Step = #".(i+1)
  set multiplot layout 2,3 
  
  i = j
  
  set xlabel "x"
  set ylabel "RHO"
  #set yrange [-0.1:1.8]
  plot "fluid.txt" every ::i*nx::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x"
  set ylabel "By"
  #set yrange [-2:2]
  plot "Bfield.txt" every ::i*nx::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x"
  set ylabel "Bz"
  #set yrange [-2:2]
  plot "Bfield.txt" every ::i*nx::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange

  set xlabel "x"
  set ylabel "vx"
  #set yrange [-2:2]
  plot "fluid.txt" every ::i*nx::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x"
  set ylabel "vy"
  #set yrange [-2:2]
  plot "fluid.txt" every ::i*nx::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x"
  set ylabel "vz"
  #set yrange [-2:2]
  plot "fluid.txt" every ::i*nx::nx+(i*nx)-1 using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  unset multiplot
  unset output
    
  
  
}

reset