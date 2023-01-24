reset
set key font ",12"

nx=100
nt=80
tskip=1


#set ytics font "Verdana,6" 



do for [j=0:nt/tskip -1] {

  set terminal png size 2560,1440
  
  i = j*tskip
  
  filename="Frames/fields_diagcut_t=".i.".png"
  set output filename
  
  set title "Time Step = #".(i+1)
  set multiplot layout 2,3 
  
  file1 = "fluid_diagcut.txt"
  file2 = "Bfield_diagcut.txt"
  
  i = j
  
  set xlabel "x=y"
  set ylabel "RHO"
  #set yrange [-0.1:1.8]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x=y"
  set ylabel "Bpar"
  #set yrange [-2:2]
  plot file2 every ::i*nx::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x=y"
  set ylabel "Bz"
  #set yrange [-2:2]
  plot file2 every ::i*nx::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange

  set xlabel "x=y"
  set ylabel "vperp"
  #set yrange [-2:2]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x=y"
  set ylabel "vpar"
  #set yrange [-2:2]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x=y"
  set ylabel "vz"
  #set yrange [-2:2]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  unset multiplot
  unset output
    
  
  
}


EXIT


do for [j=0:nt/tskip -1] {

  set terminal png size 2560,1440
  
  i = j*tskip
  
  filename="Frames/fields_xcut_t=".i.".png"
  set output filename
  
  set title "Time Step = #".(i+1)
  set multiplot layout 2,3 
  
  file1 = "fluid_xcut.txt"
  file2 = "Bfield_xcut.txt"
  
  i = j
  
  set xlabel "x"
  set ylabel "RHO"
  #set yrange [-0.1:1.8]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x"
  set ylabel "By"
  #set yrange [-2:2]
  plot file2 every ::i*nx::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x"
  set ylabel "Bz"
  #set yrange [-2:2]
  plot file2 every ::i*nx::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange

  set xlabel "x"
  set ylabel "vx"
  #set yrange [-2:2]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x"
  set ylabel "vy"
  #set yrange [-2:2]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x"
  set ylabel "vz"
  #set yrange [-2:2]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  unset multiplot
  unset output
    
  
  
}


EXIT

do for [j=0:nt/tskip -1] {

  set terminal png size 2560,1440
  
  i = j*tskip
  
  filename="Frames/fields_ycut_t=".i.".png"
  set output filename
  
  set title "Time Step = #".(i+1)
  set multiplot layout 2,3 
  
  file1 = "fluid_ycut.txt"
  file2 = "Bfield_ycut.txt"
  
  i = j
  
  set xlabel "y"
  set ylabel "RHO"
  #set yrange [-0.1:1.8]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "y"
  set ylabel "Bx"
  #set yrange [-2:2]
  plot file2 every ::i*nx::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "y"
  set ylabel "Bz"
  #set yrange [-2:2]
  plot file2 every ::i*nx::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange

  set xlabel "y"
  set ylabel "vx"
  #set yrange [-2:2]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "y"
  set ylabel "vy"
  #set yrange [-2:2]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "y"
  set ylabel "vz"
  #set yrange [-2:2]
  plot file1 every ::i*nx::nx+(i*nx)-1 using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
  unset yrange
  
  unset multiplot
  unset output
    
  
  
}



reset