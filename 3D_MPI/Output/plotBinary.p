reset

set key font ",12"

nx=32
ny=16
nz=16
nt=120
tskip=20


E_max = 0.1

set terminal png size 1920,1080
#set size ratio -1

#################
# x cut
#################



do for [j=0:nt/tskip] {

  set terminal png size 2560,1440
  
  i = j*tskip
  
  filename="Frames/fields_x_cut_t=".i.".png"
  set output filename
  
  set title "Time Step = #".(i+1)
  
  set multiplot layout 2,3 
    
  file1="Snapshots/xcut_rank=011_t=".i.".dat"
  file2="Snapshots/xcut_rank=111_t=".i.".dat"
  file3="Snapshots/xcut_rank=211_t=".i.".dat"
  
  i = j
  
  set xlabel "x"
  set ylabel "RHO"
  #set yrange [-0.1:1.8]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:2 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:2 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:2 with linespoint pointtype 4 lc rgb "green" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x"
  set ylabel "By"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:7 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:7 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:7 with linespoint pointtype 4 lc rgb "green" notitle
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x"
  set ylabel "Bz"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:8 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:8 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:8 with linespoint pointtype 4 lc rgb "green" notitle
  unset xlabel
  unset ylabel
  unset yrange

  set xlabel "x"
  set ylabel "vx"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:3 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:3 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:3 with linespoint pointtype 4 lc rgb "green" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "x"
  set ylabel "vy"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:4 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:4 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:4 with linespoint pointtype 4 lc rgb "green" notitle
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "x"
  set ylabel "vz"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:5 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:5 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:5 with linespoint pointtype 4 lc rgb "green" notitle
  unset xlabel
  unset xlabel
  unset ylabel
  unset yrange
  

  unset title
  unset multiplot
  unset output
  
  
}


exit


#################
# z cut
#################



do for [j=0:nt/tskip] {

  set terminal png size 2560,1440
  
  i = j*tskip
  
  filename="Frames/fields_z_cut_t=".i.".png"
  set output filename
  
  set title "Time Step = #".(i+1)
  
  set multiplot layout 2,3 
    
  file1="Snapshots/zcut_rank=000_t=".i.".dat"
  file2="Snapshots/zcut_rank=001_t=".i.".dat"
  file3="Snapshots/zcut_rank=002_t=".i.".dat"
  
  i = j
  
  set xlabel "z"
  set ylabel "RHO"
  #set yrange [-0.1:1.8]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:2 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:2 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:2 with linespoint pointtype 4 lc rgb "green" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "z"
  set ylabel "By"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:7 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:7 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:7 with linespoint pointtype 4 lc rgb "green" notitle
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "z"
  set ylabel "Bx"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:6 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:6 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:6 with linespoint pointtype 4 lc rgb "green" notitle
  unset xlabel
  unset ylabel
  unset yrange

  set xlabel "z"
  set ylabel "vx"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:3 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:3 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:3 with linespoint pointtype 4 lc rgb "green" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "z"
  set ylabel "vy"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:4 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:4 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:4 with linespoint pointtype 4 lc rgb "green" notitle
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "z"
  set ylabel "vz"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:5 with linespoint pointtype 8 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:5 with linespoint pointtype 4 lc rgb "red" notitle,\
  file3 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:5 with linespoint pointtype 4 lc rgb "green" notitle
  unset xlabel
  unset xlabel
  unset ylabel
  unset yrange
  

  unset title
  unset multiplot
  unset output
  
  
}

exit






#################
# y cut
#################



do for [j=0:nt/tskip -1] {

  set terminal png size 2560,1440
  
  i = j*tskip
  
  filename="Frames/fields_y_cut_t=".i.".png"
  set output filename
  
  set title "Time Step = #".(i+1)
  set multiplot layout 2,3 
  
  file1="Snapshots/ycut_rank=000_t=".i.".dat"
  file2="Snapshots/ycut_rank=010_t=".i.".dat"
  
  i = j
  
  set xlabel "y"
  set ylabel "RHO"
  #set yrange [-0.1:1.8]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:2 with linespoint pointtype 7 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "y"
  set ylabel "Bx"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:6 with linespoint pointtype 7 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:6 with linespoint pointtype 7 lc rgb "red" notitle
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "y"
  set ylabel "Bz"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:8 with linespoint pointtype 7 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:8 with linespoint pointtype 7 lc rgb "red" notitle
  unset xlabel
  unset ylabel
  unset yrange

  set xlabel "y"
  set ylabel "vx"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:3 with linespoint pointtype 7 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "y"
  set ylabel "vy"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:4 with linespoint pointtype 7 lc rgb "red" notitle
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "y"
  set ylabel "vz"
  #set yrange [-2:2]
  plot file1 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle,\
  file2 binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:5 with linespoint pointtype 7 lc rgb "red" notitle
  unset xlabel
  unset ylabel
  unset yrange
  

  unset title
  unset multiplot
  unset output
  
  
}


exit








 




#################
# yz diagonal cut
#################



do for [j=0:nt/tskip -1] {

  set terminal png size 2560,1440
  
  i = j*tskip
  
  filename="Frames/fields_yz_diagcut_t=".i.".png"
  set output filename
  
  set title "Time Step = #".(i+1)
  set multiplot layout 2,3 
  
  file="Snapshots/yzdiagcut_rank=000_t=".i.".dat"
  file2="Snapshots/yzdiagcut_rank=001_t=".i.".dat"
  
  i = j
  
  set xlabel "y=z"
  set ylabel "RHO"
  #set yrange [-0.1:1.8]
  plot file binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle,\

  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "y=z"
  set ylabel "Bpar"
  #set yrange [-2:2]
  plot file binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:7 with linespoint pointtype 7 lc rgb "blue" notitle
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "y=z"
  set ylabel "Bx"
  #set yrange [-2:2]
  plot file binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:8 with linespoint pointtype 7 lc rgb "blue" notitle
  unset xlabel
  unset ylabel
  unset yrange

  set xlabel "y=z"
  set ylabel "vperp"
  #set yrange [-2:2]
  plot file binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
  set xlabel "y=z"
  set ylabel "vpar"
  #set yrange [-2:2]
  plot file binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:4 with linespoint pointtype 7 lc rgb "blue" notitle
  unset xlabel
  unset ylabel
  unset yrange
 
  set xlabel "y=z"
  set ylabel "vx"
  #set yrange [-2:2]
  plot file binary format="%*1int%double%double%double%double%double%double%double%double%*1int" using 1:5 with linespoint pointtype 7 lc rgb "blue" notitle
  unset xlabel
  unset ylabel
  unset yrange
  

  unset title
  unset multiplot
  unset output
  
  
}


exit
