reset
# set terminal pdfcairo dashed enhanced colour size 7cm,7cm \
# 	font "Helvetica,12"
set terminal pngcairo dashed enhanced colour size 20cm,23cm \
	font "Helvetica,10"
	
set output outputfilename

set pm3d map

set multiplot layout 2,3

set tmarg 0
set bmarg 0
set lmarg 2
set rmarg 0
set size 1.12, 0.95
set origin -0.055, 0.15

set xrange [0:gridsize-1]
unset xtics 

set yrange [0:gridsize-1]
unset ytics

set cbrange[mincb:maxcb]

unset xlabel
unset ylabel
set format x ""
set format y ""

#set palette rgb 10,13,33
#set palette model HSV rgbformulae 3,2,2
#set palette model XYZ rgbformulae 7,5,15
#set palette model RGB defined (0 "brown", 1 "brown", 1 "orange", 2 "orange", 2 "yellow", 3 "yellow" )
#set palette model RGB defined (0 "white", 2 "white", 5 "blue", 10 "blue", 20 "black", 30 "brown", 50 "orange", 59 "yellow", 60 "green")

set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000', 8 '#7f0000')

# set palette defined ( 0 '#7f0000',1 '#ee0000',2 '#ff7000',3 '#ffee00',4 '#90ff70',5 '#0fffee',6 '#0090ff',7 '#000fff', 8 '#000090')

# set palette defined (0 '#0fffee', 1 '#0090ff', 2 '#000fff', 3 'black', 3 'black', 4 '#ee0000', 5 '#ff7000', 6 '#ffee00', 7 'white', 7 'white', 8 '#0fffee', 9 '#0090ff', 10 '#000fff', 11 'black', 11 'black', 12 '#ee0000', 13 '#ff7000', 14 '#ffee00')

set cblabel clabelname

set colorbox horizontal back user origin 0.065, 0.15 size 0.86, 0.05

splot filename u 1:2:(log10(1.-$3*column))

unset multiplot
