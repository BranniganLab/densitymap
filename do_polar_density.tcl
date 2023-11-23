source ./TCL/polarDensity.tcl

set HeadNames [atomselect top "name PO4"] ;#B2 is DDM
set lipids [lsort -unique [$HeadNames get resname]]
puts $lipids
$HeadNames delete
puts $lipids
set Rmax 40.
set Rmin 5.
set dr 2.
set Ntheta 30
set dt 1
set sample_frame 1
foreach lip $lipids {
	polarDensityBin $lip.dat $lip $Rmin $Rmax $dr $Ntheta $dt $sample_frame	
}
