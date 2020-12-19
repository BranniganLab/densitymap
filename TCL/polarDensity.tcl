
 package require pbctools
 set UTILS "/u2/home_u2/lms464/github/JPC_Special/common/utils" 
 set QWRAP "/u1/home/lms464/lms464/github/qwrap"
source $UTILS/BinTools.tcl
load ${QWRAP}/qwrap.so
#
#  Developed using Martini CG'ed membranes
# 
# ;#Lipid_Saturation_HeadG are a series of macros to parse Martini lipids
# source ${UTILS}/Lipid_Saturation_HeadG.tcl


#Grace Brannigan 7/2018

#Sample Use:
#Assuming trajectory of interest is loaded AND top.  All frames will be used so any frames to be ignored would need to be unloaded
#set HeadNames [atomselect top "name PO4 ROH B2"] ;#B2 is DDM
#set lipids [lsort -unique [$HeadNames get resname]]
#$HeadNames delete
#set RMax 40.
#set RMin 5.
#set dr 2.
#set Ntheta 30
#foreach lip in $lipids {
#polarDensityBin $lip.dat $lip $Rmin $Rmax $dr $Ntheta	
#}

###############################################Sample plotting within python
#
#Ntheta = 30
#data = np.loadtxt('DPPC.dat',skiprows=2)
#rad = data[:,1] + (data[:,1]-data[:,0])/2.0
#the = np.linspace(0,2*np.pi,Ntheta +1)
#theta,radius=np.meshgrid(the,rad)
#density = data[:,3:]/radius 
#plt.figure(figsize = (5,5))
#plt.subplot(projection="polar")
#plt.pcolormesh(theta,radius,density,cmap="RdBu",zorder=0,edgecolors='k',lw=.001)
#plt.show()

proc get_avg_area {molid} {
	;#set nframes [molinfo $molid get numframes]
    set box [pbc get -all]
    set xbox [list]
    set ybox [list]
    foreach r $box {
        lappend xbox [lindex $r 0]
        lappend ybox [lindex $r 1]
    }
    set x [expr 1.0 * [vecsum $xbox] / [llength $xbox]]
    set y [expr 1.0 * [vecsum $ybox] / [llength $ybox]]
	set avg [expr 1.0 * $x * $y]
	return $avg
}

;# histograms lists of data
proc lcount list {
    foreach x $list {lappend arr($x) {}}
    set res1 {}	
	set res2 {}
    foreach name [array names arr] {
		lappend res1 $name
		lappend res2 [llength $arr($name)]
    }
	set res [list $res1 $res2]
    return $res
 }
 
 proc RtoD {r} {
    set pi 3.14159265358979323846
    return [expr $r*180.0/$pi]
}


proc get_theta {x y} {
    set pi 3.14159265358979323846
    set tmp  [expr {atan2($y,$x)}]
    if {$tmp < 0} {
        set theta [expr 2*$pi + $tmp]    
    } else {
        set theta $tmp
    }
    return [RtoD $theta]
}

;# sums over a list
proc Sum_list {list_in} {
    set list_out 0
    foreach li $list_in {
        set list_out [expr 1.0*$list_out+$li]
    }
    return $list_out
}

;# calculates average center of mass of a flat membrane
proc z_mid {init_frm nframes} {
    set z_list {}
    for {set frm ${init_frm}} {${frm} < ${nframes}} {incr frm} {
        set mid [atomselect top "name PO4 ROH C3 PO41" frame $frm]
        lappend z_list [lindex [measure center $mid weight mass] 2]
        $mid delete
    }
    return [expr 1.0*[vecsum $z_list]/([llength $z_list]) ]
}

;# Ouputs position of the centered protein in a membrane
;# accross both leaflets
proc Protein_Position {{a ""}} {
	;# list for the chain names
    set chain_names [list "A" "B" "C" "D" "E"]
    ;# finds the center of the membranes
    set zed [z_mid 0 20]
	set occupancy [list 1 2 3 4]
	;# calculates the center of mass for subunit alpha helices in both leaflets
	foreach eq {"<" ">"} eqtxt {"lwr" "upr"} {
		set fout [open "/u2/home_u2/lms464/github/JPC_Special/tasks/17_Aim1/Data/Protein${a}_coords_${eqtxt}.dat" w]
        puts $fout  "# chain A ooc1r occ1the occ2r occ2the... "
        foreach chnm $chain_names {
            foreach occ $occupancy {
                set sel [atomselect top "(chain ${chnm}) and (occupancy $occ and name BB) and (z ${eq} $zed)" frame 0]
                set com [measure center $sel weight mass]
                $sel delete
                set x [lindex $com 0]
                set y [lindex $com 1]
                set r [expr sqrt($x*$x+$y*$y)]
                set theta [get_theta $x $y]
                puts "chain ${chnm} and occupancy $occ $r $theta"

                puts -nonewline $fout "$r $theta "
            }
            puts $fout ""
        }
        close $fout
    }
}
;# calculates the average chain length per lipid species
proc avg_acyl_chain_len {species} {
    
    set acyl_num 0
    set sel [atomselect top "$species"]
    set sel_resname [lsort -unique [$sel get resname]]
    $sel delete
    ;# if a group of lipids are selected (only PUFAs), 
    ;# determines which lipid species are there (ie,
    ;# PUPS, DOPC....)
    foreach res $sel_resname {
        set sel [atomselect top "${species} and (resname $res) and (not name NH3 NC3 GL1 GL2 AM1 AM2 PO4 CNO CN0 C1 C2 C3)"]
        set sel_len [llength [lsort -unique [$sel get name]]]

        ;# signals counting full lipids not acyl chains
        if {$sel_len > 6} {
            lappend sel_resname "${res}"
        }

        $sel delete
        set acyl_num [expr $acyl_num + $sel_len]
    }
    ;# calculates average chain lengths
    set avg_acyl_chain [ expr (1.0 * $acyl_num / [llength $sel_resname]) ]
    ;# no lipids, "error" return
    if {$avg_acyl_chain < 1} {
        return 1
    }
    return $avg_acyl_chain

}

proc Center_System {inpt} {
    puts "${inpt}"
    ;# confirms your box is either square or paraelleogram-ish
    ;# will preform qwrap or pbc wrap depending
    ;# TODO, does not properly scan to confirm
    ;#		protein is centered for every SINGLE frame
    ;#		Coarse Grained proteins that move accross
    ;#		PBC tend to have an issue and need to be
    ;#		centered mulitple times
    set pbc_angles [molinfo top get {alpha beta gamma}]
    
    set sel [atomselect top "$inpt" frame 0]
    set com [measure center $sel weight mass]
    
    set counter_i 0
    ;# continues to try and recenter's box until ~ 0'ed out
    while {[expr abs([lindex $com 0])] > 1.0 &&  [expr abs([lindex $com 1])] > 1.0} {
        
        if {$counter_i > 5} {
            puts "Script was unable to converge system to (0,0,0)"
            puts "Please check your system vissually, there may be"
            puts "unintended artifacts"
            $sel delete
            return
        }
        
        if {([lindex $pbc_angles 0]!=90.0) && ([lindex $pbc_angles 0]!=90.0) && ([lindex $pbc_angles 0]!=90.0)} {
            puts "qwrap may not be optimal for your system...\n"
            puts "Running pbc wrap. To verify proper centering"
            puts "pbc wrap will be run multiple times" ; after 100
            foreach i {0 1 2 3} {
                pbc wrap -centersel "$inpt" -all
            }
        } else {
            qwrap centersel "$inpt" ;#center entire system at ~0,0,0
        }
        set com [measure center $sel weight mass]
        incr counter_i
    }
    $sel delete
}

proc resnamer {input} {

    ;# adds resname if the input is DPPC, CHOL, PUPI...
    set out ""
    if {[string length $input] == 4 && $input != "chol"} { 
        set out "resname $input"
    } else {
        set out "$input"
    }
    return $out
}

;# writes line to file
proc output_bins {fl  ri rf dtheta bins} {
    puts -nonewline $fl "[format {%0.2f} $ri]  [format {%0.2f} $rf] [format {%0.2f} $dtheta]  " 
    puts $fl "$bins" 
}

;# runs through a bin over desired frames
proc bin_over_frames {shell species dtheta sample_frame nframes Ntheta dt ri rf fupper flower} {

    set theta_bin_high [lrepeat [expr $Ntheta+1] 0]
    set theta_bin_low [lrepeat [expr $Ntheta+1] 0]
    for {set frm $sample_frame} {$frm < ${nframes}} {incr frm $dt} {

        $shell frame $frm
        $shell update 

        set singleFrame_counts [bin_frame $shell $species $dtheta $frm ]

        set singleFrame_upper [lindex $singleFrame_counts 0] 

        set singleFrame_lower [lindex $singleFrame_counts 1]
        set theta_bins [theta_histogram $singleFrame_upper $singleFrame_lower $Ntheta]

        if { [llength $theta_bin_high] != [llength [lindex $theta_bins 0]] } {
            error "theta_bin_high/low and theta_bins do not have the same length."
        }
        set theta_bin_high [vecadd $theta_bin_high [lindex $theta_bins 0] ]
	    set theta_bin_low [vecadd $theta_bin_low [lindex $theta_bins 1]]

    }
  return [list ${theta_bin_high} ${theta_bin_low}]
}


;# Determines if the lipid is in the outer or iner leaflet,
;#  based on the lipid resid
proc local_mid_plane {atsel_in frame_i} {


    set sel_resid [atomselect top "$atsel_in" frame $frame_i]
    set ind 1
    if { [string range [lsort -unique [$sel_resid get resname]] end-1 end] == "PA" } {
    	set ind 0
    }
    set sel_Z [${sel_resid} get z] 
    $sel_resid delete
	if {[lindex ${sel_Z} $ind] < [lindex ${sel_Z} end] } { 
		return 1 
	} else { 
		return 0
	}
}



;# Collects data for a bin in a single frame
proc bin_frame {shell species dtheta frm } {
    set indexs [$shell get index]
    set resids [$shell get resid]
    set nShell [$shell num]
    set theta_high_out [list]
    set theta_low_out [list]
    set resd_old 0
    set high_low 0
    foreach indx $indexs resd $resids {
        #loop over lipids in the shell
        ;# selection to count
        set a "($species and index $indx)"
        ;# selection to determine leaflet
        set b "(resid $resd)" 
        set thislipid [atomselect top $a frame $frm]
        if {[string length ${species}] == 2} {
	    	if {[$thislipid get name] == "PO4"} {
	        	continue
	    	}
    	}
    	;# don't update
        if {${resd_old} != ${resd}} {
        	set high_low [local_mid_plane $b  $frm]
        }
        set x [$thislipid get x]
        set y [$thislipid get y]
        $thislipid delete
        set theta [get_theta $x $y]
        set ti [expr int($theta/$dtheta)] 
        #determine theta bin
        if {$high_low > 0} {
            lappend theta_low_out $ti
        } else {
            lappend theta_high_out $ti
        }

    }
    
    return [list $theta_high_out $theta_low_out] 
}

;# Builds histogram over theta bins
proc theta_histogram {singleFrame_upper singleFrame_lower Ntheta } {
    
    set theta_bin_out [list]

    foreach ud [list $singleFrame_upper $singleFrame_lower] {
        #cleanup and output 
        set theta_bin_counts [lcount $ud]
        #Shell_Test $shel_count $theta_bin_counts
        set theta_bins {}
        ;# make this into the new lcount? better Idea TEST lcount
        for {set ti 0} { $ti<=$Ntheta} {incr ti 1} {
            set tindex [lsearch [lindex $theta_bin_counts 0]  $ti]
            if { $tindex >= 0} {
                set frame_count [expr 1.0 * [lindex [lindex $theta_bin_counts 1] $tindex]] 
            } else { 
                set frame_count 0.0
            }
            lappend theta_bins $frame_count
        }
        lappend theta_bin_out $theta_bins
    }
    return $theta_bin_out
}


;########################################################################################
;# polarDensity Funciton

proc polarDensityBin { outfile species Rmin Rmax dr Ntheta} {

	;# funciton that sets CG'ed chains and selects alpha helecies
	source /u2/home_u2/lms464/github/JPC_Special/common/grace/assign_helices_2BG9_CG_lms2.tcl

	;# see resnamer
    set species [resnamer ${species}]
    
    set sel [atomselect top "$species"]
	set sel_num [$sel num]
	
	if {$sel_num == 0} {
        error "No lipid saturation set exists"
	}
	
	;# Center's system (weak hack)
 	Center_System "occupancy 1 to 4 and name BB"
    Center_System "occupancy 1 to 4 and name BB"
    Center_System "occupancy 1 to 4 and name BB"
    ;# aligns protein
 	Align "occupancy 1 to 4 and name BB"
 	;# outputs protein positions
    Protein_Position   
    ;# initialize some constants
	set dt 1
	set nframes [molinfo top get numframes]
  	set sample_frame 100
	set delta_frame [expr ($nframes - $sample_frame) / $dt]
	set area [get_avg_area top]
    $sel delete
	puts "Acyl Chain:\t$species"
	set low_f [open "${outfile}.low.dat" w]
    set upp_f [open "${outfile}.upp.dat" w]
	set dtheta [expr 360.0/(1.0*($Ntheta))]

	;# builds header for output file (extimates expected density, membrane area...)	
    foreach lu [list $low_f $upp_f] zed [list "(z<0)" "(z>0)"] {
        set sel [ atomselect top "(($species) and $zed) and (name PO4 ROH)"  frame 0]
        set sel_num [llength [lsort -unique [$sel get resid] ] ]
        if {$sel_num < 1} {
            set num_beads 0
        	set expected 0
        } else {
	        set sel_resid [lsort -unique [$sel get resid] ]
	        $sel delete
	        set beads [atomselect top "${species} and (resid $sel_resid) and (not name PO4)" frame 0]
	        set num_beads [$beads num]
	        set expected [expr 1.0 * $num_beads/$area]
	        $beads delete
    	}
        puts "#Lipid species $species : ${sel_num} molecules, Num beads : ${num_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Density : [format {%0.5f} [expr $expected]]/A^2, Average Chain : [avg_acyl_chain_len ${species}] beads, dr*dtheta : [format {%0.5f} [expr $dr*[DtoR $dtheta]]] "
	    puts $lu "#Lipid species $species : ${sel_num} molecules, Num beads : ${num_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Density : [format {%0.5f} [expr $expected]]/A^2, Average Chain : [avg_acyl_chain_len ${species}] beads, dr*dtheta : [format {%0.5f} [expr $dr*[DtoR $dtheta]]] "
    }

	#loop over radial shells
	for {set ri $Rmin} { $ri<=${Rmax}} { set ri [expr $ri + $dr]} {
		puts "Ring {$ri [expr ${ri}+${dr}]}"
		set rf [expr $ri + $dr]
		set rf2 [expr $rf*$rf]
		set ri2 [expr $ri*$ri]
		set shell [atomselect top "($species) and ((x*x + y*y < $rf2) and  (x*x + y*y > $ri2)) and (not name PO4)"]
		#selects lipids in the radial shell
		
        set theta_bin [bin_over_frames $shell $species $dtheta $sample_frame $nframes $Ntheta $dt $ri $rf $upp_f $low_f]
        set theta_bin_high [lindex $theta_bin 0]
        set theta_bin_low [lindex $theta_bin 1]

        $shell delete	
        set time_avg_upper [vecscale $theta_bin_high [expr 1.0 / (1.0 * $delta_frame)]]
        set time_avg_lower [vecscale $theta_bin_low [expr 1.0 / (1.0 * $delta_frame)]]
        output_bins $upp_f $ri $rf $dtheta "$time_avg_upper" 
        output_bins $low_f $ri $rf $dtheta "$time_avg_lower" 
	}
	close $low_f
	close $upp_f
}
