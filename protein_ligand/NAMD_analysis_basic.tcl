#!/usr/bin/tclsh

# ########################################################################## #
# This script runs NAMD basic analysis.                                      #
# It calculates rmsd_ptn, rmsd_lig, rmsf_ptn, rg_ptn, hbonds.                #
# It requires psf and dcd files.                                             #
# You can change the file names in psf and dcd variables at lines 16 and 17. #
# ########################################################################## #

#if { $argc != 2 } {
#    puts "This script requires two files to be inputed: <psf> and <dcd> files"
#    puts "For example, source new1.tcl protein.psf trajectory.dcd"
#    puts "Please try again"
#}

set psf "complex_wb_neutralized.psf"
set dcd "out_eq_npt_01.dcd"
 
# Wrap
set mol [mol new $psf type psf waitfor all]
#set mol [mol new [lindex $argv 1] type psf waitfor all]

mol addfile $dcd type dcd step 1 waitfor all 
#mol addfile [lindex $argv 2] type dcd step 10 waitfor all 

pbc wrap -centersel "protein" -center com -compound fragment -sel "protein or resname LIG" -all

# RMSD_ptn

set ref [atomselect $mol "name CA" frame 0]
set all [atomselect $mol all]
set sel [atomselect $mol "name CA"]


set nframes [molinfo top get numframes]

proc myrmsd { frame } {
	global ref sel all
	$all move [measure fit $sel $ref]
	set fout [open "rmsd_ptn.csv" a+]
	puts "$frame: [measure rmsd $sel $ref]"
	puts $fout [format "%d %s %f" $frame , [measure rmsd $sel $ref]]
close $fout
}



set fout [open "rmsd_ptn.csv" w]
	puts $fout "frame, rmsd"
close $fout

for {set i 0} {$i <$nframes} {incr i} {
	animate goto $i
	myrmsd $i
}

# RMSD_lig

set ref [atomselect $mol "resname LIG" frame 0]
set all [atomselect $mol all]
set sel [atomselect $mol "resname LIG"]

set nframes [molinfo top get numframes]

proc rmsd_lig { frame } {
	global ref sel all
	set fout [open "rmsd_lig.csv" a+]
	puts "$frame: [measure rmsd $sel $ref]"
	puts $fout [format "%d %s %f" $frame , [measure rmsd $sel $ref]]
close $fout
}

set fout [open "rmsd_lig.csv" w]
	puts $fout "frame, rmsd"
close $fout

for {set i 0} {$i <$nframes} {incr i} {
	animate goto $i
	rmsd_lig $i
}


# RMSF
set num [expr {$nframes - 1}]

set fout [open "rmsf_ptn.csv" w]
	puts $fout "residue, rmsf"

set sel [atomselect top "name CA"]
set rmsf [measure rmsf $sel first 0 last $num step 1]
for {set i 0} {$i < [$sel num]} {incr i} {
#  puts $fout [format "%d %s %f" $frame , $rmsf]
  puts $fout "[expr {$i+1}], [lindex $rmsf $i]"
#  puts $fout [format "%d %s %f" [expr{$i+1}] , [lindex $rmsf $i]]

} 
close $fout

# Rog
proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
	#  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}

proc gyr_radius {sel} {
  # make sure this is a proper selection and has atoms
  if {[$sel num] <= 0} {
    error "gyr_radius: must have at least one atom in selection"
  }
  # gyration is sqrt( sum((r(i) - r(center_of_mass))^2) / N)
  set com [center_of_mass $sel]
  set sum 0
  foreach coord [$sel get {x y z}] {
    set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
  }
  return [expr sqrt($sum / ([$sel num] + 0.0))]
}

set outfile [open "rg_ptn.csv" w]
puts $outfile "frame, rad_of_gyr"
set nf [molinfo top get numframes] 
set i 0

set prot [atomselect top "protein"] 
while {$i < $nf} {

    $prot frame $i
    $prot update

    set i [expr {$i + 1}]
    set rog [gyr_radius $prot]

    puts $outfile "$i, $rog"

} 

close $outfile

# H-bonds
package require hbonds
hbonds -sel1 [atomselect top protein] -sel2 [atomselect top "resname LIG"] -writefile yes -plot no

set fout [open hbonds.csv w]
	puts $fout "frame, hbonds"

set f [open "hbonds.dat"]
while {[gets $f line] >= 0} {
    set new_line [regsub -all {\s+} $line {, } ]
    puts $fout "$new_line"
}
close $f

close $fout
