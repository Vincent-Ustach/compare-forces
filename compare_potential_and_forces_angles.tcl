#
# print the maximal error in the numerical derivation of a potential
# against the displacement th. This is a central derivate, so we
# expect at least a dth^3 scaling (???), if the forces are continuous.
#
# This can be used to check whether the forces and energy functions of
# a potential fit.
#
#######################################################################
#
# Copyright (C) 2010,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 

#Modified by Vincent Ustach July 24 2014

#for some reason the system would not recognize veccross unless I added it here!
proc veccross {b1 b2} {
    lassign $b1 b11 b12 b13
    lassign $b2 b21 b22 b23

    set b31 [expr $b12*$b23 - $b13*$b22]
    set b32 [expr $b13*$b21 - $b11*$b23]
    set b33 [expr $b11*$b22 - $b12*$b21]
    return [list $b31 $b32 $b33 ]
}
	set data [open "data.file" "w"]
# box size, should be larger than the potentials cutoff
set box_l 10.0
setmd box_l $box_l $box_l $box_l

# number of interaction partners, 2 for nonbonded
set n_part 3

# create the particles so that you can create bonds
for {set i 0} {$i < $n_part} {incr i} {
    part $i pos 0 0 0
}

# The potential to test, including parameters.
#
# - can be bonded or nonbonded
# - all particles have type 0
#inter 1 harmonic 10 0.5 2
inter 1 angle_harmonic 10 [expr atan(1)*2.0]
#inter 1 tabulated bond "./harmonic.txt"
#inter 1 tabulated angle "./harmonic_angle.txt"
part 1 bond 1 0 2

# minimal distance of the particles, to avoid areas with large forces
# where the convergence is bad
set min_dist 0

# particles are always in a sphere of this diameter
# choose bigger than the cutoff to test the shifting
set max_dist 1.0

# number of tests
set tries 10000

# below you typically don't need to change anything
##################################################

setmd periodic 1 1 1
thermostat off
setmd time_step 1
setmd skin 0

puts "# dth maximal deviation"
puts "# should scale like dth^3 for sufficiently small dth??"
for {set dth 4.0} {$dth > 0.0001} {set dth [expr $dth*0.5]} {
    set maxDev 0
    for {set r 0} {$r < $tries} {incr r} {
	# 1. position
	while {1} {
	    # particles in a sphere around $bp
	    set bp "0.0 0.0 0.0"
	    #set bp "[expr $box_l*rand()] [expr $box_l*rand()] [expr $box_l*rand()]"
	    for {set i 0} {$i < $n_part} {incr i} {
		while {1} {
		    set vec ""
		    for {set j 0} {$j<3} {incr j} {
			lappend vec [expr (2.0*rand() - 1.0)*$max_dist/2] }
		    if { [veclen $vec] <= $max_dist/2 } { break }
		}
		set p1($i) [vecadd $bp $vec]
		eval part $i pos $p1($i)
	    }
	    if {[analyze mindist] > $min_dist} { break }
	}

	set e1 [analyze energy total]
	integrate 0
	for {set i 0} {$i < $n_part} {incr i} {
	    set f1($i) [part $i print force]
	}

	# 2., slightly displaced position for particle 0
	#Rodrigues rotation formula
	#vector k normal to the plane containing the three beads is the axis of rotation if the inner bead is at the origin
	#vector k is defined by the cross product of two lines within the plane containing the coords of all three beads 
        set k [ veccross [vecadd $p1(2) -$p1(0)] [vecadd $p1(1) -$p1(0) ] ] 
		set first [vecscale [expr cos($dth)] $p1(0)]
		set second [vecscale [expr sin($dth)] [veccross $k $p1(0)  ] ]
		set factor [expr [vecdot_product $k $p1(0) ] * [expr 1.0 - cos($dth)] ]
		set third [vecscale $factor $k] 
	    set p2(0) [vecadd $first [vecadd $second $third]] 
	    eval part 0 pos $p2(0)


	set e2 [analyze energy total]
	integrate 0
	for {set i 0} {$i < $n_part} {incr i} {
	    set f2($i) [part $i print force]
	}

	# Check forces
	set dE [expr $e2 - $e1]

	set fdotdx 0
	for {set i 0} {$i < $n_part} {incr i} {
	    # average force
	    set force [vecscale 0.5 [vecadd $f1($i) $f2($i)]]
	    # no need to dot the force with with the change in direction because the force already points in that direction!
	    set fdotx [expr $fdotdx + [veclen $force] ]
	}
	puts $data "$r $p1(0) $p2(0) $dE $fdotx"
	set dev [expr abs($dE + $fdotdx)]
	if {$dev > $maxDev} { set maxDev $dev }
    }
    puts "$dth $maxDev"
}

exit
