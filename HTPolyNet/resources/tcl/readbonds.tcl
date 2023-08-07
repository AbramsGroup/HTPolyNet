
# This script is intended to be loaded in an interactive Tcl session in
# VMD.  If "fn" is set to the name of a bonds CSV file (see example),
# then this will generate the two lists "ai" and "aj" that contain
# the serial numbers of each atom in each bond.  The loop at the
# end just draws them.
#
# Cameron F Abrams cfa22@drexel.edu

set fn "proj-0/systems/init-2/2-cure_update-bonds.csv"
set fp [open $fn "r"]
set hdr [gets $fp]
set cols [regexp -all -inline {\S+} $hdr]
puts "$cols"
set idx [lsearch $hdr "ai"]
set jdx [lsearch $hdr "aj"]
# puts "$idx"
set ais {}
set ajs {}
while {[gets $fp line]>=0} {
    set tok [regexp -all -inline {\S+} $line]
    set ai [lindex $tok $idx]
    lappend ais $ai
    set aj [lindex $tok $jdx]
    lappend ajs $aj
}
close $fp
foreach i $ai j $aj {
    set ri [lindex [[[atomselect top "serial $ai"] get {x y z}]] 0]
    set rj [lindex [[[atomselect top "serial $aj"] get {x y z}]] 0]
    graphics 0 cylinder $ri $rj diameter 0.2
}
