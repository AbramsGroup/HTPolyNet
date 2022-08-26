# This script is intended to be sourced in an interactive Tcl session
# in VMD.  If "fn" is set to a grx file, it reads this file in and
# generates the "molname" and "molid" lists for each atom.  These
# can then be used to set beta or occupancy or whatever in a
# full-system atom selection to allow "color by" to color
# by molecule name or molecule id.
#
# Cameron F. Abrams cfa22@drexel.edu

#set fn "proj-0/systems/final-results/final.grx"
set fn "proj-0/systems/init/init.grx"
set fp [open $fn "r"]
set hdr [gets $fp]
set cols [regexp -all -inline {\S+} $hdr]
puts "$cols"
set idx [lsearch $hdr "molecule_name"]
set jdx [lsearch $hdr "molecule"]
set molnames {}
set molids {}
while {[gets $fp line]>=0} {
    set tok [regexp -all -inline {\S+} $line]
    set molname [lindex $tok $idx]
    lappend molnames $molname
    set molid [lindex $tok $jdx]
    lappend molids $molid
}
close $fp
