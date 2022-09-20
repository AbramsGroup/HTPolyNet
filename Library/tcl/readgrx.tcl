# This proc returns a TcL dictionary keyed by attribute name with values that are atom-by-atom
# Source this file in your vmd script or a TcL interactive session
# IMPORTANT: If you are using a TcL interactive session, you must SQUASH the output echo or vmd will HANG!
#
# One way to do this is to chain the "set" command with creation of an empty list:
#
# % set D [get_grx "myfile.grx"]; list
#
# That little "list" command chained after the set squashes the console output of the set command
# You do not want to see the console outuput of this set command if your system has more than
# a few thousand atoms, trust me.
#
# But having this dictionary is useful.  Say you want to color your system by molecule index.  You can use 
# any editable tag in the atom selection, like "beta":
#
# % set a [atomselect top all]
# % set D [get_grx "myfile.grx"]; list
# % $a set beta [dict get $D "molecule"]
#
# Cameron F. Abrams cfa22@drexel.edu

proc get_grx { fn } {
    set D [dict create]
    set fp [open $fn "r"]
    set data [read -nonewline $fp]
    close $fp
    set lines [split $data "\n"]
    set nlines [llength $lines]
    set length_of_last_line [string length [lindex $lines end]]
    while { $length_of_last_line == 0 } {
        set lines [lrange $lines 0 end-1]
        set length_of_last_line [string length [lindex $lines end]]
    }
    set hdr [lindex $lines 0]
    set lines [lrange $lines 1 end]
    set ndatalines [llength $lines]
    set cols [regexp -all -inline {\S+} $hdr]
    for {set i 0} {$i < [llength $cols]} {incr i} {
        set mylist {}
        foreach {ln} $lines {
            lappend mylist [lindex $ln $i]
        }
        set D [dict append D [lindex $cols $i] $mylist]
    }
    return $D
}