# This script makes still images of systems that include only residues that
# have reacted
#
# 
# inside of a "systems" subdirectory, issue "vmd -dispdev text -e render.tcl"
# Cameron F Abrams

package require pbctools

# read a column-oriented data file with a header row
# into a dictionary keyed by column name with value list
# from that column
# Assumes rows are whitespace-delimited
proc read_cod { filename } {
    set pfi [open $filename "r"]
    # read header
    set cnt [gets $pfi header]
    if {$cnt < 0} {return}
    set columns [regexp -all -inline {\S+} $header]
    set D [dict create]
    foreach h $columns {
        dict append D $h [list]
    }
    while {1 == 1} {
        set cnt [gets $pfi line]
        if {$cnt < 0} {break}
        set columns [regexp -all -inline {\S+} $line]
        for {set i 0} {$i < [llength $columns]} { incr i } {
            set key [lindex $header $i]
            dict lappend D $key [lindex $columns $i]
        }
    }
    close $pfi
    return $D
}

proc render_by_reactants { iter_no image_file current_resids } {
    set pre "iter-"
    append pre $iter_no
    set dfn $pre
    append dfn "/2-cure_update-bonds.csv"
    set cfg $pre
    append cfg "/4-cure_equilibrate-npt.gro"
    set D [read_cod $dfn]
    set R [lsort -integer -unique [concat $current_resids [dict get $D ri] [dict get $D rj]]]

    set viewplist {}
    mol new $cfg type gro first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
    mol delrep 0 top
    mol representation Licorice 0.300000 12.000000 12.000000
    mol color ResName
    mol selection "resid $R"
    mol material Opaque
    mol addrep top
    
    # set viewpoints {{{1 0 0 -25.9876} {0 1 0 -26.254} {0 0 1 -26.024} {0 0 0 1}} {{0.779795 -0.0110901 0.625937 0} {0.0650928 0.99586 -0.0634485 0} {-0.622642 0.0902208 0.777288 0} {0 0 0 1}} {{0.0264722 0 0 0} {0 0.0264722 0 0} {0 0 0.0264722 0} {0 0 0 1}} {{1 0 0 -0.06} {0 1 0 0.07} {0 0 1 0} {0 0 0 1}}}
    # molinfo top set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints
    
    axes location off
    color Display {Background} white
    pbc box
    render TachyonInternal $image_file
    mol delete top
    return $R
}

set current_resids [list]
for {set i 1} { $i < 36 } { incr i } {
  set imname $i 
  append imname ".tga"
  set current_resids [render_by_reactants $i $imname $current_resids]
}
quit

