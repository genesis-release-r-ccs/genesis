# read pdb and select molecule
mol load pdb pdb/dna_original.pdb
set all [atomselect top all]

# remove center of mass
$all moveby [vecinvert [measure center $all weight mass]]

# save chains separately
foreach chain {A B} {
    [atomselect top "chain $chain"] writepdb dna_$chain.pdb
}

# generate psf
package require psfgen
resetpsf
topology ../../param/top_all27_na.rtf
foreach chain {A B} {
    set sel [atomselect top "chain $chain and name C1'"]
    set seg DNA$chain
    # create segments apply the patches for terminals
    segment $seg {
        first 5TER
        last 3TER
        pdb dna_$chain.pdb
    }
    # apply the patched DEO1 for pyrimidines and DEO2 for purines
    # otherwise the resulting psf becomes RNA not DNA
    foreach resid [$sel get resid] resname [$sel get resname] {
        if { $resname eq "THY" || $resname eq "CYT" } {
            patch DEO1 $seg:$resid
        } elseif { $resname eq "ADE" || $resname eq "GUA" } {
            patch DEO2 $seg:$resid
        }
    }
    coordpdb dna_$chain.pdb $seg
}

# guess coordinates of missing hydrogens
guesscoord

# write
writepsf dna_only.psf
writepdb dna_only.pdb

# solvation
mol delete all
package require solvate
solvate dna_only.psf dna_only.pdb -minmax {{-31.8 -31.8 -31.8} {31.8 31.8 31.8}} -o dna_solv

# add counter ions
mol delete all
package require autoionize
autoionize -psf dna_solv.psf -pdb dna_solv.pdb -neutralize -cation SOD -o dna_ion

# check boxsize
mol delete all
mol new dna_ion.psf
mol addfile dna_ion.pdb
set all [atomselect top all]
set minmax [measure minmax $all]
foreach {min max} $minmax { break }
foreach {xmin ymin zmin} $min { break }
foreach {xmax ymax zmax} $max { break }
puts "cellBasisVector1 is: "
puts "{[expr abs($xmin)+ $xmax + 0.3] 0.0 0.0}"
puts "cellBasisVector2 is: "
puts "{0.0 [expr abs($ymin)+ $ymax + 0.3] 0.0}"
puts "cellBasisVector3 is: "
puts "{0.0 0.0 [expr abs($zmin)+ $zmax + 0.3]}"
puts "the center of the box is (rough estimation): "
puts "{[expr ($xmin + $xmax)/2.0] [expr ($ymin + $ymax)/2.0] [expr ($zmin + $zmax)/2.0]}"

# make a constraint reference file
# B-value will be replaced with force constants of restraint
mol delete all
mol new dna_ion.psf
mol addfile dna_ion.pdb
set all [atomselect top "nucleic"]
set fixed 1.0
set nofixed 0
set fixed_list {}
foreach i [$all get "nucleic"] {
    if {$i == 1} {
        lappend fixed_list $fixed
    } else {
        lappend fixed_list $nofixed
    }
}
$all set beta $fixed_list
$all writepdb constraint.pdb

# make a restraint reference file
# B-value will be replaced with force constants of restraint
mol delete all
mol new dna_ion.psf
mol addfile dna_ion.pdb
set all [atomselect top "all"]
set kforce 10.0
set kforce_list {}
foreach i [$all get "backbone"] {
    if {$i == 1} {
        lappend kforce_list $kforce
    } else {
        lappend kforce_list 0.0
    }
}
$all set beta $kforce_list
$all writepdb restraint.pdb

# make both charmm-type and xplor-type psf files
mol delete all
resetpsf
#topology ../../param/top_all27_na.rtf
readpsf dna_ion.psf
coordpdb dna_ion.pdb

writepsf x-plor dna_xplor.psf
writepsf charmm dna_charmm.psf

