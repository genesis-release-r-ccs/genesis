# read pdb
mol load pdb pdb/256l_edited.pdb

# remove center of mass
set all [atomselect top all]
$all moveby [vecinvert [measure center $all weight mass]]

# select water and save pdb
#set water [atomselect top water]
set water [atomselect top "water and (within 5 of protein)"]
$water writepdb pdb/256l_water.pdb

# select protein and save pdb
set protein [atomselect top protein]
$protein writepdb pdb/256l_protein.pdb

# generate psf
package require psfgen
resetpsf
topology ../../param/top_all27_prot_lipid.rtf

pdbalias residue HIS HSE
segment T4 {
 first nter
 last cter
 pdb pdb/256l_protein.pdb
}

pdbalias atom ILE CD1 CD
coordpdb pdb/256l_protein.pdb T4

guesscoord

# write
writepsf t4_only.psf
writepdb t4_only.pdb

# crystal water
pdbalias residue HOH TIP3
segment XTAL {
 auto none
 pdb pdb/256l_water.pdb
}

pdbalias atom HOH O OH2
coordpdb pdb/256l_water.pdb XTAL

guesscoord

# write
writepsf t4_xtal.psf
writepdb t4_xtal.pdb

# solvation
mol delete all
package require solvate
solvate t4_xtal.psf t4_xtal.pdb -rotate -t 12.0 -o t4_solv

# add counter ions
mol delete all
package require autoionize
autoionize -psf t4_solv.psf -pdb t4_solv.pdb -neutralize -cation SOD -anion CLA -seg ION -o t4_ion

# check boxsize
mol delete all
mol new t4_ion.psf
mol addfile t4_ion.pdb
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
mol new t4_ion.psf
mol addfile t4_ion.pdb
set all [atomselect top "all"]
set fixed 1.0
set nofixed 0
set fixed_list {}
foreach i [$all get protein] {
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
mol new t4_ion.psf
mol addfile t4_ion.pdb
set all [atomselect top "all"]
set kforce 10.0
set kforce_list {}
foreach i [$all get protein] {
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
readpsf t4_ion.psf
coordpdb t4_ion.pdb

writepsf x-plor t4_xplor.psf
writepsf charmm t4_charmm.psf

