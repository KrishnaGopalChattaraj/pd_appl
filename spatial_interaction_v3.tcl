#######################################################################################################################
## Revision: vDW radii are fixed
#######################################################################################################################

#######################################################################################################################
## A procedure which returns the index of CA atom. The argument is the index of any atom of the residue
proc CAindex {idl} {
    set selca [list]
    foreach id $idl {
	set selt [atomselect top "index $id"]
	set selres [$selt get resid]
	lappend selca [[atomselect top "resid $selres and (index > $id - 30 ) and (index < $id + 30) and name CA"] get index]
    }
    return $selca
}
#######################################################################################################################

## Set the effective hydrophobicity cut-off
puts "Enter the hydrophobicity cutoff value (use 0.3 for high accuracy and 0.1 for high coverage)"
gets stdin eff_hydro

puts "Do you want to retain cluster of only one hydrophobic reside (0=no, 1=yes)"
gets stdin cluster_one

## Set the radius R over which neighbors are identified
set R 10.0
set Rcharge 5.0

puts "Performing the SIM calculations..."

## Set the output file names and path
set my_path [molinfo top get filename]
set pdbfiletarg [file rootname $my_path]_sim_$eff_hydro$cluster_one.pdb
set hydrofiletarget [file rootname $my_path]_hydrophobicity_$eff_hydro$cluster_one.log
set simfiletarget [file rootname $my_path]_sim_$eff_hydro$cluster_one.log
set hydrotarget [open $hydrofiletarget w]
set simtarget [open $simfiletarget w]

puts $simtarget "RESNAME\tRESID\tSIM\tHYDRO\tCHAIN"

## ASSIGN Hydrophobicity FOR 20 AMINO ACIDS
source hydrophobicity_scale.tcl

## Assign Lennard-Jones radii parameters from charmm force-field
source vdw_radii.tcl

# READ REFERENCE SAA FOR 20 AMINO ACIDS FROM TRIPEPTIDE VALUES
source saa_fully_exposed.tcl

# Select all the atoms
set allsel [atomselect top all]
set allCA [atomselect top "protein and name CA"]

set protsel [atomselect top protein]

# For static, we set iframe = 0, for dynamics, we iterate over all frames
set iframe 0

# Initialize arrays to hold averaged sim value
$protsel frame 0
$protsel set beta 0
set sap_tot [$protsel get beta]
    
$allsel frame $iframe
$allCA frame $iframe

set alllist [$allsel get index]
set CAlist [$allCA get index]   

$allsel set user 0.0
$allsel set beta 0.0
$allsel set occupancy 0

# set the User field with the effective hydrophobicity value for the selected residue 
    foreach id $CAlist { 
	set sel1 [atomselect top "index $id" frame $iframe]
	# Homodimers: Different chain should have either different chain id or different segname
	set sel1_resid [$sel1 get resid]
	set sel1_resname [$sel1 get resname]
	set sel1_chain [$sel1 get chain]
	set sel1_segname [$sel1 get segname]
       	set sel [atomselect top "resid $sel1_resid and chain $sel1_chain and segname $sel1_segname" frame $iframe] 
	set sel_sc [atomselect top "resid $sel1_resid and chain $sel1_chain and segname $sel1_segname and sidechain" frame $iframe] 

	# Measure SASA of the side-chain atoms residue
	set resasa [measure sasa 1.4 $allsel -restrict $sel_sc -samples 50000]
	# -samples 50000 is used for high accuracy. Default value is only 5000. If program becomes very slow, decrease this value.
	# Whenever this value is changed, rerun the saa_tripeptide.tcl with the new value.

	# Mulitply the sasa of each atom by "(residue hydrophobicity - GLY hydrophobicity)/residue fully exposed sasa"
	set sapNORMfactor 0   
	set sap 0
        switch [$sel1 get resname] {
	    "HSD" {set sapNORMfactor [expr ($hphobHSD - $hphobGLY)/$REFHSD]}
	    "HSE" {set sapNORMfactor [expr ($hphobHSE - $hphobGLY)/$REFHSE]}
	    "HSP" {set sapNORMfactor [expr ($hphobHSD - $hphobGLY)/$REFHSP]}
	    "ARG" {set sapNORMfactor [expr ($hphobARG - $hphobGLY)/$REFARG]}
	    "LYS" {set sapNORMfactor [expr ($hphobLYS - $hphobGLY)/$REFLYS]}
	    "ILE" {set sapNORMfactor [expr ($hphobILE - $hphobGLY)/$REFILE]}
	    "PHE" {set sapNORMfactor [expr ($hphobPHE - $hphobGLY)/$REFPHE]}
	    "LEU" {set sapNORMfactor [expr ($hphobLEU - $hphobGLY)/$REFLEU]}
	    "TRP" {set sapNORMfactor [expr ($hphobTRP - $hphobGLY)/$REFTRP]}
	    "ALA" {set sapNORMfactor [expr ($hphobALA - $hphobGLY)/$REFALA]}
	    "MET" {set sapNORMfactor [expr ($hphobMET - $hphobGLY)/$REFMET]}
	    "PRO" {set sapNORMfactor [expr ($hphobPRO - $hphobGLY)/$REFPRO]}
	    "CYS" {set sapNORMfactor [expr ($hphobCYS - $hphobGLY)/$REFCYS]}
	    "ASN" {set sapNORMfactor [expr ($hphobASN - $hphobGLY)/$REFASN]}
	    "VAL" {set sapNORMfactor [expr ($hphobVAL - $hphobGLY)/$REFVAL]}
	    "GLY" {set sapNORMfactor [expr ($hphobGLY - $hphobGLY)/$REFGLY]}
	    "SER" {set sapNORMfactor [expr ($hphobSER - $hphobGLY)/$REFSER]}
	    "GLN" {set sapNORMfactor [expr ($hphobGLN - $hphobGLY)/$REFGLN]}
	    "TYR" {set sapNORMfactor [expr ($hphobTYR - $hphobGLY)/$REFTYR]}
	    "ASP" {set sapNORMfactor [expr ($hphobASP - $hphobGLY)/$REFASP]}
	    "GLU" {set sapNORMfactor [expr ($hphobGLU - $hphobGLY)/$REFGLU]}
	    "THR" {set sapNORMfactor [expr ($hphobTHR - $hphobGLY)/$REFTHR]}
	}

	set sap [expr {$resasa*$sapNORMfactor}]
	$sel set user $sap 
	$sel set user2 $resasa
	puts $hydrotarget [format "%d %s %f %f" $sel1_resid $sel1_resname $sap $resasa]
	$sel delete 
	$sel_sc delete
	$sel1 delete
    }
close $hydrotarget


## Perform the SIM calculation
foreach id $CAlist { 
    set sel1 [atomselect top "index $id" frame $iframe]
    if {[$sel1 get user] < $eff_hydro} { continue }
    # Homodimers: Different chain should have either different chain id or different segname
    set sel1_resname [$sel1 get resname]
    set sel1_resid [$sel1 get resid]
    set sel1_chain [$sel1 get chain]
    set sel1_segname [$sel1 get segname]
    set sel [atomselect top "resid $sel1_resid and chain $sel1_chain and segname $sel1_segname" frame $iframe]
    set seln [atomselect top "(within $R of (resid $sel1_resid and chain $sel1_chain and segname $sel1_segname)) and (user > $eff_hydro)" frame $iframe]
    set selnlist [$seln get index]

    # For each of the selected atom, we find the CA of its residue
    set CAnlist [CAindex $selnlist]
    # Remove duplicated so that each CA occurs only once in selnlist
    set CAnlistu [lsort -unique $CAnlist]

    set sim 0
    set num_neigh 0
    foreach idn $CAnlistu {
	set selt [atomselect top "index $idn" frame $iframe]
	set simt [$selt get user]
	set sim [expr {$sim + $simt}]
	incr num_neigh
	$selt delete
    }

    if {$cluster_one == 0} {
    if {$num_neigh == 1 } { set sim 0 }
    }
    $sel set beta $sim
    
    ## Find the charged residues in vicinity of this residue. We select only those with SAA > 10 A2
    if {$sim > 0} {
	set selc [atomselect top "(within $Rcharge of (resid $sel1_resid and chain $sel1_chain and segname $sel1_segname)) and (resname LYS or resname ARG or resname ASP or resname GLU) and user2 > 10" frame $iframe]
	set selc_list [$selc get index]
	set selc_ca [CAindex $selc_list]
	foreach idp $selc_ca {
	    set sel_ct [atomselect top "index $idp"]
	    set sel_ct_resid [$sel_ct get resid]
	    set sel_ct_chain [$sel_ct get chain]
	    set sel_ct_segname [$sel_ct get segname]
	    set sel_ct1 [atomselect top "resid $sel_ct_resid and chain $sel_ct_chain and segname $sel_ct_segname" frame $iframe]
	    $sel_ct1 set occupancy 1
	    $sel_ct delete
	    $sel_ct1 delete
	}
	$selc delete
	unset selc_list
	unset selc_ca
    }
    
    $sel delete
    $sel1 delete
    $seln delete
    unset CAnlist
    unset CAnlistu
    unset sim
}


## Print the results
foreach id $CAlist { 
    set sel1 [atomselect top "index $id" frame $iframe]
    if {[$sel1 get beta] > 0} {
	puts $simtarget [format "%s\t%d\t%.2f\t%.2f\t%s" [$sel1 get resname] [$sel1 get resid]  [$sel1 get beta]  [$sel1 get user] [$sel1 get chain]]
    }
    if {[$sel1 get occupancy] == 1} {
	puts $simtarget [format "%s\t%d\t%.2f\t%.2f\t%s"  [$sel1 get resname] [$sel1 get resid]  [$sel1 get beta]  [$sel1 get user] [$sel1 get chain]]
    }
}

## Save the SIM to a protein_sim_mapped.pdb
$allsel writepdb $pdbfiletarg


close $simtarget

puts "Program successfully over !!"
puts "SIM values of residues are stored in $pdbfiletarg"
puts "Residue hydrophobicity values are in $hydrofiletarget"
puts "Residue SIM values are in $simfiletarget"



