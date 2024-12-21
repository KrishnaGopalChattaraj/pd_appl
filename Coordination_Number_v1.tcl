# Load the PSF and DCD files
mol new mAb_ion.psf
mol addfile mAb_min_all_5ns.dcd type dcd first 0 step 1 waitfor all

# Define the output file
set output "average_calcium_beta_H_I_L_M.pdb"

# Get the total number of frames
set num_frames [molinfo top get numframes]
puts "Number of frames: $num_frames"


#################################
# Calculation for segment H
#################################
set segment_H "segname H"
set residues_H [atomselect top "$segment_H and name CA"]
set residue_ids_H [$residues_H get resid]

array unset cal_sums_H
foreach r $residue_ids_H {
    set cal_sums_H($r) 0
}

# Loop over all frames for segment H
for {set i 0} {$i < $num_frames} {incr i} {
    animate goto $i
    foreach r $residue_ids_H {
        set ca_atom [atomselect top "$segment_H and resid $r and name CA"]
        if {[$ca_atom num] > 0} {
            # Count CAL within 6 Å of this residue in segment H
            set nearby_cal [atomselect top "resname CAL and within 6 of ($segment_H and resid $r)"]
            set count [llength [$nearby_cal get index]]
            set cal_sums_H($r) [expr {$cal_sums_H($r) + $count}]
            $nearby_cal delete
        }
        $ca_atom delete
    }
}

array unset cal_avgs_H
foreach r $residue_ids_H {
    set cal_avgs_H($r) [expr {double($cal_sums_H($r)) / $num_frames}]
}

$residues_H delete


#################################
# Calculation for segment I
#################################
set segment_I "segname I"
set residues_I [atomselect top "$segment_I and name CA"]
set residue_ids_I [$residues_I get resid]

array unset cal_sums_I
foreach r $residue_ids_I {
    set cal_sums_I($r) 0
}

# Loop over all frames for segment I
for {set i 0} {$i < $num_frames} {incr i} {
    animate goto $i
    foreach r $residue_ids_I {
        set ca_atom [atomselect top "$segment_I and resid $r and name CA"]
        if {[$ca_atom num] > 0} {
            # Count CAL within 6 Å of this residue in segment I
            set nearby_cal [atomselect top "resname CAL and within 6 of ($segment_I and resid $r)"]
            set count [llength [$nearby_cal get index]]
            set cal_sums_I($r) [expr {$cal_sums_I($r) + $count}]
            $nearby_cal delete
        }
        $ca_atom delete
    }
}

array unset cal_avgs_I
foreach r $residue_ids_I {
    set cal_avgs_I($r) [expr {double($cal_sums_I($r)) / $num_frames}]
}

$residues_I delete


#################################
# Calculation for segment L
#################################
set segment_L "segname L"
set residues_L [atomselect top "$segment_L and name CA"]
set residue_ids_L [$residues_L get resid]

array unset cal_sums_L
foreach r $residue_ids_L {
    set cal_sums_L($r) 0
}

# Loop over all frames for segment L
for {set i 0} {$i < $num_frames} {incr i} {
    animate goto $i
    foreach r $residue_ids_L {
        set ca_atom [atomselect top "$segment_L and resid $r and name CA"]
        if {[$ca_atom num] > 0} {
            # Count CAL within 6 Å of this residue in segment I
            set nearby_cal [atomselect top "resname CAL and within 6 of ($segment_L and resid $r)"]
            set count [llength [$nearby_cal get index]]
            set cal_sums_L($r) [expr {$cal_sums_L($r) + $count}]
            $nearby_cal delete
        }
        $ca_atom delete
    }
}

array unset cal_avgs_L
foreach r $residue_ids_L {
    set cal_avgs_L($r) [expr {double($cal_sums_L($r)) / $num_frames}]
}

$residues_L delete


#################################
# Calculation for segment M
#################################
set segment_M "segname M"
set residues_M [atomselect top "$segment_M and name CA"]
set residue_ids_M [$residues_M get resid]

array unset cal_sums_M
foreach r $residue_ids_M {
    set cal_sums_M($r) 0
}

# Loop over all frames for segment M
for {set i 0} {$i < $num_frames} {incr i} {
    animate goto $i
    foreach r $residue_ids_M {
        set ca_atom [atomselect top "$segment_M and resid $r and name CA"]
        if {[$ca_atom num] > 0} {
            # Count CAL within 6 Å of this residue in segment I
            set nearby_cal [atomselect top "resname CAL and within 6 of ($segment_M and resid $r)"]
            set count [llength [$nearby_cal get index]]
            set cal_sums_M($r) [expr {$cal_sums_M($r) + $count}]
            $nearby_cal delete
        }
        $ca_atom delete
    }
}

array unset cal_avgs_M
foreach r $residue_ids_M {
    set cal_avgs_M($r) [expr {double($cal_sums_M($r)) / $num_frames}]
}

$residues_M delete



#################################
# Assign Beta Values and Write Output
#################################
# Go to the last frame (or any chosen frame) before writing
animate goto [expr {$num_frames - 1}]

# Update CA atoms in segment H
foreach r $residue_ids_H {
    set ca_atom [atomselect top "$segment_H and resid $r"]
    if {[$ca_atom num] > 0} {
        $ca_atom set beta $cal_avgs_H($r)
    }
    $ca_atom delete
}

# Update CA atoms in segment I
foreach r $residue_ids_I {
    set ca_atom [atomselect top "$segment_I and resid $r"]
    if {[$ca_atom num] > 0} {
        $ca_atom set beta $cal_avgs_I($r)
    }
    $ca_atom delete
}

# Update CA atoms in segment L
foreach r $residue_ids_L {
    set ca_atom [atomselect top "$segment_L and resid $r"]
    if {[$ca_atom num] > 0} {
        $ca_atom set beta $cal_avgs_L($r)
    }
    $ca_atom delete
}

# Update CA atoms in segment M
foreach r $residue_ids_M {
    set ca_atom [atomselect top "$segment_M and resid $r"]
    if {[$ca_atom num] > 0} {
        $ca_atom set beta $cal_avgs_M($r)
    }
    $ca_atom delete
}



# Write out all atoms with updated beta values for CA atoms in H and I
set all_atoms [atomselect top "all"]
$all_atoms writepdb $output
puts "Average CAL coordination values for segments H and I written to CA atoms in $output"

$all_atoms delete

