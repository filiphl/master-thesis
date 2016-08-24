#type <logfile console> in extensions->tk console to display
#gui input as commands
#-------#
# Input #
#-------#
#read from file
set filename "./vmd.dat";
set fp [open $filename r];
set filedata [read $fp];
close $fp;
set data [split $filedata "\n"];

set path [lindex $data 0]; #"/home/$USER/simulations/H2ORun";
set infile [lindex $data 1]; #"sorted.xyz";
set last [lindex $data 2]; #1000
set step [lindex $data 3]; #1

#-------------------#
# derived variables #
#-------------------#
set last [expr $last + 1];
set length [expr $last - 1]
mol new [format "%s/00000/%s" $path $infile];
puts "\033\[1A\033\[K";
#------#
# Main #
#------#


puts "\n\n";
for {set i $step} {$i < $last} {set i [expr $i + $step]} {
	# puts "\033\[0m\033\[1A\033\[K\033\[1A\033\[K\033\[1A\033\[K\033\[1A\033\[KReading file $i of $length\033\[8m"
	mol addfile [format "%s/%05d/%s" $path $i $infile];
}
# puts"\033\[0m\033\[1A\033\[K\033\[1A\033\[K\033\[1A\033\[K\033\[1A\033\[K";

#Delete line drawing
mol delrep 0 0


# Changing graphic representations of atoms
mol selection atomicnumber 0;
mol color colorid 4;
mol material Diffuse;
mol representation VDW 1.0 12;
# mol material opaque;
mol addrep top;

mol selection atomicnumber 1;
mol color colorid 1;
mol material Diffuse;
mol representation VDW 0.7 12;
# mol material opaque;
mol addrep top;

# Changing graphic representations of atoms
mol selection atomicnumber 2;
mol color colorid 0;
mol material Diffuse;
mol representation VDW 0.4 12;
# mol material opaque;
mol addrep top;

mol selection atomicnumber 3;
mol color colorid 1;
mol material Diffuse;
mol representation VDW 0.7 12;
# mol material opaque;
mol addrep top;

set sel [atomselect top "atomicnumber 0"];
$sel set radius 1;
set sel [atomselect top "atomicnumber 1"];
$sel set radius 1;
set sel [atomselect top "atomicnumber 2"];
$sel set radius 1;
set sel [atomselect top "atomicnumber 3"];
$sel set radius 1;

unset last;
unset step;
unset path;
unset infile;
unset i;
unset sel;