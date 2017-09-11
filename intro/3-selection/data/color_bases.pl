#! /usr/bin/perl -w

use strict;

my $USAGE = "$0 pdb_file positions_list threshold_red threshold_yellow output_prefix";
die $USAGE unless @ARGV==5;
my ($pdb_filename, $list_filename, $th1, $th2, $output_basename) = @ARGV;
die $USAGE unless -r $pdb_filename && -r $list_filename && $th2<=$th1;

my $out;
open $out, ">", "$output_basename.pml" or die;
print $out "load $pdb_filename\n";
print $out "hide all
show cartoon
#bg_color white
color green
";


my $in;
open $in, "<", $list_filename or die;
while(my $line = <$in>) {
    my @parts = split ' ', $line;
    die unless @parts==2;
    if($parts[1]>=$th1) { 
	print $out "color red, (resi $parts[0])\n";
	print $out "show spheres, $parts[0]/ca\n";
    }
    elsif($parts[1]>=$th2) {
	print $out "color yellow, (resi $parts[0])\n";
	print $out "show spheres, $parts[0]/ca\n";
    }
}
close $in or die;

print $out "zoom complete=1
save $output_basename.pse\n";
close $out or die;

my_run("pymol $output_basename.pml -c");

############################
sub my_run
{
    my ($run, $die) = @_;
    if(!defined($die)) { $die = 1; }

    my $short = substr($run, 0, 20);

    print STDERR $run, "\n";
    my $res = system("bash", "-c", $run);
    if($res<0) {
        die "Error in program '$short...' '$!'";
    }
    if($? && $die) {
        my $exit  = $? >> 8;
        my $signal  = $? & 127;
        my $core = $? & 128;

        die "Error in program '$short...' "
            . "(exit: $exit, signal: $signal, dumped: $core)\n\n ";
    }
}


