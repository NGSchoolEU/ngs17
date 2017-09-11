#! /usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Temp qw/tempdir/;

my $USAGE = "
$0 source.fa target.fa remap.list 

Remap the sites from the file remap.list where the sites are
referenced to source.fa with coordinates starting from 1 to target.fa
based on muscle alignment.
";

my $sourcef = shift or die $USAGE;
my $targetf = shift or die $USAGE;
my $remapf = shift or die $USAGE;

my $source=read_single_fa($sourcef);
my $target=read_single_fa($targetf);

my $dir = tempdir("remap_XXXXXXX", TMPDIR => 0, CLEANUP => 1); 
print STDERR "Creating temporary directory $dir\n";

open TOALN,">$dir/toaln.fa";
print TOALN ">source\n$source\n>target\n$target\n";
close TOALN;
my_run("muscle -in $dir/toaln.fa -out $dir/muscle.fa 2>/dev/null");
my $alignment = read_fa("$dir/muscle.fa");

print STDERR $alignment->{'source'},"\n",$alignment->{'target'},"\n";

my %translate;
my $srcc = 0;
my $tgtc = 0;
my $len = length $alignment->{'source'};
for (my $i=0; $i < $len; $i++) {
    my $cs = substr($alignment->{'source'},$i,1);
    my $ct = substr($alignment->{'target'},$i,1);
    $srcc++ unless $cs eq '-';
    $tgtc++ unless $ct eq '-';
    $translate{$srcc}=$tgtc unless ($cs eq '-' || $ct eq '-');
}

open IN,"<$remapf" or die "Cannot open $remapf";
while (my $line = <IN>) {
    chomp $line;
    my @parts = split " ",$line;
    if (defined $translate{$parts[0]}) {
	$parts[0] = $translate{$parts[0]};
	print join(" ",@parts),"\n";
    } else {
	warn("Cannot translate position ".$parts[0].", no alignment");
    }
}
close IN;
exit 0;

sub read_single_fa {
    my ($fname) = @_;
    my $result = read_fa($fname);
    my @names = keys %{$result};
    die "$fname should contain a single sequence" 
	unless (scalar(@names)==1);
    return $result->{$names[0]};
}

sub read_fa {
    my ($fname) = @_;
    my %results;

    open my $in, "<$fname" or die "Cannot open $fname";
    my $seqname;
    my $seq;
   
    while (my $line = <$in>) {
	chomp $line;
	if ($line =~ /^>(.*)$/) {
	    my $newseqname = $1;
	    $results{$seqname} = $seq if (defined $seqname);
	    $seq = "";
	    $seqname = $newseqname;
	} else {
	    $seq .= $line;
	}
    }
    $results{$seqname} = $seq if (defined $seqname);
    close $in;
    return \%results;
}

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
