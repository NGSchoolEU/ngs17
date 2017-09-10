#! /usr/bin/perl -w

use strict;
use Bio::Tools::CodonTable;
use File::Temp qw/tempdir/;
use Data::Dumper;
use File::Copy;

my $PAML = "/usr/bin/codeml";

my $USAGE = "USAGE:
  $0 alignment.ph tree.nh reference_species threshold output_prefix

Output files:
output_prefix.fa     protein for reference species
output_prefix.list   list of significant positions (based 1) wrt protein above 
output_prefix.paml   paml output
output_prefix.log    standard output from paml run 
output_prefix.codons list of significant codons, human readable

Paml executable: $PAML
";
die $USAGE unless @ARGV==5;
my ($aln_filename, $tree_filename, $ref_species, $threshold, $output_basename) = @ARGV;
die $USAGE unless -r $aln_filename && -r $tree_filename && -x $PAML;

my $dir = tempdir("paml_XXXXXXX", TMPDIR => 0, CLEANUP => 0); 
print STDERR "Creating temporary directory $dir\n";

copy($aln_filename, "$dir/aln.ph") or die "Copy failed: $!";
copy($tree_filename, "$dir/tree.nh") or die "Copy failed: $!";
my $out;
open $out, ">", "$dir/codeml.ctl" or die;
print_control_file($out);
close $out;

my_run("cd $dir; $PAML codeml.ctl > paml.log");
die unless -r "$dir/paml.out";
copy("$dir/paml.out", $output_basename . ".paml") or die "Copy failed: $!";
copy("$dir/paml.log", $output_basename . ".log") or die "Copy failed: $!";

my $in;
open $in, "<", "$dir/aln.ph" or die;
my $seqs = read_ph($in);
close $in or die;

open $in, "<", "$dir/paml.out" or die;
my $sites = read_paml($in, $threshold);
close $in or die;

#print list of codons and amino acids at significant positions
open $out, ">", "$output_basename.codons" or die;
foreach my $site (@$sites) {
    # compute nucleotide position, assuming aa position numbered from 1
    my $pos = $site->{'ntpos'};
    print $out $site->{'aa'};
    for(my $j=0; $j<@$seqs; $j++) {
	my $codon = substr($seqs->[$j]{'seq'}, $pos, 3);
	printf $out " %s:(%s,%s)", $seqs->[$j]{'name'}, $codon, translate($codon, 0);
    }
    printf $out " %d %s\n", $site->{'pos'}, $site->{'prob'};
}
close $out;

#find reference sequence
my $ref;
foreach my $seq (@$seqs) {
    if($seq->{'name'} eq $ref_species) {
	die if defined $ref;
	$ref = $seq;
    }
}
die unless defined $ref;

#print protein translation of reference sequence (without gaps)
my $ref_seq = $ref->{'seq'};
$ref_seq =~ s/-//g;
my $prot_seq = translate($ref_seq);
open $out, ">", "$output_basename.fa" or die;
write_fasta($out, ">" . $ref->{'name'}, \$prot_seq);
close $out or die;

#print protein positions wrt protein written to fasta file
open $out, ">", "$output_basename.list" or die;
foreach my $site (@$sites) {
    my $pos = $site->{'ntpos'};  
    my $str = substr($ref->{'seq'}, 0, $pos);
    $str =~ s/-//g;
    my $len = length($str);
    die unless $len%3==0;
    printf $out "%d %s\n", $len/3+1, $site->{'prob'};
}
close $out;

exit 0;

##########
sub read_paml {
    my ($in, $threshold) = @_;
    my @result;

    while(my $line = <$in>) {
	if($line=~ /^Bayes Empirical Bayes \(BEB\) analysis/) {
	    last;
	}
    }
    my $line = <$in>;
    die unless defined $line;
    die unless $line=~/^Positive sites for foreground lineages/;
    while(my $line = <$in>) {
	my @parts = split ' ', $line;
	if(@parts==0) { last; }
	die unless @parts==3 && $parts[0]=~/^\d+$/ && $parts[1]=~/^[A-Z*]$/;
	$parts[2] =~ s/\**$//;
	die unless $parts[2]=~/^[0-9.]+$/;
	my $prob = $parts[2];
	if($prob>=$threshold) {
	    # compute nucleotide position within alignment based-0
            # if input is aa-position within alignment based-1
	    my $pos = ($parts[0]-1)*3;  
	    push @result, {'pos'=>$parts[0], 
			   'ntpos'=>$pos, 'aa'=>$parts[1], 'prob'=>$prob};
	}
    }
    return \@result;
}

##############
sub read_ph {
    my ($in) = @_;
    my @result;
    my $line = <$in>;
    die unless defined $line;
    my @parts = split ' ', $line;
    die unless @parts==2;
    my ($num_species, $len) = @parts;
    for(my $sp = 0; $sp<$num_species; $sp++) {
	my $name = <$in>;
	die unless defined $name;
	$name=~ s/\s*$//;
	my $seq = "";
	while(length($seq)<$len) {
	    my $line = <$in>;
	    die unless defined $line;
	    $line =~ s/\s//g;
	    $seq .= $line;
	}
	die unless length($seq)==$len;
	push @result, {'name'=>$name, 'seq'=>$seq};
    }
    $line = <$in>;
    die if defined $line;
    return \@result;
}

############################
sub write_fasta {
    my ($file, $name, $seq) = @_;

    print $file $name, "\n";
    my $n = length($$seq);
    my $linelen = 60;
    my $i=0;
    while($i<$n) {
        print $file (substr($$seq, $i, $linelen), "\n");
        $i+=$linelen;
    }
}

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

#################
sub translate
{
    my ($seq, $code) = @_;
    my $CodonTable = Bio::Tools::CodonTable->new( -id => $code);
    my $result = $CodonTable->translate($seq);

    return $result;
}

################
sub print_control_file
{
    my ($out) = @_;

    print $out '
      seqfile = aln.ph        * sequence data filename
     treefile = tree.nh       * tree structure file name
      outfile = paml.out      * main result file name

        noisy = 0  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
   aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

        model = 2
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 2  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0
                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
                   * AA: 0:rates, 1:separate

    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
        omega = .4 * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 3  * # of categories in dG of NSsites models

        getSE = 0  * 0: dont want them, 1: want S.E.s of estimates
 RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6
    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
*  fix_blength = 1  * 0: ignore, -1: random, 1: initial, 2: fixed
        method = 0   * 0: simultaneous; 1: one branch at a time


* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
';

}
