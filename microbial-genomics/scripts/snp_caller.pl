#! /usr/local/bin/perl -w

use strict;
use warnings;

use lib 'scripts';

use POSIX;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

###############################################################################
# Process command line options
#

my ( $qual_cutoff, $dp4_cutoff, $dp_cutoff, $dpmax_cutoff, $af_cutoff, $mq_cutoff, @vcfs, $chrom, 
		$ignorefile, $noindels, $phylip, $noheader, $showfiltered, $refilter, $uncertain, $cpus, $behaviour,
		$gt, $site, $outdir, $verbose, $help );

GetOptions(
           'qual|q=i'     => \$qual_cutoff,
           'dp4=i'        => \$dp4_cutoff,
           'dp=i'         => \$dp_cutoff,
           'dpmax=i'      => \$dpmax_cutoff,
           'site=i'       => \$site,
           'af=f'         => \$af_cutoff,
           'mq=f'         => \$mq_cutoff,
           'vcf=s'        => \@vcfs,
           'chrom=s'      => \$chrom,
           'gt'           => \$gt,
           'ignore=s'     => \$ignorefile,
           'showfiltered' => \$showfiltered,
           'noindels'     => \$noindels,
           'phylip=s'     => \$phylip,
           'dir=s'        => \$outdir,
           'noheader'     => \$noheader,
           'refilter'     => \$refilter,
           'uncertain'    => \$uncertain,
           'cpus=i'       => \$cpus,
           'behaviour|b=i' => \$behaviour,
           'verbose=i'    => \$verbose,
           'help|?'       => \$help,
          ) or pod2usage(2);

# Defaults
 $qual_cutoff  ||= 30;
 $dp_cutoff    ||=  0;
 $dpmax_cutoff ||=  100000;
 $dp4_cutoff   ||=  0;
 $af_cutoff    ||=  0;
 $mq_cutoff    ||=  0;
 $showfiltered ||=  0;
 $phylip       ||=  0;
 $outdir       ||=  './';
 $verbose      ||=  0;
 $gt           ||=  0;
 $refilter     ||=  0;
 $cpus         ||=  1;
 $behaviour    ||=  1;
 $uncertain    ||=  0;
 $site         ||=  0;
 
 # Force INDEL filtering if we want phylip format
 $noindels = 1 if $phylip;
 
pod2usage(-exitstatus => 0, -verbose => 2) if $help;

if ( $help || !@vcfs || !$chrom || ( $phylip && !$noindels ) ) { pod2usage(1) }

###############################################################################

my @ignore_temp;
my @ignore;

if ( $ignorefile ) {
	print  "Processing regions to ignore: ".$ignorefile."\n";
	open(IGNORE, $ignorefile) || die "\nCannot open $ignorefile: $!\n";

	while (<IGNORE>) {
		my $line = $_;
		chomp $line;
	
		next if $line =~ m/^#/;
		next if $line eq '';

		my ($start, $end, $name) = split("\t", $line);

#		print $name."\t".$start."\t".$end."\n";
		
		@ignore_temp[$start..$end] = map { 1 } ( $start..$end );
	}	
	close(IGNORE);
}

@ignore = map { defined $_ ? 1 : 0 } @ignore_temp;

## NEED TO CHECK THESE ARE IN ORDER
#print join("\n", @ignore)."\n";
#exit(0);

###############################################################################
## Process VCF files
###############################################################################

# First split labels from filenames in --vcf option

my @vcffiles;

# set up labels hash with reference sequence
my $vcf_labels = { $chrom => $chrom };

foreach my $vcf ( @vcfs ) {
	my ($label, $file ) = split(",", $vcf);
	
	die "\n --vcf label,filename\n is required\n" if ! $label || ! $file;
	die "\n repeat --vcf file : ".$file."\n" if defined $vcf_labels->{$file};
	
	push(@vcffiles, $file);
	$vcf_labels->{$file} = $label;
}

##########################################
## Set up multi core forking
##########################################

my $pm = new Parallel::ForkManager($cpus); 

$pm->run_on_finish(\&fork_run_on_finish);
$pm->run_on_start(\&fork_run_on_start);
#$pm->run_on_wait(\&fork_run_on_wait, 1);

##########################################

my $uncertain_outfile = $outdir.'/'.$phylip.'.uncertain.txt';
open (UNCERTAIN, ">".$uncertain_outfile) || die "\nCannot open $uncertain_outfile: $!\n";

my @filtered_tags = ('AF', 'DP', 'DPMAX', 'DP4', 'MQ', 'QUAL', 'SBREF', 'SBALT', 'SB', 'GT', 'INDEL', 'HET');
print UNCERTAIN "VCF\t".join("\t", @filtered_tags)."\tUncertain\tVariant\tTotal\tSangerHET\n";

our $data;
our @uncertain_sites;
our $non_core;

my $snp_filter = $outdir.'/'.$phylip.'-snp_filter';

unless ( -d $snp_filter ) {
	mkdir $snp_filter;
}
	
foreach my $vcffile ( @vcffiles ) {
#	print  "Reading: ".$vcffile."\n" if $verbose;
	
	my $vcf_outfile = $vcffile;
	$vcf_outfile =~ s/.+\///;
	
	my $vcf_out = $snp_filter.'/'.$vcf_outfile.'.txt';
	open (VCFOUT, ">".$vcf_out) || die "\nCannot open $vcf_out: $!\n";

	$pm->start($vcffile) and next;
	
	my $vcf_data;
	my $vcf_uncertain_sites;
	my $vcf_filtered_sites;
	my $missing_sites;

	if ( $vcffile =~ m/\.gz$/ ) {
		open(VCF, "gunzip -c $vcffile | ") || die "\nCannot open gzipped $vcffile: $!\n";
	} else {
		open(VCF, $vcffile) || die "\nCannot open $vcffile: $!\n";
	}

	my $flag = 0;
	my $uncertain_total = 0;
	my $variant_total = 0;
	my $total = 0;
	
	my $het_counter = 0;
	
	my $position_index = 0;

	while (<VCF>) {
		my $line = $_;
		chomp $line;
	
		next if $line eq '';
	
		if ( $flag == 0 ) {
			if ( $line =~ m/CHROM\tPOS\tID/ ) {
				print $line."\n" unless $noheader;
				$flag++;
				next;
			} else {
				print $line."\n" unless $noheader;
				next;
			}
		}
	
		my @filtered;
		my ($seqname, $pos, $id, $ref, $alt, $qual, $filter, $info, $format1, $format2) = split("\t", $line);
		
		# Check that we haven't skipped any sites, these are gaps in mapping, if so set to '-'
		if ( $pos - $position_index > 1 ) {
		
			for (my $i = ($position_index + 1); $i < $pos; $i++ ) {
	#			$vcf_uncertain_sites->[$i] = 1;
				$missing_sites->{$i}++;

#die "Site ".$i." already present in ".$vcffile."\n" if $vcf_data->{$seqname}->{$i}->{$vcffile};
#
#				if ( ! $vcf_data->{$seqname}->{$i}->{$vcffile} ) {
#					$vcf_data->{$seqname}->{$i}->{$vcffile} = { alt => '-' };
#				}
#
#				if ( ! $vcf_data->{$seqname}->{$i}->{$chrom} ) {
#					$vcf_data->{$seqname}->{$i}->{$chrom} = { ref => '-' };
#				}
#				
				print  "MISSING1: ".$vcffile.": ".$position_index."..".$pos."\t".$i."\n" if $verbose > 2;
				
				print VCFOUT $seqname."\t".$i."\t[MISSING]\n";
			}		

			print  "MISSING2: ".$vcffile.": ".$position_index."..".$pos."\n" if $verbose > 2;
		}
		$position_index = $pos;	
		
		# Only look at sites from this chromosome
		next if $seqname ne $chrom;
		
		# remove anything in the ignore regions
		next if $ignorefile && $ignore[$pos];

		my $type = $info =~ m/INDEL/ ? 'INDEL' : 'SNP';

		push(@filtered, 'INDEL') if $noindels && $type eq 'INDEL'; # have to skip these for behaviour 2
		next if $noindels && $type eq 'INDEL';
		push(@filtered, 'HET') if $alt =~ m/,/;
		
		# QUAL filtering
		push(@filtered, 'QUAL') if $qual < $qual_cutoff;

		# MQ filtering
		if ( $info =~ m/MQ=([\.\d]+);/ ) {
			my $mq = $1;
			
			push(@filtered, 'MQ') if $mq < $mq_cutoff; 
		} else {
			push(@filtered, 'MQ')
		}

		# Site depth of coverage filtering
		if ( $info =~ m/DP=(\d+);/ ) {
			my $dp = $1;
			
			push(@filtered, 'DP') if $dp < $dp_cutoff; 
			push(@filtered, 'DPMAX') if $dp > $dpmax_cutoff; 
		} else {
			push(@filtered, 'DP')
		}
		
		if ( $af_cutoff ) {
			if ( $info =~ m/AF1=([\.\d]+);/ ) {
				my $af = $1;

				push(@filtered, 'AF') if ( $alt eq '.' && $af > 0 ); 
				push(@filtered, 'AF') if ( $alt ne '.' && $af < $af_cutoff ); 
#			} else {
#				push(@filtered, 'AF');
			}
		}

		# Strand bias filtering
		if ( $info =~ m/DP4=(\d+),(\d+),(\d+),(\d+)/ ) {
			my $reff = $1;
			my $refr = $2;
			my $altf = $3;
			my $altr = $4;
		
			my $total = $reff + $refr + $altf + $altr;

			my $strand_min_depth = $dp_cutoff / 2;
		
			if ( $alt eq '.' ) {
			
				# strand bias for all bases
				push(@filtered, 'SBREF') if ( $reff < $strand_min_depth || $refr < $strand_min_depth );

				$het_counter++ if @filtered == 0 && ( $altf + $altr ) >= $dp_cutoff && (( $altf + $altr ) / $total > 0.05) && $altf >= $strand_min_depth && $altr >= $strand_min_depth;

				push(@filtered, 'DP4') if ( ( $reff + $refr ) / $total < ( $dp4_cutoff / 100 ) );

			} else {
			
				# strand bias for variants
				push(@filtered, 'SBALT') if ( $altf < $strand_min_depth || $altr < $strand_min_depth );				

				$het_counter++ if @filtered == 0 && ( $reff + $refr ) >= $dp_cutoff && (( $reff + $refr ) / $total > 0.05) && $reff >= $strand_min_depth && $refr >= $strand_min_depth ;

				push(@filtered, 'DP4') if ( ( $altf + $altr ) / $total < ( $dp4_cutoff / 100 ) ); 
			}
			
		} else {
			push(@filtered, 'SB');
		}

		# Strand bias filtering
		if ( $gt && $format2 && $format2 =~ m/0\/1/ ) {
			push(@filtered, 'GT');
		}
		
		if ( $site ) {
			if ( $site == $pos ) {
				print "\n\n".$line."\n";
				print join(", ", @filtered)."\n\n";
				exit(0);
			} else {
				next;
			}
		}
		
		if ( @filtered ) {
				print $line."\t".join(",", @filtered)."\n" if $showfiltered;
		} else {
				print $line."\n" if ! $showfiltered && ! $phylip;
		}
		
		# If site is filtered do not store it
	
		
# 		if ( @filtered ) {
# 			$uncertain_sites[$pos] = 1;
# 			$uncertain_total++;
# 		} else {
# 			# Collect site information for SNPs
# 			
# 			if ( $alt ne '.' ) {
# 				$data->{$seqname}->{$pos}->{$vcffile} = { ref => $ref, alt => $alt };
# 				$data->{$chrom}->{$pos}->{$chrom} = $ref;
# 			}
# 		}
		
		# Count total number of sites in file
		$total++;
		
		
		if ( $behaviour == 1 ) {           ## Conservative approach, ignore all sites that are filtered in any one strain
		
			if ( @filtered ) {
				$vcf_uncertain_sites->[$pos] = 1;
				$uncertain_total++;
			}
#			} else {			
				if ( $alt ne '.' ) {
				$alt =~ s/,.//;
					$vcf_data->{$seqname}->{$pos}->{$vcffile} = { ref => $ref, alt => $alt };
					$vcf_data->{$chrom}->{$pos}->{$chrom}     = { ref => $ref };
					$variant_total++;
				} else {
					# these will be set to ref base later, too much memory if set here
				}
#			}
		
		} elsif ( $behaviour == 2 ) {      ## if a site is filtered then use an N for that site

			if ( @filtered ) {
				$alt = 'N';
				$uncertain_total++;
			}		
			
			if ( $alt ne '.' ) {
				$vcf_data->{$seqname}->{$pos}->{$vcffile} = { ref => $ref, alt => $alt };
				$vcf_data->{$chrom}->{$pos}->{$chrom}     = { ref => $ref };
				$variant_total++;
			} else {
				# these will be set to ref base later, too much memory if set here
			}
			
		} elsif ( $behaviour == 3 ) {      ## if a site is filtered then 

			if ( @filtered ) {
				$vcf_uncertain_sites->[$pos] = 1;
				$uncertain_total++;
			} else {			
				if ( $alt ne '.' ) {
					$vcf_data->{$seqname}->{$pos}->{$vcffile} = { ref => $ref, alt => $alt };
					$vcf_data->{$chrom}->{$pos}->{$chrom}     = { ref => $ref };
					$variant_total++;
				} else {
					$vcf_data->{$chrom}->{$pos}->{$chrom}     = { ref => $ref }; # watch memory??
				}
			}

		} else {
			die "\nUnknown behaviour\n"
		}
		
		foreach ( @filtered ) {
			$vcf_filtered_sites->{$vcffile}->{$_}++;
		}

		print VCFOUT $line."\t[".join(", ", @filtered)."]\n";
	}

	close(VCFOUT);
		
	my $hets = $vcf_filtered_sites->{$vcffile}->{'HET'} ? $vcf_filtered_sites->{$vcffile}->{'HET'} : 0;
	
	print  "\t".$vcffile."\tUncertain sites = ".$uncertain_total." / ".$total."\tVariant sites = ".$variant_total."\n" if $verbose > 1;
	print  "\t".$vcffile."\tvcf_data = ".join(", ", map { $_." = ".(keys %{$vcf_data->{$_}}) } sort keys %{$vcf_data})."\n" if $verbose > 1;
	
	my $counts = {
		uncertain_total => $uncertain_total,
		variant_total => $variant_total,
		total => $total,
		het_counter => $het_counter
	};
	
	# Generate summary file for filtered data
#	print UNCERTAIN $vcffile;
#	foreach (@filtered_tags) {
#		print UNCERTAIN "\t".($vcf_filtered_sites->{$vcffile}->{$_} ? $vcf_filtered_sites->{$vcffile}->{$_} : 0);
#	}
#	print UNCERTAIN "\t".$uncertain_total."\t".$variant_total."\t".$total."\t".$het_counter."\n";

	$pm->finish(0, { file => $vcffile, missing => $missing_sites, data => $vcf_data, uncertain_sites => $vcf_uncertain_sites, vcf_filtered_sites => $vcf_filtered_sites, counts => $counts, filtered_tags => \@filtered_tags });
}


## Wait for child processes to finish
$pm->wait_all_children;

foreach ( keys %$data ) {
	print "\nCombined total variant sites: ".$_."\t".(keys %{$data->{$_}})." / ".@uncertain_sites."\n" if $verbose;
}

#my $window = 100;
#my $counter = 0;
#my $counter_total = 0;
#for (my $i = 0; $i < @uncertain_sites; $i++ ) {
#	
#	if ( $counter == ($window-1) ) {
#		print "".($i+1)."\t".$counter_total."\n";
#		$counter_total = 0;
#		$counter = 0;
#	} else {
#		$counter_total += $uncertain_sites[$i];
#		$counter++;
#	}
#}

#my @temp;
#my $temp_count = 0;
#for (my $i = 0; $i < @uncertain_sites; $i++ ) {
##	print $i."\t".$uncertain_sites[$i]."\n";
#	$temp[$uncertain_sites[$i]]++;
#
#	if ( $uncertain_sites[$i] > 0 && $uncertain_sites[$i] < 2 ) { #(0.1 * 86) ) {
#		$uncertain_sites[$i] = 0;
#		$temp_count++;
#	}
#}
#print "Rescued sites: ".$temp_count."\n";

#print "\nTEMP: \n";
#for (my $i = 0; $i < @temp; $i++ ) {
#	print $i."\t".$temp[$i]."\n";
#}

print "\nMissing sites: ".(keys %{$non_core})."\n" if $verbose;

###############################################################################
## Process data collected from VCF
###############################################################################

# fill in rest of array positions with a 0
my @uncertain_sites2 = map { $_ ? 1 : 0 } @uncertain_sites;

my $uncertain_sites2_count = 0;
map { $uncertain_sites2_count++ if $_ > 0 } @uncertain_sites2;
print  "Combined uncertain sites: ".$uncertain_sites2_count." / ".@uncertain_sites2."\n" if $verbose;

# if we don't need to generate Phylip output format for clustering then exit here
exit(0) unless $phylip;

my $posfile = $outdir.'/'.$phylip.'.positions.txt';
open(PHYLIP, ">".$outdir.'/'.$phylip) || die "\nCannot open $outdir/$phylip: $!\n";
open(PHYLIPPOS, ">".$posfile) || die "\nCannot open $posfile: $!\n";

my $combined_data;
my $position_data;

my @combined_pos;   # used to highlight SNPs in positions file

foreach my $seq ( sort keys %{$data} ) {
	print  "\nReference sequence: ".$seq."\n" if $verbose;
	
	my $data_counter1 = 0;
	my $data_counter2 = 0;

	foreach my $pos ( sort { $a <=> $b } keys %{$data->{$seq}} ) {
		
		$data_counter1++;
		
		# skip missing sites (non-core)
		next if $non_core->{$pos};

		# skip uncertain sites
		next if $uncertain_sites2[$pos] && !$uncertain;
		
		$data_counter2++;
				
		print  $seq."\t".$pos if $verbose > 2;

		my $ref = $data->{$chrom}->{$pos}->{$chrom}->{ref};
		$position_data->{$pos}->{$chrom} = $ref;
		
		foreach my $vcf ( @vcffiles ) {
		
			my $base_call;
			
			if ( $data->{$seq}->{$pos}->{$vcf} ) {			
				my $alt = $data->{$seq}->{$pos}->{$vcf}->{alt};
				print  "\t".$vcf."\t".$alt if $verbose > 2;
				
				die "\nThis shouldn't happen: ".$pos." :: ".$alt."\n" if $alt =~ m/,.+$/;

				if ( defined $ref && $ref ne $data->{$seq}->{$pos}->{$vcf}->{ref} ) {
					warn "Ref has changed : ".$ref." --> ".$data->{$seq}->{$pos}->{$vcf}->{ref}."\t".$vcf." : ".$pos."\n";
				} else {
					$ref = $data->{$seq}->{$pos}->{$vcf}->{ref};
				}
				
				$base_call = $alt;
			} else { # if not called as variant, uncertain or gap, must be ref 
				print  "\tBASE CALL: ".$seq.','.$pos.','.$vcf."\t".$data->{$seq}->{$pos}->{$vcf}."\n" if $verbose > 2;
				$base_call = $ref;
			}

			push(@{$combined_data->{$vcf}}, $base_call);				
			$position_data->{$pos}->{$vcf} = $base_call;
		}

		push(@{$combined_data->{$chrom}}, $ref);				
		push(@combined_pos, $pos);
		
		print  "\t".$ref."\n" if $verbose > 2;
	}

	print  "\tPre uncertain = ".$data_counter1."\tPost uncertain = ".$data_counter2."\n\n" if $verbose;
}

# Count number of N's in each pseudo sequence
my $counter1;
my $temp_total = @{$combined_data->{$chrom}} || 0;
foreach ( my $i = 0; $i < @{$combined_data->{$chrom}}; $i++ ) {
	foreach my $vcf ( keys %{$combined_data} ) {
		next if $vcf eq $chrom;
	
		my $alt = $combined_data->{$vcf}->[$i];

		if ( $alt eq 'N' ) {
			$counter1->{$vcf}++;
		}		
	}
}

foreach my $vcf ( keys %{$combined_data} ) {
	next if $vcf eq $chrom;
	print "DEBUG: ".$vcf." N count = ".($counter1->{$vcf} ? $counter1->{$vcf} : 0)." / ".$temp_total."\n" if $verbose == 2;
}

## Filter out ref only SNPs if required

if ( $refilter ) {
	print  "Filtering out Ref only SNPs: ".$chrom." :: ".@{$combined_data->{$chrom}}."\n" if $verbose;
	
	my $no_taxa = (keys %{$combined_data});
	print  "Taxa found: ".$no_taxa."\n" if $verbose;

	my $core_sites;
	
	foreach ( my $i = 0; $i < @{$combined_data->{$chrom}}; $i++ ) {
			
		my $ref = $combined_data->{$chrom}->[$i];
		my $alt = '#';
		my $flag = 0;
		my $num_n = 0;
		
		foreach my $vcf ( keys %{$combined_data} ) {
			next if $vcf eq $chrom;
			
			# Check if this is an N, if so don't compare with other sites
			if ( $combined_data->{$vcf}->[$i] eq 'N' ) {
				$num_n++;
				next;
			}
			
			# comparing sites
			
if( !$combined_data->{$vcf}->[$i] ) {
	print "DD: ".$vcf."\t".$i."\n";
}

#print "DD: ".$vcf."\t".$i."\t".$alt."\n";

			$alt = $combined_data->{$vcf}->[$i] if $alt eq '#';

			if ( $alt ne $combined_data->{$vcf}->[$i] ) {
				$flag++;
			}
		}
		
		# Remove site if all non-ref bases are N; -1 as reference is included in %{$combined_data}
		if ( $num_n == ( $no_taxa - 1 ) ) {
			$flag = 0; 

#print "EE: ".$i."\t".$alt."\n";

			# count non-variable sites for BEAST output below
			$core_sites->{$alt}++;
		}
	
		# mark site with 'X' for removal in split below
		if ( $flag == 0 ) {
			foreach my $vcf ( keys %{$combined_data} ) {
				$combined_data->{$vcf}->[$i] = 'X';
			}
		}
	}
	
	print  "Counting invarient core sites (for BEAST)\n" if $verbose; # --".$core_sites->{'A'}."--\n";
	foreach ( keys %{$core_sites} ) {
		print "\t".$_."\t".$core_sites->{$_}."\n";
	}
	
	## Remove the sites marked as 'X' above
	
	my $splice_once = 0;
	foreach my $vcf ( keys %{$combined_data} ) {
	
		# find indicies, but reverse so that we don't have renumbering problems in splice
		# (splice will remove an element and renumber all elements after it, reverse avoids this problem)
	
		my @del_indexes = reverse(grep { $combined_data->{$vcf}->[$_] eq 'X' } 0..$#{$combined_data->{$vcf}});

		foreach my $index (@del_indexes) {
		   	splice (@{$combined_data->{$vcf}}, $index, 1);
		   	splice (@combined_pos, $index, 1) unless $splice_once > 0;
		}
		
		$splice_once++;
	}
}

# Count number of N's in each pseudo sequence
my $counter2;
$temp_total = @{$combined_data->{$chrom}} || 0;
foreach ( my $i = 0; $i < @{$combined_data->{$chrom}}; $i++ ) {
	foreach my $vcf ( keys %{$combined_data} ) {
		next if $vcf eq $chrom;
	
		my $alt = $combined_data->{$vcf}->[$i];

		if ( $alt eq 'N' ) {
			$counter2->{$vcf}++;
		}		
	}
}

foreach my $vcf ( keys %{$combined_data} ) {
	next if $vcf eq $chrom;
	print "DEBUG: ".$vcf." N count = ".($counter2->{$vcf} ? $counter2->{$vcf} : 0)." / ".$temp_total."\n" if $verbose == 2;
}

## Generate PHYLIP file

print PHYLIP " ".(keys %{$combined_data})." ".@{$combined_data->{$chrom}}."\n";
foreach my $vcf ( keys %{$combined_data} ) {
	
	my $name = $vcf_labels->{$vcf};
	
#	printf PHYLIP ("%- 15s", $name);
	printf PHYLIP $name."  ";

	for ( my $i = 0; $i < @{$combined_data->{$vcf}}; $i++) {
		print PHYLIP $combined_data->{$vcf}->[$i]; 		
	}
	
#print  "\n";
	print PHYLIP "\n";
}

## Generate PHYLIP positions file

my @line;
foreach my $vcffile ( sort @vcffiles ) {
	push(@line, $vcf_labels->{$vcffile}); 
}
print PHYLIPPOS "Site\tGap\t".join("\t", @line)."\t".$chrom."\n";

my $previous_pos = 0;

my @kept_sites;
foreach ( @combined_pos ) {
	$kept_sites[$_] = 1;
}

foreach my $pos ( sort { $a <=> $b } keys %{$position_data} ) {
	
	my @line;
	push(@line, $pos);
	push(@line, ($pos - $previous_pos));
	$previous_pos = $pos;
	
	foreach my $vcf ( sort keys %{$position_data->{$pos}} ) {
		next if $vcf eq $chrom;
		push(@line, $position_data->{$pos}->{$vcf});
	}
	
	push(@line, $position_data->{$pos}->{$chrom});
	push(@line, 'SNP') if $kept_sites[$pos];
	
	print PHYLIPPOS join("\t", @line)."\n";
}

###############################################################################
## Subroutines
###############################################################################

##############################
# ForkManager subroutines
##############################

sub fork_run_on_finish {
	my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result) = @_;
	
	if ( $exit_code == 0 && $exit_signal == 0 && $core_dump == 0 ) {
		print "\tSUCCESS ".$ident." [pid: ".$pid."]\n" if $verbose;
	} else {
		print "\tERROR ".$ident." [pid: ".$pid."]";
		print "\t[EXIT_CODE: ".$exit_code."\tEXIT_SIGNAL: ".$exit_signal."\tCORE_DUMP: ".$core_dump."]\n";
	}
	
	foreach my $seqname ( keys %{$result->{data}} ) {
		foreach my $pos ( keys %{$result->{data}->{$seqname}} ) {

			if ( $data->{$seqname}->{$pos} ) {
				# Merge hashes for each pos
				@{$data->{$seqname}->{$pos}}{keys %{$result->{data}->{$seqname}->{$pos}}} = values %{$result->{data}->{$seqname}->{$pos}};
			} else {
				$data->{$seqname}->{$pos} = $result->{data}->{$seqname}->{$pos};
			}
		}
	}
	
	if ( defined $result->{uncertain_sites} ) {
		for ( my $i = 0; $i < @{$result->{uncertain_sites}}; $i++ ) {

			$uncertain_sites[$i] = 0 unless defined $uncertain_sites[$i];
		
			if ( defined $result->{uncertain_sites}->[$i] && $result->{uncertain_sites}->[$i] ) {
				$uncertain_sites[$i]++ if $result->{uncertain_sites}->[$i];
			}
		}
	}

	if ( defined $result->{missing} ) {
		foreach my $site ( keys %{$result->{missing}} ) {
			$non_core->{$site} = 1;

#			die "HHH: ".$site."\t".$non_core->{$site}."\n";
		}
	}
	
	# Generate summary file for filtered data
	print UNCERTAIN $result->{file};
	foreach (@{$result->{filtered_tags}}) {
		print UNCERTAIN "\t".($result->{vcf_filtered_sites}->{$result->{file}}->{$_} ? $result->{vcf_filtered_sites}->{$result->{file}}->{$_} : 0);
	}
	print UNCERTAIN "\t".$result->{counts}->{uncertain_total}."\t".$result->{counts}->{variant_total}."\t".$result->{counts}->{total}."\t".$result->{counts}->{het_counter}."\n";
}

##############################

sub fork_run_on_start {
	my ($pid, $ident) = @_;
	print "Reading ".$ident." [pid: ".$pid."]\n"; 
}

##############################

sub fork_run_on_wait {
	print "** Have to wait for one children ...\n";
}




###############################################################################

__END__

=head1 NAME

filter_vcf.pl - Filters VCF file based on various criteria

=head1 SYNOPSIS

filter_vcf.pl --vcf vcffile --chrom NC_009777

filter_vcf.pl --vcf vcffile1 --vcf vcffile2 --chrom NC_009777 --qual 30 --dp4 75 --dp 10 --dpmax 100 -af 0.75 --showfiltered --noindels --noheader --phylip 
	              --ignorefile Mobiles.txt --verbose 1

=head1 OPTIONS

  --vcf             Input VCF file and label. Format --vcf label,filename (no spaces between , and text)
  --chrom           Reference sequence name to analyse (In case VCF is mapped against multiple sequences)
  --qual            Minimum QUAL score [30]
  --dp4             % reads supoprting SNP [0]
  --dp              Minimum read depth at site [0]
  --dpmax           Maximum read depth at site [100000]
  --af              Allele frequency cutoff e.g. 0.75 [0]
  --showfiltered    Also show all the sites removed (in a Filtered subsection at end),
                       useful to see how well filtering is performing
  --noindels        Filter out INDEL variants
  --phylip          Export data in phylip alignement format (may have to adjust sequence names)
                    (forces --noindels)
  --ignorefile      Specify a text file containing coordinates of regions to ignore e.g. mobile elements
  --noheader        Do not print VCF header
  
  --verbose         print more information (not much at the moment!) 
                    1 - logging
                    2 - debugging
  
=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=head1 AUTHOR

Adam Witney <awitney@sgul.ac.uk>
BuG@S group, Deptartment of Cellular and Molecular Medicine,
St George's, University of London,
London, UK

=head1 COPYRIGHT

bugasbase_pars.pl is Copyright (c) 2009 Adam Witney. UK. All rights reserved.

You may distribute under the terms of either the GNU General Public License or the Artistic License, as specified in the Perl README file.

=cut

=head1 DISCLAIMER

This software comes with no warranty and you use it at your own risk. There may be bugs!

=cut
