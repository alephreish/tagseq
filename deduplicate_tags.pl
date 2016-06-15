#!/usr/bin/perl

### Written by Andrey Rozenberg (jaera at yandex.com), Ruhr-Universität Bochum

### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
### GNU General Public License for more details.

### You should have received a copy of the GNU General Public License
### along with this program. If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
                              
use File::Basename;
use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my $version = "deduplicate_tags.pl v. 0.2\n";
my $help    = "Use as:
  deduplicate_tags.pl [arguments] <reads.fq(.gz)>

  Input file:

       reads.fq(.gz) : the fastq file with the sequences (optionally gzipped)

  Arguments:
       -m : Matching length                                                       [default: 10]
       -i : Index read length                                                     [8]
       -x : allowed mismatches when comparing mDBRs for PCR-duplication detection [1]
       -t : Trim this number of bases from the beginning of every sequence        [2]
       -o : Output file                                                           [- (STDOUT)]
       -f : use only this Fraction of the reads                                   [1]
       -d : Dry run (no output) (switch)
";

our($opt_m, $opt_i, $opt_x, $opt_t, $opt_o, $opt_f, $opt_d, $opt_h, $opt_v);
my $fqfile = pop;
die $help    if !defined $fqfile || !getopts('m:i:x:t:o:f:dhv') || $opt_h;
die $version if $opt_v;

sub VERSION_MESSAGE {
	my ($fp) = @_;
	print { $fp } $version;
}

sub HELP_MESSAGE {
	my ($fp) = @_;
	print { $fp } $help;
}

-f $fqfile || die "'$fqfile' file not found\n";

## defaults

$opt_m    = 10  if !defined $opt_m;
$opt_x    = 1   if !defined $opt_x;
$opt_t    = 2   if !defined $opt_t;
$opt_i    = 8   if !defined $opt_i;
$opt_f    = 1   if !defined $opt_f;
$opt_o    = "-" if !defined $opt_o;

## check input values

die "-m must be a postive integer\n" if $opt_m =~ /\D/;
die "-x must be a postive integer\n" if $opt_x =~ /\D/;
die "-t must be a postive integer\n" if $opt_t =~ /\D/;
die "-i must be a postive integer\n" if $opt_i =~ /\D/;
die "-f must be a float number (0;1]\n" if $opt_f <= 0 || $opt_f > 1;

my $wet_run = !defined $opt_d;
my $m2 = $opt_m * 2;

my $start_time = time();

my $fqp = openin($fqfile);
my (@loci, @uniq);
my (%pieces12, %pieces13, %pieces23);

# open output files if needed
if ($wet_run) {
	my $pipe = ($opt_o =~ /\.gz$/ ? "| gzip > $opt_o" : "> $opt_o");
	open FP, $pipe or die "$!\n";
}

my $c_total = 0;
my $c_duplic = 0;

my @fq;
my $i5;
my @header;
my $i5start = -1-$opt_i;

# iterate over the fastq records
while (1) {
	$fq[$_] = <$fqp> foreach (0..3);
	last if !defined($fq[3]);
	next if $opt_f < 1 && rand > $opt_f;
	report_time() if ++$c_total % 1000000 == 0;
	@header = split(/ |:/, $fq[0]);
	die "Indices couldn't be identified on line ".($.-3).":\n$fq[0]\n" if $#header < 10;
	$i5 = substr($header[10], $i5start, $opt_i);
	if (is_duplicate($fq[1], $i5)) {
		$c_duplic++;
		next;
	}
	if ($opt_t) {
		chomp($fq[1], $fq[3]);
		$fq[1] = substr($fq[1], $opt_t)."\n";
		$fq[3] = substr($fq[3], $opt_t)."\n";
	}
	print FP @fq if $wet_run;
}

close FP if $wet_run;

print STDERR "\nJob finished in ";

my $time = time() - $start_time;
if ($time < 3600) { printf STDERR "%d sec\n", $time;      }
else              { printf STDERR "%d min\n", $time / 60; }

# print out the summaries

printf STDERR "All/duplicate reads: %d/%d\n", $c_total, $c_duplic;

## BEGIN SUBS ##

sub report_time {
	$time = time() - $start_time;
	printf STDERR "Read #%d, %d sec, %d reads/sec\n", $c_total, $time, $c_total / $time if $time > 0;
}

# detect if this is a duplicate
sub is_duplicate {
	my $seq = shift;
	my $i5  = shift;
	my $subseq1 = substr($seq, 0,      $opt_m);
	my $subseq2 = substr($seq, $opt_m, $opt_m);
	my $frag = $pieces12{$subseq1}{$subseq2};
	if (!defined($frag)) {
		my $subseq3 = substr($seq, $m2, $opt_m);
		$frag = $pieces13{$subseq1}{$subseq3};
		if (!defined($frag)) {
			$frag = $pieces23{$subseq2}{$subseq3};
			if (!defined($frag)) {
				$frag = push(@loci, $subseq1.$subseq2.$subseq3) - 1;
				$pieces23{$subseq2}{$subseq3} = $frag;
			}
			$pieces13{$subseq1}{$subseq3} = $frag;
		}
		$pieces12{$subseq1}{$subseq2} = $frag;
	}
	if ($uniq[$frag]{$i5}++) {
		return 1;
	}
	elsif (!$opt_x) {
		return 0;
	}
	foreach my $dbr (keys %{$uniq[$frag]}) {
		next if $i5 eq $dbr;
		my $mask = $i5 ^ $dbr;
		$mask =~ tr/\0//d;
		if (length($mask) <= $opt_x) {
			++$uniq[$frag]{$dbr};
			return 1;
		}
	}
	return 0;
}

# open (optionally gzipped) input file
sub openin {
	my $fname = shift;
	my $pipe = ($fname =~ /\.gz$/ ? "gzip -cd $fname |" : "< $fname");
	open my $fp, $pipe or die "$!\n";
	return $fp;
}
