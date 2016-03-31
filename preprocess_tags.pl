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

my $version = "process_tags.pl v. 0.2\n";
my $help    = "Use as:
  process_tags.pl [arguments] <reads.fq(.gz)> <specimens.txt>

  Input files:

       reads.fq(.gz) : the fastq file with the sequences (optionally gzipped)
       specimens.txt : a list of specimen names with associated P5 and P7 index sequences

  Arguments:
       -m : Matching length                                                       [default: 10]
       -i : Index read length                                                     [8]
       -x : allowed mismatches when comparing mDBRs for PCR-duplication detection [1]
       -y : allowed mismatches when matching P7 index                             [1]
       -t : Trim this number of bases from the beginning of every sequence        [2]
       -Z : don't gZip output files (flag)
       -d : Dry run (no output) (flag)
       -o : Output folder                                                         [.]
       -f : use only this Fraction of the reads                                   [1]
";

our($opt_m, $opt_i, $opt_x, $opt_y, $opt_t, $opt_o, $opt_f, $opt_Z, $opt_d, $opt_h, $opt_v);
my $list   = pop;
my $fqfile = pop;
die $help    if !defined $list || !defined $fqfile || !getopts('m:i:x:y:t:o:f:Zdhv') || $opt_h;
die $version if $opt_v;

sub VERSION_MESSAGE {
	my ($fp) = @_;
	print { $fp } $version;
}

sub HELP_MESSAGE {
	my ($fp) = @_;
	print { $fp } $help;
}

-f $list   || die "'$list' file not found\n";
-f $fqfile || die "'$fqfile' file not found\n";

## defaults

$opt_m    = 10  if !defined $opt_m;
$opt_x    = 1   if !defined $opt_x;
$opt_y    = 1   if !defined $opt_y;
$opt_t    = 2   if !defined $opt_t;
$opt_i    = 8   if !defined $opt_i;
$opt_f    = 1   if !defined $opt_f;
$opt_o    = "." if !defined $opt_o;

## check input values

die "-m must be a postive integer\n" if $opt_m =~ /\D/;
die "-x must be a postive integer\n" if $opt_x =~ /\D/;
die "-y must be a postive integer\n" if $opt_y =~ /\D/;
die "-t must be a postive integer\n" if $opt_t =~ /\D/;
die "-i must be a postive integer\n" if $opt_i =~ /\D/;
die "-f must be a float number (0;1]\n" if $opt_f <= 0 || $opt_f > 1;
die "Failed to create the output folder: $!\n" if ! -d $opt_o && !mkdir($opt_o);

my $wet_run = !defined $opt_d;
my $m2 = $opt_m * 2;
my $m3 = $opt_m * 3;
my (%loci, %combs);

my $start_time = time();

my $fqp = openin($fqfile);
my (@specimens, @i7s, @i5s);
my %i5s_rescue;
my @fps;
my (@loci, @uniq);
my (@pieces12, @pieces13, @pieces23);
my %i5s_rescue_bas;
my %i5s_rescue_pos;
my @atgc = ("A","T","G","C");

open(LIST, "<", $list) || die "$!\n";
my $s = 0;

# read info about the specimens
while (<LIST>) {
	chomp;
	next if m/^\s*$/;
	my @line = split;
	next if !$line[0] || $line[0] =~ /^\s*#/;
	die "Wrong number of fields on line $. of $list\n" if $#line < 2;
	my $i7    = uc($line[1]);
	my $i5raw = uc($line[2]);
	die "i7 index '$i7' has incorrect format\n"    if $i7    !~ /^[ATGC]+$/;
	die "i5 index '$i5raw' has incorrect format\n" if $i5raw !~ /^[ACGTMRWSYKVHDBN]+$/;
	die "i7 index must be exactly $opt_i bases\n"  if length($i7)    != $opt_i;
	die "i5 index must be exactly $opt_i bases\n"  if length($i5raw) != $opt_i;
	my $i5 = prepare_dbr($i5raw);
	for (my $j = 0; $j < length($line[2]); $j++) {
		next if substr($i5raw, $j, 1) eq "N";
		my $rescue = prepare_dbr(substr($i5raw, 0, $j)."N".substr($i5raw, $j + 1));
		$i5s_rescue{$i7}{$rescue}   = $i5;
		$i5s_rescue_pos{$rescue}    = $j;
		@{$i5s_rescue_bas{$rescue}} = degenerate_base(substr($i5raw, $j, 1));
	}
	die "Index combination $i7/$i5 met two times\n" if exists $combs{$i7}{$i5};
	die "Specimen '$line[0]' met two times\n" if grep($_ eq $line[0], @specimens);
	push(@i7s, $i7);
	push(@i5s, $i5);
	push(@specimens, $line[0]);
	$combs{$i7}{$i5} = $s++;
}
close LIST;
die "No info could be found in the specimens file\n" if !@specimens;

# open output files if needed
if ($wet_run) {
	my $basename = basename($fqfile);
	$basename =~ s/(\.f(ast)?q)?(\.gz)?$//;
	push(@fps, openout("$basename.noi7.fq"));
	push(@fps, openout("$basename.noi5.fq"));
	foreach my $specimen (@specimens) {
		push(@fps, openout("$specimen.fq"));
	}
}

my $c_total = 0;
my @c_duplic;
my @c_all;
$c_all[0] = 0;
$c_all[1] = 0;

my (@fq, @i7, @i5);
my ($f, $i5r);
my ($i5, $i7, $i5m, $i7m);
my $i;
my $mask;
my $i5pos;
my @i5bas;
my @header;
my $i5start = -1-$opt_i;
my $f2;

# iterate over the fastq records
while (1) {
	$fq[$_] = <$fqp> foreach (0..3);
	last if !defined($fq[3]);
	next if $opt_f < 1 && rand > $opt_f;
	report_time() if ++$c_total % 1000000 == 0;
	@header = split(/ |:/, $fq[0]);
	die "Indices couldn't be identified on line ".($.-3).":\n$fq[0]\n" if $#header < 10;
	$i7 = substr($header[10],        0, $opt_i);
	$i5 = substr($header[10], $i5start, $opt_i);
	($f, $i5r) = detect_indices($i7, $i5);
	$f2 = $f + 2;
	$c_all[$f2]++;
	if (is_duplicate($f, $fq[1], $i5r)) {
		$c_duplic[$f2]++;
		next;
	}
	if ($opt_t) {
		chomp($fq[1], $fq[3]);
		$fq[1] = substr($fq[1], $opt_t)."\n";
		$fq[3] = substr($fq[3], $opt_t)."\n";
	}
	print { $fps[$f2] } @fq if $wet_run;
}

foreach my $fp (@fps) {
	close $fp if defined $fp;
}

print STDERR "\nJob finished in ";

my $time = time() - $start_time;
if ($time < 3600) { printf STDERR "%d sec\n", $time;      }
else              { printf STDERR "%d min\n", $time / 60; }

# print out the summaries
printf STDERR "Total number of input reads: %d\nAmong them: %d i7 orphans, %d i5 orphans\n",
	$c_total, $c_all[0], $c_all[1];

print STDERR "Duplicate/all reads per library:\n";

for $f (0..$#specimens) {
	$f2 = $f + 2;
	printf STDERR "%s: %d/%d\n", $specimens[$f], defined($c_duplic[$f2]) ? $c_duplic[$f2] : 0, defined($c_all[$f2]) ? $c_all[$f2] : 0;
}

## BEGIN SUBS ##

sub report_time {
	$time = time() - $start_time;
	printf STDERR "Read #%d, %d sec, %d reads/sec\n", $c_total, $time, $c_total / $time if $time > 0;
}

# convert IUPAC-coded mDBR string into a regex pattern
sub prepare_dbr {
	my $pattern = shift;
	$pattern =~ s/R/[AG]/g;
	$pattern =~ s/Y/[CT]/g;
	$pattern =~ s/S/[GC]/g;
	$pattern =~ s/W/[AT]/g;
	$pattern =~ s/K/[GT]/g;
	$pattern =~ s/M/[AC]/g;
	$pattern =~ s/B/[CGT]/g;
	$pattern =~ s/D/[AGT]/g;
	$pattern =~ s/H/[ACT]/g;
	$pattern =~ s/V/[ACG]/g;
	$pattern =~ s/N/[ACGTNX]/g;
	return $pattern;
}

# return a list of non-degenerate bases corresponding to
# the input IUPAC code
sub degenerate_base {
	my $b = shift;
	return ("A")         if $b eq "A";
	return ("T")         if $b eq "T";
	return ("G")         if $b eq "G";
	return ("C")         if $b eq "C";
	return ("A","G")     if $b eq "R";
	return ("C","T")     if $b eq "Y";
	return ("G","C")     if $b eq "S";
	return ("A","T")     if $b eq "W";
	return ("G","T")     if $b eq "K";
	return ("A","C")     if $b eq "M";
	return ("C","G","T") if $b eq "B";
	return ("A","G","T") if $b eq "D";
	return ("A","C","T") if $b eq "H";
	return ("A","C","G") if $b eq "V";
	return ("A","T","G","C");
}

my ($len, $mylen);

# recognized indices from the expected set
sub detect_indices {
	$i7  = shift;
	$i5  = shift;
	$i7m = 0;
	$len = $opt_y;
	foreach $i (@i7s) {
		$mask = $i7 ^ $i;
		$mask =~ tr/\0//d;
		$mylen = length($mask);
		if ($mylen <= $len) {
			$i7m = $i;
			last if $mylen == 0;
			$len = $mylen;
		}
	}
	return (-2, $i5) if !$i7m;
	foreach $i5m (keys %{$combs{$i7m}}) {
		return ($combs{$i7m}{$i5m}, $i5) if $i5 =~ /^($i5m)/;
	}
	my $i;
	while (my ($i, $i5m) = each %{$i5s_rescue{$i7m}}) {
		if ($i5 =~ /^($i)/) {
			$i5pos =   $i5s_rescue_pos{$i};
			@i5bas = @{$i5s_rescue_bas{$i}};
			$i5    = substr($i5, 0, $i5pos).$i5bas[rand @i5bas].substr($i5, $i5pos + 1);
			die("$i7m, $i5m") if !defined($combs{$i7m}{$i5m});
			return ($combs{$i7m}{$i5m}, $i5);
		}
	}
	return (-1, $i5);
}

my $seq;
my ($subseq0, $subseq1, $subseq2);
my ($b0, $b1, $b2, $h0, $h1);
my $pos;

# detect if this is a duplicate
sub is_duplicate {
	my $f   = shift;
	return 0 if $f < 0;
	my $seq = shift;
	my $i5  = shift;
	my $subseq1 = substr($seq, 0,      $opt_m);
	my $subseq2 = substr($seq, $opt_m, $opt_m);
	my $frag = $pieces12[$f]{$subseq1}{$subseq2};
	if (!defined($frag)) {
		my $subseq3 = substr($seq, $m2, $opt_m);
		$frag = $pieces13[$f]{$subseq1}{$subseq3};
		if (!defined($frag)) {
			$frag = $pieces23[$f]{$subseq2}{$subseq3};
			if (!defined($frag)) {
				$frag = push(@{$loci[$f]}, $subseq1.$subseq2.$subseq3) - 1;
				$pieces23[$f]{$subseq2}{$subseq3} = $frag;
			}
			$pieces13[$f]{$subseq1}{$subseq3} = $frag;
		}
		$pieces12[$f]{$subseq1}{$subseq2} = $frag;
	}
	if ($uniq[$f][$frag]{$i5}++) {
		return 1;
	}
	elsif (!$opt_x) {
		return 0;
	}
	foreach my $dbr (keys %{$uniq[$f][$frag]}) {
		next if $i5 eq $dbr;
		my $mask = $i5 ^ $dbr;
		$mask =~ tr/\0//d;
		if (length($mask) <= $opt_x) {
			++$uniq[$f][$frag]{$dbr};
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

# open (optionally gzipped) output file
sub openout {
	my $fname = shift;
	$fname    = "$opt_o/$fname" if $opt_o ne '.';
	my $pipe = defined $opt_Z ? "> $fname" : "| gzip -c > $fname.gz";
	open my $fp, $pipe or die "$!\n";
	return $fp;
}
