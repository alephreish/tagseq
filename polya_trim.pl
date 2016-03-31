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
use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my $version = "polya_trim.pl v. 0.2\n";
my $help    = "Use as:
polya_trim.pl [arguments] < input.fq > output.fq
or
gzip -cd input.fq.gz | polya_trim.pl [arguments] | gzip > output.fq.fq

  Arguments:
    -L : minimum read Length before or after trimming [default: 77]
    -A : trim polyA-tails of at least this length     [1]
    -G : trim polyG-tails of at least this length     [1]
    -5 : trim this Number of bases from 5' end        [0]
    -I : discard reads with polyN and 'slipped' indexes (flag)
";

our($opt_L, $opt_A, $opt_G, $opt_5, $opt_I, $opt_h, $opt_v);
die $help    if !getopts('L:A:G:5:Ihv') || $opt_h;
die $version if $opt_v;

sub VERSION_MESSAGE {
	my ($fp) = @_;
	print { $fp } $version;
}

sub HELP_MESSAGE {
	my ($fp) = @_;
	print { $fp } $help;
}

my @fq;
my ($len1, $len2);
my $b;
my %short;
my $ns1     = 0;
my $ns2     = 0;
my $polys1  = 0;
my $polys2  = 0;
my $all     = 0;
my $success = 0;
my $slip    = 0;

my $len_min   = defined $opt_L ? $opt_L : 77;
my $trimpolyA = defined $opt_A ? $opt_A : 1;
my $trimpolyG = defined $opt_G ? $opt_G : 1;
my $trim5     = defined $opt_5 ? $opt_5 : 0;
my $trimI     = defined $opt_I;

my @header;
my ($n1, $n2);
my ($poly1, $poly2);

while (1) {
	$fq[$_] = <STDIN> foreach (0..3);
	last if !defined $fq[3];
	$all++;
	chomp($fq[0], $fq[1]);
	$len1 = length($fq[1]);
	if ($len1 < $len_min) {
		$short{'too short input'}++;
		next;
	}
	if ($trimI) {
		@header = split(/\s+|\:|\+/, $fq[0]);
		$n1 = () = $header[10] =~ /N/g;
		$n2 = () = $header[11] =~ /N/g;
		$ns1 += $n1 > 2;
		$ns2 += $n2 > 1;
		next if $n1 > 2 || $n2 > 1;
		$poly1 = $header[10] =~ /([ATGC])\1{7}/;
		$poly2 = $header[11] =~ /([ATGC])\1{7}/;
		$polys1 += $poly1;
		$polys2 += $poly2;
		next if $poly1 || $poly2;
		if ($header[11] eq "AGATCTCG") {
			$slip++;
			next;
		}
	}
	$b = substr($fq[1], -1, 1);
	if ($trimpolyG) {
		$fq[1] =~ s/G{$trimpolyG,}$//;
	}
	if ($trimpolyA) {
		$fq[1] =~ s/A{$trimpolyA,}$//;
	}
	$len2 = length($fq[1]);
	if ($len1 != $len2) {
		if ($len2 < $len_min) {
			$short{$b}++;
			next;
		}
		$fq[3] = substr($fq[3], 0, $len2)."\n";
	}
	if ($trim5) {
		$fq[1] = substr($fq[1], $trim5);
		$fq[3] = substr($fq[3], $trim5);
	}
	$fq[0] .= "\n";
	$fq[1] .= "\n";
	$success++;
	print @fq;
}

die "No input provided\n" if !$all;

printf STDERR "Total: %d\nSuccess: %d (%s%%)\nDiscarded reads:\n", $all, $success, $success*100/$all;

for my $base (keys %short) {
	printf STDERR "  %s: %d\n", $base, $short{$base};
}

if ($trimI) {
	printf STDERR "  N's left: %d\nN's right: %d\nPolyN left: %d\nPolyN right: %d\nSlippery: %d", $ns1, $ns2, $polys1, $polys2, $slip;
}
