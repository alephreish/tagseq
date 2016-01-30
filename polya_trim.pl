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

if (defined($ARGV[0]) && $ARGV[0] eq "-h") {
	die "polya_trim.pl [min_len] [trim polyA] [trim polyG] [trim N bases from 5'] [trim Indices] < [input] > [output]
";
}

my $len_min   = 77;
   $len_min   = $ARGV[0] if defined($ARGV[0]);
my $trimpolyG = 1;
   $trimpolyG = $ARGV[1] if defined($ARGV[1]);
my $trimpolyA = 1;
   $trimpolyA = $ARGV[2] if defined($ARGV[2]);
my $trim5     = 0;
   $trim5     = $ARGV[3] if defined($ARGV[3]);
my $trimi     = 1;
   $trimi     = $ARGV[4] if defined($ARGV[4]);

my @header;
my ($n1, $n2);
my ($poly1, $poly2);

while (!eof(STDIN)) {
	$all++;
	$fq[$_] = <STDIN> foreach (0..3);
	chomp($fq[0], $fq[1]);
	$len1 = length($fq[1]);
	if ($len1 < $len_min) {
		$short{'start'}++;
		next;
	}
	if ($trimi) {
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

printf STDERR "Total: %d\nSuccess: %d (%s%%)\npolyN's:\n", $all, $success, $success*100/$all;

for my $base (keys %short) {
	printf STDERR "%s: %d\n", $base, $short{$base};
}

if ($trimi) {
	printf STDERR "N's left: %d\nN's right: %d\nPolyN left: %d\nPolyN right: %d\nSlippery: %d", $ns1, $ns2, $polys1, $polys2, $slip;
}
