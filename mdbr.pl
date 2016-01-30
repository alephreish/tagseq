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

use warnings;
use strict;
use IPC::Open3;

my @dbrs;

if ($#ARGV != 2) {
	die "Use as mdbr.pl <DBR length> <difference threshold> <alphabet>\n";
}

my $len  = shift @ARGV;
my $thre = shift @ARGV;


if ($len eq 'own') {
	@dbrs = uniq(@ARGV);
}
else {
	my $abjad = uc shift @ARGV;
	die "The alphabet $abjad is not supported: use only RY or ATGC symbols\n" if $abjad !~ /^[RY]+$/ && $abjad !~ /^[ATGC]+$/;
	my @ab = split //, $abjad;
	recurs("", $len, @ab);
}

`cliquer -h 2> /dev/null` or die "cliquer is needed to run this script\n";
open3 \*CLIQ_IN, \*CLIQ_OUT, \*CLIQ_ERR, "cliquer - -aq" or die "Failed to open cliquer: $!\n";

my $n = $#dbrs + 1;
printf CLIQ_IN "p\tmygraph\t%d\t0\n", $n;
for (my $i = 0; $i < $n; $i++) {
	my $dbr1 = $dbrs[$i];
	for (my $j = $i + 1; $j < $n; $j++) {
		my $dbr2 = $dbrs[$j];
		my $mask = mymask($dbr1, $dbr2);
		my $mlen = masklen($mask);
		printf CLIQ_IN "e\t%d\t%d\n", $i + 1, $j + 1 if $mlen >= $thre;
	}
}

close CLIQ_IN;
while (<CLIQ_ERR>) { print STDERR if !/^Reading/ }
my @parts;
my @nodes;
{
	local $, = ' ';
	while (<CLIQ_OUT>) {
		chomp;
		@parts = split /:\s*/;
		@nodes = ();
		my $first = 1;
		foreach my $i (split /\s+/, $parts[1]) {
			next if !$i;
			die "Unidentified node: $i\n" if !defined($dbrs[$i-1]);
			push @nodes, $dbrs[$i-1];
		}
		print @nodes;
		print "\n";
	}
}
printf STDERR "Clique size: %d; number of cliques: %d\n", $#nodes + 1, $.;
close CLIQ_OUT;
close CLIQ_ERR;


sub masklen {
	my $mask = shift;
	$mask =~ tr/\0//d;
	return length($mask);
}

sub mymask {
	my $a = shift;
	my $b = shift;
	my $mask = $a ^ $b;
	return $mask;
}

sub uniq {
	my %seen;
	grep !$seen{$_}++, @_;
}

sub recurs {
	my $mer = shift;
	my $len = shift;
	my $ret = length($mer) == $len - 1;
	foreach my $a (@_) {
		my $dbr = $mer.$a;
		if ($ret) {
			push(@dbrs, $dbr);
		}
		else {
			recurs($dbr, $len, @_);
		}
	}
}
