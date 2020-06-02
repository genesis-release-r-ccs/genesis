#!/usr/bin/perl
#
use strict;
use warnings;
my $command='git log --date=short --pretty=format:"%h [%ai]" -1';
my $string=`$command`;
while(<>) {
	$_=~ s/\$GENESIS_VERSION\$/\$GENESIS_VERSION\$: $string/;
} continue {
	print $_;
}
