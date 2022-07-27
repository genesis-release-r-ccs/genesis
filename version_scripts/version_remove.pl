#!/usr/bin/perl
#
use strict;
use warnings;
while(<>) {
	if (/GENESIS_VERSION/) {
		$_=~ s/\$GENESIS_VERSION\$: \S* \S* \S* \S*/\$GENESIS_VERSION\$/;
	} else {
		$_=~s/$_//;
	}
} continue {
	print $_;
}
