use strict;
use warnings;

# how much overlap is there between the wishlist of enzymes and pathways Jack sent ("Gene_targets.txt") and
# the Piphillin metagenomics ("enzymes.txt" and "pathways.txt"); are the identifiers the same or do we need
# to do some painful matching???

my $print_enzymes_pathways = 1;

my @list_enzymes_piphillin = ();
my @list_pathways_piphillin = ();
my @targets_jack = ();

# read in enzymes.txt -------------------------------------------------------------------------------------

my $token = "";
open(DATA, "original_data/Piphillin_20190222/enzymes.txt");
while(<DATA>) {
	if($_ =~ /^OTU_ID/) {
	} elsif($_ =~ /^"(.*?)"/) {
		$token = $1;
		push @list_enzymes_piphillin, $token;
	} elsif($_ =~ /^(.*?)\s+\d+\s+\d+/) {
		$token = $1;
		chomp($token);
		push @list_enzymes_piphillin, $token;
	}
}
close(DATA);
@list_enzymes_piphillin = sort(@list_enzymes_piphillin);

# read in pathways.txt ------------------------------------------------------------------------------------

open(DATA, "original_data/Piphillin_20190222/pathways.txt");
while(<DATA>) {
	if($_ =~ /^OTU_ID/) {
	} elsif($_ =~ /^(.*?)\s+[\d\.]+\s+[\d\.]+/) {
		$token = $1;
		chomp($token);
		push @list_pathways_piphillin, $token;
	}
}
close(DATA);
@list_pathways_piphillin = sort(@list_pathways_piphillin);

# write out -----------------------------------------------------------------------------------------------

foreach my $e (@list_enzymes_piphillin) {
	print($e."\n");
}
foreach my $p (@list_pathways_piphillin) {
	print($p."\n");
}

# read in Gene_targets.txt --------------------------------------------------------------------------------

=pod
$token = "";
open(DATA, "original_data/Gene_targets.txt");
while(<DATA>) {
	if($_ =~ /^Genes\tFunctional_categories\t\n/) {
	} elsif($_ =~ /^(.*?)\t/) {
		$token = $1;
		chomp($token);
		push @targets_jack, $token;
	}
}
close(DATA);
@targets_jack = sort(@targets_jack);

my $i = 0;
while($i < 5) {
	print("Enzyme: ".$list_enzymes_piphillin[$i]."\n");
	print("Pathway: ".$list_pathways_piphillin[$i]."\n");
	print("Gene list: ".$targets_jack[$i]."\n\n");
	$i++;
}
=cut
