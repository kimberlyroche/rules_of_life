use strict;
use warnings;
use POSIX;
use List::MoreUtils qw(uniq);

my $count_table_file = "original_data/Abundances_16S_unrarefied_OTU_table";

# ========================================================================================================
# HOW MANY TAXA MISSING AT EACH LEVEL?
# ========================================================================================================

=pod
my @levels = qw(0 0 0 0 0 0 0);

open(DATA, $count_table_file.".csv");
while(<DATA>) {
	if($_ =~ /^.*,k__(.*?): p__(.*?): c__(.*?): o__(.*?): f__(.*?): g__(.*?): s__(.*?)$/) {
		if($7 eq "") { $levels[6] += 1; }
		if($6 eq "") { $levels[5] += 1; }
		if($5 eq "") { $levels[4] += 1; }
		if($4 eq "") { $levels[3] += 1; }
		if($3 eq "") { $levels[2] += 1; }
		if($2 eq "") { $levels[1] += 1; }
		if($1 eq "") { $levels[0] += 1; }
	}
}
close(DATA);

for(my $i = 0; $i < 7; $i++) {
	print("Missing at ".($i+1).": ".$levels[$i]."\n");
}
=cut

# ========================================================================================================
# HOW MANY UNIQUE GENUS LEVELS?
# ========================================================================================================

=pod
my @tax_strings = [];
my $k = ""; my $p = ""; my $c = ""; my $o = ""; my $f = ""; my $g = ""; my $s = "";
open(DATA, $count_table_file.".csv");
while(<DATA>) {
	if($_ =~ /^(.*),k__(.*?): p__(.*?): c__(.*?): o__(.*?): f__(.*?): g__(.*?): s__(.*?)$/) {
		$k = $2; $p = $3; $c = $4; $o = $5; $f = $6; $g = $7;
	} elsif($_ =~ /^(.*),k__(.*?): p__(.*?): c__(.*?): o__(.*?): f__(.*?): g__(.*?)$/) {
		$k = $2; $p = $3; $c = $4; $o = $5; $f = $6; $g = $7;
	} elsif($_ =~ /^(.*),k__(.*?): p__(.*?): c__(.*?): o__(.*?): f__(.*?)$/) {
		$k = $2; $p = $3; $c = $4; $o = $5; $f = $6; $g = "NA";
	} elsif($_ =~ /^(.*),k__(.*?): p__(.*?): c__(.*?): o__(.*?)$/) {
		$k = $2; $p = $3; $c = $4; $o = $5; $f = "NA"; $g = "NA";
	} elsif($_ =~ /^(.*),k__(.*?): p__(.*?): c__(.*?)$/) {
		$k = $2; $p = $3; $c = $4; $o = "NA"; $f = "NA"; $g = "NA";
	} elsif($_ =~ /^(.*),k__(.*?): p__(.*?)$/) {
		$k = $2; $p = $3; $c = "NA"; $o = "NA"; $f = "NA"; $g = "NA";
	} elsif($_ =~ /^(.*),k__(.*?)$/) {
		$k = $2; $p = "NA"; $c = "NA"; $o = "NA"; $f = "NA"; $g = "NA";
	} elsif($_ =~ /^(.*),Unassigned$/) {
		$k = "NA"; $p = "NA"; $c = "NA"; $o = "NA"; $f = "NA"; $g = "NA";
	}
	if($k eq "") { $k = "NA"; }
	if($p eq "") { $p = "NA"; }
	if($c eq "") { $c = "NA"; }
	if($o eq "") { $o = "NA"; }
	if($f eq "") { $f = "NA"; }
	if($g eq "") { $g = "NA"; }
	push @tax_strings, $k.".".$p.".".$c.".".$o.".".$f.".".$g;
}
close(DATA);

my @uniq_tax = uniq @tax_strings;
print("Unique taxonomic strings (to genus): ".($#uniq_tax + 1)."\n");
=cut

# ========================================================================================================
# AGGLOMERATE TO GENUS
# ========================================================================================================

my %hash = ();
my @tax_strings = [];

open(DATA, $count_table_file.".csv");
my $first_line = "";
my $key = "";
my $missing_count = 0;
my $idx = 0;
my $count_str = "";
my $k = ""; my $p = ""; my $c = ""; my $o = ""; my $f = ""; my $g = ""; my $s = "";
while(<DATA>) {
	if($idx == 0) {
		$first_line = $_;
	} else {
		if($_ =~ /^(.*),k__(.*?): p__(.*?): c__(.*?): o__(.*?): f__(.*?): g__(.*?): s__(.*?)$/) {
			$count_str = $1;
			$k = $2; $p = $3; $c = $4; $o = $5; $f = $6; $g = $7;
		} elsif($_ =~ /^(.*),k__(.*?): p__(.*?): c__(.*?): o__(.*?): f__(.*?): g__(.*?)$/) {
			$count_str = $1;
			$k = $2; $p = $3; $c = $4; $o = $5; $f = $6; $g = $7;
		} elsif($_ =~ /^(.*),k__(.*?): p__(.*?): c__(.*?): o__(.*?): f__(.*?)$/) {
			$count_str = $1;
			$k = $2; $p = $3; $c = $4; $o = $5; $f = $6; $g = "NA";
		} elsif($_ =~ /^(.*),k__(.*?): p__(.*?): c__(.*?): o__(.*?)$/) {
			$count_str = $1;
			$k = $2; $p = $3; $c = $4; $o = $5; $f = "NA"; $g = "NA";
		} elsif($_ =~ /^(.*),k__(.*?): p__(.*?): c__(.*?)$/) {
			$count_str = $1;
			$k = $2; $p = $3; $c = $4; $o = "NA"; $f = "NA"; $g = "NA";
		} elsif($_ =~ /^(.*),k__(.*?): p__(.*?)$/) {
			$count_str = $1;
			$k = $2; $p = $3; $c = "NA"; $o = "NA"; $f = "NA"; $g = "NA";
		} elsif($_ =~ /^(.*),k__(.*?)$/) {
			$count_str = $1;
			$k = $2; $p = "NA"; $c = "NA"; $o = "NA"; $f = "NA"; $g = "NA";
		} elsif($_ =~ /^(.*),Unassigned$/) {
			$count_str = $1;
			$k = "NA"; $p = "NA"; $c = "NA"; $o = "NA"; $f = "NA"; $g = "NA";
		} else {
			print("Missed something with regex!\n");
		}
		if($k eq "") { $k = "NA"; }
		if($p eq "") { $p = "NA"; }
		if($c eq "") { $c = "NA"; }
		if($o eq "") { $o = "NA"; }
		if($f eq "") { $f = "NA"; }
		if($g eq "") { $g = "NA"; }
		$key = $k.".".$p.".".$c.".".$o.".".$f.".".$g;
		push @tax_strings, $key;
		my @temp = split(/,/, $count_str);
		if(exists $hash{$key}) {
			# to push a new index: push(@{$hash{$key}}, $idx);
				for(my $j = 1; $j <= $#temp; $j++) {
					${$hash{$key}}[$j-1] += $temp[$j];
				}
		} else {
			# to push an initial hit index: $hash{$key} = [$idx];
				$hash{$key} = [];
				for(my $j = 1; $j <= $#temp; $j++) {
					push(@{$hash{$key}}, $temp[$j]);
				}
		}
	}
	$idx++;
}
close(DATA);

# my @uniq_tax = uniq @tax_strings;
# print("Unique taxonomic strings (to genus): ".($#uniq_tax + 1)."\n");

=pod
# print for sanity check
my $num_keys = 0;
foreach my $k_it (keys %hash) {
	$num_keys++;
	print($k_it."\n\t");
	for(my $j = 0; $j < $#{$hash{$k_it}}; $j++) {
	 	print(${$hash{$k_it}}[$j]." ");
	}
	print("\n");
}
print($first_line);
print("Keys in hash: ".$num_keys."\n");
=cut

# ========================================================================================================
# WRITE TO FILE, FILTERING TO AT LEAST A 3-COUNT IN 20% OF INDIVIDUALS (9 / 47)
# ========================================================================================================

my $total_taxa = 0;
my $filtered_taxa = 0;
my $original_reads = 0;
my $total_reads = 0;
my $fh_counts = "";
my $k_count = 0;
$k = 3;
my $indiv = floor(0.9*47);
open(my $fh, ">", $count_table_file."_AGG.csv");
print $fh $first_line;
foreach my $k_it (keys %hash) {
	$total_taxa++;
	$k_count = 0;
	$fh_counts = ${$hash{$k_it}}[0];
	my $subtotal_reads = $fh_counts;
	if($fh_counts >= $k) { $k_count++; }
	for(my $j = 1; $j <= $#{$hash{$k_it}}; $j++) {
		my $temp_count = ${$hash{$k_it}}[$j];
		if($temp_count >= $k) { $k_count++; }
	 	$fh_counts = $fh_counts.",".$temp_count;
	 	$subtotal_reads += $temp_count;
	}
	$original_reads += $subtotal_reads;
	if($k_count >= $indiv) {
		$filtered_taxa += 1;
		$total_reads += $subtotal_reads;
		print $fh $fh_counts;
		my @temp = split(/\./, $k_it);
		print $fh ",k__".$temp[0].": ";
		print $fh "p__".$temp[1].": ";
		print $fh "c__".$temp[2].": ";
		print $fh "o__".$temp[3].": ";
		print $fh "f__".$temp[4].": ";
		print $fh "g__".$temp[5]."\n";
	}
}
close $fh;

print("Agglomerated taxa count: ".$total_taxa."\n");
print("Filtered taxa count: ".$filtered_taxa."\n\n");

print("Agglomerated total reads: ".$original_reads."\n");
print("Filtered total reads: ".$total_reads."\n");
print("Percent retained: ".($total_reads/$original_reads)."\n");

