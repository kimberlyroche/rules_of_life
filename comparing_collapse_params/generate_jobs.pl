use strict;
use warnings;

my @filters = qw(0.2 0.5 0.9);
#my @agglomerations = qw(genus species);
my @agglomerations = qw(species);

for my $a (@agglomerations) {
  for my $f (@filters) {
    open(my $fh, '>', "temp.slurm");
    print $fh '#!/bin/bash'."\n";
    print $fh '#SBATCH -J VC_'.substr($a, 0, 1).'_'.$f."\n";
    print $fh '#SBATCH --mem=32GB'."\n";
    print $fh '#SBATCH --time=24:00:00'."\n";
    print $fh '#'."\n\n";
    print $fh 'module add R/3.4.2-fasrc01'."\n";
    print $fh 'module add gcc/5.3.0-fasrc01'."\n\n";
    print $fh 'cd /data/mukherjeelab/rulesoflife'."\n\n";
    print $fh 'srun Rscript compare_collapsing.R '.$a.' '.$f."\n\n";

    my $call_str = "sbatch temp.slurm";
    print("Calling: ".$call_str."\n");
    my $sys_response = `$call_str`;
    sleep(1);
    `rm temp.slurm`;
  }
}
