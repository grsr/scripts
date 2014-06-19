use strict;
use warnings;

my ($in_dir, $out_dir, $perms) = @ARGV;

my $idx = $ENV{'LSB_JOBINDEX'};

my @files = `ls $in_dir`;

my $file = $files[$idx - 1];

chomp $file;

my $name = $file;

$name =~ s/\.[^\.]*$//;

my $cmd = "./GSEA_new $in_dir/$file -p $perms -o $out_dir/$name.out -s $out_dir/$name.txt";

`$cmd`;

