# merge the transcript level output of gene_regions.pl to a BED file

use strict;
use warnings;

my ($in , $out) = @ARGV;

open my $OUT, ">$out";

my @regions = qw(EXON INTRON CDS START STOP DONOR ACCEPTOR UTR5 UTR3);

for my $reg (@regions) {
    `grep $reg $in > $reg.bed`;
    `mergeBed -nms -i $reg.bed > ${reg}_merged.bed`;
    open my $MF, "<${reg}_merged.bed";
    while (<$MF>) {
        chomp; 
        print $OUT "$_\t$reg\n"; 
    }
    close $MF;
    `rm $reg.bed`;
    `rm ${reg}_merged.bed`;
}

close $OUT;

`sort -k1,1 -k2,2n $out -o $out`;
