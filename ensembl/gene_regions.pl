# a script to dump all genic regions from ensembl to a BED format file
# use the merge_gene_regions to collapse the transcript level information
# down to overlapping regions

use strict;
use warnings;

use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
    -host => 'ens-livemirror',
    -user => 'ensro'
);

my $ga = $reg->get_adaptor('human', 'core', 'gene');

#my $genes = $ga->fetch_all_by_external_name('ZAR1L');
my $genes = $ga->fetch_all;

while (my $g = shift @$genes) {

    my $c = 'chr'.$g->seq_region_name;

    for my $t (@{ $g->get_all_Transcripts }) {

        my $id = $t->stable_id;

        if ($t->translation) {
            my $tm = $t->get_TranscriptMapper;

            for my $e ($tm->cdna2genomic($t->cdna_coding_start, $t->cdna_coding_start + 2)) {
                if ($e->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                    print join("\t", $c, $e->start - 1, $e->end, $id, 'START'), "\n";
                }
            }
            
            for my $e ($tm->cdna2genomic($t->cdna_coding_end - 2, $t->cdna_coding_end)) {
                if ($e->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                    print join("\t", $c, $e->start - 1, $e->end, $id, 'STOP'), "\n";
                }
            }
        }

        for my $e (@{ $t->get_all_Exons }) {
            print join("\t", $c, $e->seq_region_start - 1, $e->seq_region_end, $id, 'EXON'), "\n";

            if ($t->seq_region_strand == 1) {
                
                # forward strand transcript

                if (my $cds_start = $e->coding_region_start($t)) {

                    if (my $cds_end = $e->coding_region_end($t)) {
                        print join("\t", $c, $cds_start - 1, $cds_end, $id, 'CDS'), "\n";
                        
                        if ($cds_end < $e->seq_region_end) {
                            print join("\t", $c, $cds_end, $e->seq_region_end, $id, 'UTR3'), "\n";
                        }
                    }
                    
                    if ($cds_start > $e->seq_region_start) {
                        print join("\t", $c, $e->seq_region_start - 1, $cds_start - 1, $id, 'UTR5'), "\n";
                    }
                    
                }
                elsif ($t->translation) {
                    # whole exon is UTR
                    print join("\t", $c, $e->seq_region_start - 1, $e->seq_region_end, $id, ($t->coding_region_start > $e->seq_region_start ? 'UTR5' : 'UTR3')), "\n";
                }
            }
            else {
                # transcript on the reverse strand
               
                if (my $cds_start = $e->coding_region_start($t)) {

                    if (my $cds_end = $e->coding_region_end($t)) {
                        print join("\t", $c, $cds_start - 1, $cds_end, $id, 'CDS'), "\n";
                        
                        if ($cds_end < $e->seq_region_end) {
                            print join("\t", $c, $cds_end, $e->seq_region_end, $id, 'UTR5'), "\n";
                        }
                    }
                    
                    if ($cds_start > $e->seq_region_start) {
                        print join("\t", $c, $e->seq_region_start - 1, $cds_start - 1, $id, 'UTR3'), "\n";
                    }
                    
                }
                elsif ($t->translation) {
                    # whole exon is UTR
                    print join("\t", $c, $e->seq_region_start - 1, $e->seq_region_end, $id, ($t->coding_region_start > $e->seq_region_start ? 'UTR3' : 'UTR5')), "\n";
                }
            }
        }
        
        for my $i (@{ $t->get_all_Introns }) {
            my $il = join "\t", $c, $i->seq_region_start - 1, $i->seq_region_end, $id, 'INTRON';
            print "$il\n";
            my $ss1 = join "\t", $c, $i->seq_region_start - 1, $i->seq_region_start + 1, $id, ($i->strand == 1 ? 'DONOR' : 'ACCEPTOR');
            my $ss2 = join "\t", $c, $i->seq_region_end - 2, $i->seq_region_end, $id, ($i->strand == 1 ? 'ACCEPTOR' : 'DONOR');

            print "$ss1\n$ss2\n";
        }        
    }
}

