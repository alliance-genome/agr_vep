=head1 NAME

GenomePos

=head1 SYNOPSIS

 mv GenomePos.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin GenomePos

=head1

 A VEP plugin that adds the genomic position of variations.  Insertion start and end positions
 are given as the bases either side of the insertion point.

=cut

package GenomePos;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);


sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);

    $self->{initial_pid} = $$;

    return $self;
}


sub get_header_info {
    my $self = shift;

    return {
	Genomic_start_position => 'Genomic start position of variation on chromosome/contig',
	Genomic_end_position => 'Genomic end position of variation on chromosome/contig'
    };
}


sub feature_types {
    return ['Feature', 'Intergenic'];
}


sub variant_feature_types {
    return ['BaseVariationFeature'];
}


sub run {
    my ($self, $bvfoa) = @_;

    my $bvf = $bvfoa->base_variation_feature;

    my $results = {};
    $results->{Genomic_start_position} = $bvf->seq_region_start;
    $results->{Genomic_end_position} = $bvf->seq_region_end;

    return $results;
}

1;
