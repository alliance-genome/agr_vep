=head1 NAME

ProtSeq

=head1 SYNOPSIS

 mv ProtSeq.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin ProtSeq

=head1

 A VEP plugin that adds wild-type/variant protein sequences

=cut

package ProtSeq;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);


sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);
   
    return $self;
}


sub get_header_info {
    my $self = shift;

    return {
	WtSeq => 'Wild-type translated amino acid sequence',
	VarSeq => 'Variant translated amino acid sequence'
    };
}


sub feature_types {
    return ['Transcript'];
}


sub variant_feature_types {
    return ['VariationFeature'];
}


sub run {
    my ($self, $tva) = @_;

    my $results = {};

    my $tr = $tva->transcript;
    my $tv = $tva->transcript_variation;

    my $wt_translation = $tr->translate;
    return $results unless $wt_translation;
    $results->{'WtSeq'} = $wt_translation->seq;

    return $results unless $tva->peptide;

    my $tl_start = $tv->translation_start;
    my $tl_end = $tv->translation_end;

    my $var_translation = $results->{'WtSeq'};
    if ($tva->peptide =~ /X$/) {
	substr($var_translation, $tl_start - 1) = $tva->peptide;
    }
    else {
	substr($var_translation, $tl_start - 1, $tl_end - $tl_start + 1) = $tva->peptide;
    }
    $results->{'VarSeq'} = $var_translation;

    return $results;
}

1;
