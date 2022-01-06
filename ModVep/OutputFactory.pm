package ModVep::OutputFactory;

use strict;
use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');


sub run {
    my $self = shift;

    my $working_dir = $self->required_param('vep_working');
    my $chr_dir_nr = 1;

    my @chromosome_nrs;
    while (-d "${working_dir}/${chr_dir_nr}") {
	push @chromosome_nrs, $chr_dir_nr;
	$chr_dir_nr++;
    }

    $self->param('chromosome_nrs', [map {{chromosome_nr => $_}} @chromosome_nrs]);  
}

sub write_output {
    my $self = shift;

    $self->dataflow_output_id($self->param('chromosome_nrs'), 1);
}

1;
