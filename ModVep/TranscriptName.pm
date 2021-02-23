=head1 NAME

TranscriptName

=head1 SYNOPSIS

 cp TranscriptName.pm ~/.vep/Plugins/
 ./vep -i variations.vcf --plugin TranscriptName,gff=gff_file.gff

=head1

 A VEP plugin that adds the transcript name from the GFF attributes

=cut

package TranscriptName;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);


sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);
    my $param_hash = $self->params_to_hash();

    $self->{gff} = $param_hash->{gff};
    $self->{initial_pid} = $$;

    return $self;
}


sub get_header_info {
    my $self = shift;

    return {
	transcript_name => 'Transcript name obtained from GFF'
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

    my $tr = $tva->transcript;

    my $transcript_name = $self->fetch_from_cache($tr->stable_id);

    if (!$transcript_name) {
	my $id_name_map = $self->map_transcript_name_to_id();
	$self->add_to_cache($id_name_map);
	$transcript_name = $id_name_map->{$tr->stable_id};
    }

    my $results = {};
    $results->{transcript_name} = $transcript_name;

    return $results;
}


sub fetch_from_cache {
    my ($self, $transcript_id) = @_;

    return unless $self->{_transcript_id_name_map_cache};

    return $self->{_transcript_id_name_map_cache}{$transcript_id};
}


sub add_to_cache {
    my ($self, $id_name_map) = @_;

    $self->{_transcript_id_name_map_cache} = $id_name_map;

    return;
}


sub map_transcript_name_to_id {
    my $self = shift;

    my %map;
    my $gff_file = $self->{gff}
    open (GFF, '<', "gunzip -c $gff_file|");
    while (<GFF>) {
	next if $_ =~ /^#/;
	my @columns = split("\t", $_);
	my $attributes = get_attributes($columns[8]);
	my $transcript_id = $attributes->{transcript_id};
	next unless $transcript_id;
	my $transcript_name = $attributes->{Name} || $attributes->{name} || $transcript_id;
	$transcript_name =~ s/\s+$//;
	$map{$transcript_id} = $transcript_name;
    }

    return \%map;
}

sub get_attributes {
    my $info = shift;
    
    my %attributes;
    for my $attr (split(';', $info)) {
	my ($key, $value) = split('=',$attr);
	$attributes{$key} = $value;
    }
    
    return \%attributes;
}

1;
