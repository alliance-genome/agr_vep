=head1 NAME

TranscriptName

=head1 SYNOPSIS

 cp TranscriptName.pm ~/.vep/Plugins/
 ./vep -i variations.vcf --plugin TranscriptName,name=db_name,host=db_host,user=db_user,port=db_port,pass=db_pass

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

    $self->{db_name} = $param_hash->{name};
    $self->{db_host} = $param_hash->{host};
    $self->{db_user} = $param_hash->{user};
    $self->{db_port} = $param_hash->{port};
    $self->{db_pass} = $param_hash->{pass};

    $self->{query} = "SELECT transcript_name FROM transcript_map WHERE transcript_id = ?";
    
    $self->{dsn} = 'dbi:mysql:database=' . $self->{db_name} . ';host=' . $self->{db_host} . ';port=' . $self->{db_port};
    $self->{dbh} ||= DBI->connect($self->{dsn}, $self->{db_user}, $self->{db_pass}) or die $DBI::errstr;
    $self->{get_sth} = $self->{dbh}->prepare($self->{query});

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
	if ($$ != $self->{initial_pid}) { #forked, reconnect to DB
	    $self->{dbh} = DBI->connect($self->{dsn}, $self->{db_user}, $self->{db_pass});
	    $self->{get_sth} = $self->{dbh}->prepare($self->{query});

	    #set this so only do once per fork
	    $self->{initial_pid} = $$;
	}

	$self->{get_sth}->execute($tr->stable_id);

	while (my $arrayref = $self->{get_sth}->fetchrow_arrayref) {
	    $transcript_name = $arrayref->[0];
	}
	$self->add_to_cache($tr->stable_id, $transcript_name);
    }

    my $results = {transcript_name => $transcript_name};

    return $results;
}


sub fetch_from_cache {
    my ($self, $transcript_id) = @_;

    my $cache = $self->{_transcript_id_name_map_cache} ||= [];

    my ($data) = map {$_->{name}} grep {$_->{id} eq $transcript_id} @$cache;

    return $data;
}


sub add_to_cache {
    my ($self, $transcript_id, $transcript_name) = @_;

    my $cache = $self->{_transcript_id_name_map_cache} ||= [];
    push @$cache, {id => $transcript_id, name => $transcript_name};

    shift @$cache while scalar @$cache > 50;

    return;
}


1;
