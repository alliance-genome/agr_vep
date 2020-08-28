package VepProteinFunction::BaseVepProteinFunction;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::AnalysisJob;

use base qw(Bio::EnsEMBL::Hive::Process);

sub get_protein_fasta{
    my ($self, $md5) = @_;

    my $file = $self->required_param('pep_fasta');
    my $fasta = `samtools faidx $file $md5`;

    my $seq = $fasta;
    $seq =~ s/^>.*\n//m;
    $seq =~ s/\s//mg;

    die "No peptide for $md5?\n" unless length $seq > 0;

    return ($fasta, $seq);
}

sub possible_aa_changes {
    my ($self, $codontable_id, $amino_acid) = @_;

    my %possible_aa_changes;

    if ($codontable_id == 1) {
	%possible_aa_changes = (
	    'G' => ['C', 'S', 'R', '*', 'W', 'G', 'V', 'E', 'A', 'D', ],
	    'K' => ['T', 'R', '*', 'K', 'I', 'M', 'E', 'Q', 'N', ],
	    '*' => ['G', 'L', 'W', 'K', '*', 'R', 'C', 'Q', 'Y', 'E', ],
	    'F' => ['F', 'I', 'V', 'Y', 'L', 'C', ],
	    'P' => ['L', 'R', 'T', 'Q', 'H', 'P', 'A', ],
	    'A' => ['A', 'D', 'V', 'E', 'P', 'T', 'G', ],
	    'D' => ['G', 'H', 'N', 'V', 'E', 'A', 'Y', 'D', ],
	    'V' => ['G', 'L', 'A', 'D', 'V', 'M', 'E', 'I', 'F', ],
	    'M' => ['T', 'L', 'K', 'R', 'M', 'V', 'I', ],
	    'L' => ['R', '*', 'W', 'L', 'I', 'V', 'M', 'P', 'H', 'F', 'Q', ],
	    'W' => ['C', 'W', 'L', 'R', '*', 'G', ],
	    'R' => ['C', 'T', 'S', 'L', 'W', 'R', 'I', 'H', 'Q', 'K', '*', 'G', 'M', 'P', ],
	    'T' => ['R', 'K', 'S', 'T', 'P', 'N', 'A', 'I', 'M', ],
	    'S' => ['T', 'S', 'C', 'R', 'G', 'I', 'N', ],
	    'C' => ['F', 'Y', 'W', 'R', '*', 'G', 'C', 'S', ],
	    'Q' => ['E', 'P', 'H', 'Q', 'R', '*', 'L', 'K', ],
	    'H' => ['R', 'L', 'P', 'H', 'N', 'Q', 'D', 'Y', ],
	    'N' => ['H', 'N', 'Y', 'D', 'I', 'K', 'S', 'T', ],
	    'Y' => ['*', 'C', 'N', 'H', 'F', 'D', 'Y', ],
	    'E' => ['G', '*', 'K', 'Q', 'A', 'D', 'V', 'E', ],
	    'I' => ['N', 'F', 'M', 'V', 'I', 'L', 'K', 'R', 'S', 'T', ],
	    'X' => ['X'],
	    );
    }
    elsif ($codontable_id == 2) {
	%possible_aa_changes = (
	    'V' => ['D', 'V', 'A', 'F', 'I', 'L', 'M', 'E', 'G', ],
	    'D' => ['N', 'V', 'D', 'A', 'H', 'Y', 'G', 'E', ],
	    'F' => ['F', 'V', 'Y', 'I', 'L', 'C', ],
	    'P' => ['Q', 'T', 'R', 'H', 'L', 'A', 'P', ],
	    'E' => ['*', 'A', 'V', 'D', 'K', 'Q', 'G', 'E', ],
	    'M' => ['V', 'L', 'I', 'M', 'T', 'K', '*', ],
	    'G' => ['G', 'S', 'E', '*', 'A', 'V', 'D', 'R', 'W', 'C', ],
	    'Q' => ['E', 'Q', 'K', 'L', 'H', 'R', 'P', '*', ],
	    'T' => ['N', '*', 'P', 'A', 'I', 'T', 'K', 'S', 'M', ],
	    'C' => ['F', 'S', 'R', 'Y', 'C', 'W', 'G', ],
	    'W' => ['*', 'R', 'W', 'C', 'L', 'G', ],
	    'I' => ['V', 'F', 'I', 'L', 'N', 'S', 'M', 'T', ],
	    'L' => ['*', 'F', 'P', 'V', 'R', 'H', 'I', 'W', 'L', 'Q', 'M', ],
	    'Y' => ['D', 'F', 'C', 'H', 'Y', 'N', '*', ],
	    'H' => ['Q', 'N', 'D', 'P', 'L', 'Y', 'H', 'R', ],
	    'R' => ['S', 'G', 'Q', 'P', 'L', 'W', 'C', 'H', 'R', '*', ],
	    'A' => ['E', 'T', 'G', 'V', 'D', 'P', 'A', ],
	    '*' => ['*', 'R', 'L', 'Y', 'W', 'K', 'S', 'T', 'Q', 'G', 'E', 'M', ],
	    'N' => ['D', 'I', 'Y', 'H', 'N', 'S', 'T', 'K', ],
	    'S' => ['R', 'I', 'C', '*', 'N', 'S', 'G', 'T', ],
	    'K' => ['M', 'E', 'K', 'Q', 'T', '*', 'N', ],
	    'X' => ['X'],
	    );
    }
    elsif ($codontable_id == 3) {
	%possible_aa_changes = (
	    'A' => ['D', 'A', 'G', 'T', 'V', 'P', 'E', ],
	    'K' => ['M', 'E', 'R', 'K', '*', 'N', 'T', 'Q', ],
	    'Y' => ['Y', 'H', '*', 'N', 'D', 'F', 'C', ],
	    'D' => ['G', 'E', 'V', 'N', 'D', 'Y', 'H', 'A', ],
	    'F' => ['I', 'Y', 'F', 'T', 'V', 'C', 'L', ],
	    '*' => ['E', 'W', 'L', 'Q', '*', 'K', 'Y', ],
	    'N' => ['I', 'S', 'Y', 'K', 'H', 'N', 'D', 'T', ],
	    'V' => ['E', 'M', 'I', 'V', 'T', 'L', 'G', 'A', 'F', 'D', ],
	    'T' => ['P', 'I', 'S', 'R', 'H', 'L', 'M', 'K', 'A', 'N', 'F', 'V', 'T', 'Q', ],
	    'Q' => ['P', 'K', 'H', '*', 'E', 'T', 'Q', 'R', ],
	    'M' => ['M', 'K', 'I', 'L', 'R', 'T', 'V', ],
	    'W' => ['G', 'C', 'R', 'L', 'W', '*', ],
	    'E' => ['Q', 'G', 'E', 'V', '*', 'D', 'K', 'A', ],
	    'H' => ['R', 'P', 'T', 'Q', 'H', 'Y', 'D', 'N', ],
	    'G' => ['E', 'R', 'W', 'S', 'V', 'C', 'G', 'A', 'D', ],
	    'C' => ['R', 'W', 'S', 'G', 'C', 'Y', 'F', ],
	    'L' => ['L', 'T', 'V', '*', 'F', 'W', 'M', ],
	    'I' => ['S', 'I', 'M', 'F', 'N', 'T', 'V', ],
	    'P' => ['A', 'H', 'P', 'T', 'R', 'Q', ],
	    'S' => ['R', 'I', 'S', 'T', 'G', 'C', 'N', ],
	    'R' => ['H', 'K', 'G', 'C', 'Q', 'T', 'W', 'S', 'P', 'M', 'R', ],
	    );
    }
    elsif ($codontable_id == 5) {
	%possible_aa_changes = (
	    'L' => ['H', 'I', 'R', 'V', 'Q', 'W', 'L', 'P', 'F', '*', 'M', ],
	    'V' => ['D', 'E', 'F', 'I', 'A', 'G', 'V', 'L', 'M', ],
	    'Q' => ['R', 'K', 'P', 'E', 'H', 'L', '*', 'Q', ],
	    'N' => ['K', 'S', 'I', 'D', 'Y', 'T', 'H', 'N', ],
	    'W' => ['S', 'L', 'R', 'G', '*', 'C', 'W', ],
	    'G' => ['G', 'A', 'W', 'V', 'R', 'S', 'C', 'E', 'D', ],
	    'I' => ['M', 'L', 'V', 'N', 'I', 'S', 'F', 'T', ],
	    'S' => ['C', 'T', 'R', 'K', 'I', 'S', 'W', 'G', 'N', 'M', ],
	    'R' => ['C', 'H', 'R', 'S', 'P', 'G', 'W', 'Q', 'L', ],
	    'K' => ['K', 'S', 'E', 'T', 'M', '*', 'Q', 'N', ],
	    'T' => ['A', 'N', 'M', 'T', 'K', 'I', 'P', 'S', ],
	    'H' => ['R', 'P', 'D', 'Y', 'H', 'L', 'Q', 'N', ],
	    'C' => ['Y', 'C', 'W', 'G', 'S', 'R', 'F', ],
	    'D' => ['N', 'V', 'G', 'A', 'H', 'D', 'E', 'Y', ],
	    'E' => ['K', 'E', 'D', '*', 'A', 'G', 'V', 'Q', ],
	    'Y' => ['N', '*', 'F', 'H', 'Y', 'C', 'D', ],
	    'M' => ['V', 'M', 'L', 'T', 'S', 'I', 'K', ],
	    'A' => ['V', 'G', 'A', 'T', 'D', 'E', 'P', ],
	    '*' => ['Q', 'W', '*', 'L', 'Y', 'E', 'K', ],
	    'P' => ['H', 'T', 'R', 'P', 'A', 'Q', 'L', ],
	    'F' => ['I', 'F', 'C', 'Y', 'L', 'V', ],
	    'X' => ['X'],
	    );
    }
    else{
	die "Unrecognised amino acid codon table code: $codontable_id";
    }
    die "Amino acid '$amino_acid' not found in codontable $codontable_id"
	unless exists $possible_aa_changes{$amino_acid};

    my @possible_aas = @{$possible_aa_changes{$amino_acid}};

    return \@possible_aas;
}

sub param {
    my $self = shift;
    
    unless ($self->input_job) {
        # if we don't have an input job, add a dummy one (used when we're not 
        # running as part of a pipeline proper)
        $self->input_job(Bio::EnsEMBL::Hive::AnalysisJob->new);
    }

    return $self->SUPER::param(@_);
}

sub required_param {
    my $self        = shift;
    my $param_name  = shift;
    
    my $param_value = $self->param($param_name, @_);
    
    die "$param_name is a required parameter" unless defined $param_value;
    
    return $param_value;
}

sub run_cmd {
    my $self = shift;
    my $cmd = shift;
    if (my $return_value = system($cmd)) {
	$return_value >>= 8;
	die "system($cmd) failed: $return_value";
    }
}

sub _insert_error_msg {
    my $self = shift;
    my ($translation_md5, $error, $analysis) = @_;
    my $sql = "INSERT INTO failure_reason (translation_md5,error_msg,analysis) VALUES (?,?,?)";
    my $sth = $self->data_dbc->prepare($sql);
    $sth->execute($translation_md5, $error, $analysis);
    $sth->finish();

    my $pf_dsn = 'dbi:mysql:database=' . $self->required_param('pfdb_name') .
	    ';host=' . $self->required_param('pfdb_host') .
	    ';port=' . $self->required_param('pfdb_port');
    my $pf_dbh = DBI->connect($pf_dsn, $self->required_param('pfdb_user'),
			      $self->required_param('pfdb_pass')) or die $DBI::errstr;
    
    my $pf_sql = qq{
             INSERT INTO valid_failure (translation_md5_id, analysis_id) VALUES (
                 (SELECT translation_md5_id FROM translation_md5
                     WHERE translation_md5 = ?),
                 (SELECT analysis_id FROM analysis
                     WHERE analysis = ?)
             );
         };
    my $pf_sth = $pf_dbh->prepare($pf_sql);
    $pf_sth->execute($translation_md5, $analysis);
    $pf_sth->finish();
}


1;
