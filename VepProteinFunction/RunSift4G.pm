package VepProteinFunction::RunSift4G;

use strict;

use File::Path qw(make_path remove_tree);
use Data::Dumper;

use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;

use base ('VepProteinFunction::BaseVepProteinFunction');

sub run {
    my $self = shift;
    
    my $translation_md5 = $self->required_param('translation_md5');
    my $sift_dir        = $self->required_param('sift_dir');
    my $working_dir     = $self->required_param('sift_working');
    my $blastdb         = $self->required_param('blastdb');

    my $dir = substr($translation_md5, 0, 2);
    my $output_dir = "$working_dir/$dir/$translation_md5";
  
    unless (-d $output_dir) {
	my $err;
	make_path($output_dir, {error => \$err});
	die "make_path failed: ".Dumper($err) if $err && @$err;
    }
    
    chdir $output_dir or die "Failed to chdir to $output_dir";
    
    my $fasta_file = $translation_md5 . '.fa';
    my $subst_file = $translation_md5 . '.subst';
    my $res_file   = $translation_md5 . '.SIFTprediction';

    # Fetch our protein 
    my ($protein_fasta, $protein_seq) = $self->get_protein_fasta($translation_md5);
    $self->create_subst_file($subst_file, $protein_seq);

    $self->dbc and $self->dbc->disconnect_if_idle();

    open (FASTA_FILE, ">$fasta_file") or die "Failed to open protein FASTA $fasta_file: $!";;
    print FASTA_FILE $protein_fasta;
    close FASTA_FILE;

    # Run Sift4G
    $self->param('sift_job_failed', 0);
    $self->dbc->disconnect_when_inactive(1);

    my $cmd = "$sift_dir/sift4g -d $blastdb -q $fasta_file";
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);

    $self->dbc->disconnect_when_inactive(0);
  
    if ($exit_code != 0){
	if ($stderr =~ /No valid queries to process./) {
	    my $err_str = 'Invalid query: ' . $protein_seq;
	    $err_str = substr( $err_str, 0, 497) . '...' if length $err_str > 500;
	    $self->_insert_error_msg($translation_md5, $err_str, 'sift');
	    return;
	}
	else{
#	    $self->param('sift_job_failed', 1);
						  
	    my ($stderr_tail) = $stderr =~ /\*+([^\*]+\*+[^\*]*)$/;
	    $stderr_tail = $stderr unless defined $stderr_tail;
#	    $self->warning("SIFT alignment for $translation_md5 failed, exit code: $exit_code - $stderr_tail");
#	    return;
	    die ("SIFT alignment for $translation_md5 failed, exit code: $exit_code: $exit_code - $stderr_tail");
	}
    }
    
    if (! -e $res_file){
	if (length $protein_seq <= 20){
	    $self->_insert_error_msg($translation_md5,
				     "No results file created - query too short",
				     'sift');
	    return;
	}
	else{
	    die "No SIFT.predictions file created for $translation_md5";
	}
    }
				 	

    # parse and store the results 
	
    open (RESULTS, "<$res_file") or die "Failed to open $res_file: $!";

    my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
	-analysis       => 'sift',
	-peptide_length   => length $protein_seq,
	-translation_md5  => $translation_md5,
	);
    
    my %evidence_stored;
    my $results_available = 0;
    my $pos = 0;
    while (<RESULTS>) {
	chomp;
	my $line = $_;
	
	next if /WARNING/;
	next if /NOT SCORED/;
	
	my ($subst, $prediction, $score, $median_cons, $num_seqs, $blocks) = split("\t", $line);
	
	my ($ref_aa, $pos, $alt_aa) = $subst =~ /([A-Z])(\d+)([A-Z])/;
	
	next unless $ref_aa && $alt_aa && defined $pos;
	
	$results_available = 1;
	
	my $low_quality = 0;
	$low_quality = 1 if $median_cons > 3.25 || $num_seqs < 10;
	
	$pred_matrix->add_prediction(
	    $pos,
	    $alt_aa,
	    $prediction, 
	    $score,
	    $low_quality
	    );
	unless ($evidence_stored{$pos} == 1) {
	    ## add attribs by position
	    $pred_matrix->add_evidence('sequence_number',  $pos, $num_seqs);
	    $pred_matrix->add_evidence('conservation_score', $pos, $median_cons);
	    $evidence_stored{$pos} = 1;
	}
    }
    close RESULTS;
    
    if ($results_available == 1 ){  
	# Connect to results database
	my $pf_dsn = 'dbi:mysql:database=' . $self->required_param('pfdb_name') .
	    ';host=' . $self->required_param('pfdb_host') .
	    ';port=' . $self->required_param('pfdb_port');
	my $pf_dbh = DBI->connect($pf_dsn, $self->required_param('pfdb_user'),
			      $self->required_param('pfdb_pass')) or die $DBI::errstr;

	
	# avoid entering null matrices
	my $matrix_insert_query = qq{
                INSERT INTO protein_function_prediction(translation_md5_id, analysis_id,
                                                        prediction_matrix)
                VALUES(
                    (SELECT translation_md5_id FROM translation_md5
                        WHERE translation_md5 = ?),
                    (SELECT analysis_id FROM analysis
                        WHERE analysis = ?),
                    ?
                )
            };
	my $matrix_insert_sth = $pf_dbh->prepare($matrix_insert_query);
	$matrix_insert_sth->execute($translation_md5, 'sift', $pred_matrix->serialize);
	$self->dbc and $self->dbc->disconnect_if_idle();
    }

    # delete files unless in debug mode
    unless ($self->param('debug_mode')) {
	chdir $working_dir or die "Failed to chdir to $working_dir";
	my $err;
	remove_tree($output_dir, {error => \$err});
	die "remove_tree failed: ".Dumper($err) if $err && @$err;
    }
}

sub create_subst_file {
    my ($self, $subst_file, $protein_seq) = @_;

    my $codontable_id = $self->required_param('codontable_id');

    open SUBST, ">$subst_file" or die "Failed to open $subst_file for writing";
    my $pos = 0;
    for my $aa (split(//, $protein_seq)) {
	$pos++;
	for my $possible_aa(@{$self->possible_aa_changes($codontable_id, $aa)}){
	    print SUBST $aa . $pos . $possible_aa . "\n";
	} 
    }
    close SUBST;
}

sub write_output {
    my $self = shift;

#    if ($self->param('sift_job_failed')) {
#	$self->dataflow_output_id([{translation_md5 => $self->param('translation_md5'),
#				    codontable_id   => $self->param('codontable_id'),}], 3);
#   }

}


1;
