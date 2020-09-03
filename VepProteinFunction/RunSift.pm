package VepProteinFunction::RunSift;

use strict;

use File::Copy;
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
    $self->param('output_dir', $output_dir);

    unless (-d $output_dir) {
	my $err;
	make_path($output_dir, {error => \$err});
	die "make_path failed: ".Dumper($err) if $err && @$err;
    }
    
    chdir $output_dir or die "Failed to chdir to $output_dir";
    
    my $fasta_file = $translation_md5 . '.fa';
    my $aln_file = $translation_md5 . '.alignedfasta';
    my $subst_file = $translation_md5 . '.subst';
    my $res_file   = $translation_md5 . '.SIFTprediction';

    # Set necessary environment variables for SIFT
    $ENV{tmpdir} = $output_dir;
    $ENV{NCBI} = $self->required_param('ncbi_dir');
    $ENV{SIFT_DIR} = $sift_dir;
    $ENV{BLIMPS_DIR} = $sift_dir . '/blimps';

    # Fetch our protein 
    my ($protein_fasta, $protein_seq) = $self->get_protein_fasta($translation_md5);
    $protein_fasta =~ s/\*/X/g;
    $protein_seq =~ s/\*/X/g;
    $self->create_subst_file($subst_file, $protein_seq);

    $self->dbc and $self->dbc->disconnect_if_idle();

    open (FASTA_FILE, ">$fasta_file") or die "Failed to open protein FASTA $fasta_file: $!";;
    print FASTA_FILE $protein_fasta;
    close FASTA_FILE;

    copy($fasta_file, $translation_md5 . '.unfiltered');

    my @commands = (
	"psiblast -db $blastdb -query $fasta_file -out $translation_md5.out -outfmt 0 -num_iterations 1 -evalue 0.0001 -num_descriptions 399 -num_alignments 399",
	"$sift_dir/bin/psiblast_res_to_fasta_dbpairwise $translation_md5.out $translation_md5.globalX 1 $translation_md5.unfiltered",
	"$sift_dir/bin/clump_output_alignedseq $translation_md5.globalX $translation_md5.clumped .9 0",
	"$sift_dir/bin/choose_seqs_via_psiblastseedmedian $translation_md5.unfiltered $translation_md5.clumped $translation_md5.selectedclumped 1.0 $translation_md5 2.75",
	"$sift_dir/bin/consensus_to_seq $translation_md5.selectedclumped $translation_md5.clumped.consensuskey $translation_md5.selected",
	"$sift_dir/bin/seqs_from_psiblast_res $translation_md5.out $translation_md5.selected 1 $translation_md5.unfiltered $translation_md5.alignedfasta $translation_md5",
	"$sift_dir/bin/info_on_seqs $translation_md5.alignedfasta $translation_md5.subst $translation_md5.SIFTprediction"
	);

    # Run Sift
    $self->dbc->disconnect_when_inactive(1);
    my ($exit_code, $stderr, $flat_cmd);
    $exit_code = 0;
    for my $cmd (@commands) {
	last if $exit_code != 0;
	($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd, {timeout => 36000});
    }
    $self->dbc->disconnect_when_inactive(0);
  
    if ($exit_code != 0) {
	$self->handle_error($translation_md5, $exit_code, $stderr, $flat_cmd);
	return;
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
}

sub post_cleanup {
    my $self = shift;
    
    # delete files unless in debug mode
    unless ($self->param('debug_mode')) {
	my $working_dir = $self->param('sift_working');
	my $translation_md5 = $self->param('translation_md5');
	my $output_dir = $self->param('output_dir');

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
	    print SUBST $aa . $pos . $possible_aa . "\n" unless $possible_aa eq '*' or $possible_aa eq 'X';
	} 
    }
    close SUBST;
}

sub handle_error {
    my ($self, $translation_md5, $exit_code, $stderr, $flat_cmd) = @_;

    my $error_file = $translation_md5 . '.clumped.error';
    if (-s $error_file) {
	open my $error_fh, "<", $error_file or die $!;
	my $error_msg = <$error_fh>;
	close $error_fh;
	chomp $error_msg;
	if ($error_msg =~ /After clustering at 90%, only \d+ cluster\(s\) remained./) {
	    $self->_insert_error_msg($translation_md5, 'Insufficient sequence diversity', 'sift');
	    return;
	}
    }

    $error_file = $translation_md5 . '.alignedfasta.error';
    if (-s $error_file) {
	open my $error_fh, "<", $error_file or die $!;
	my $error_msg = <$error_fh>;
	close $error_fh;
	chomp $error_msg;
	if ($error_msg =~ /\d sequence\(s\) were chosen\. Less than the minimum number of sequences required/) {
	    $self->_insert_error_msg($translation_md5, 'Less than the minimum number of sequences required', 'sift');
	    return;
	} 
    }

    $error_file = $translation_md5 . ".globalX.error";
    if (-s $error_file) {
	open my $error_fh, "<", $error_file or die $!;
	my $error_msg = <$error_fh>;
	close $error_fh;
	chomp $error_msg;
	if ($error_msg =~ /Not enough sequences \(only \d\) found by the PSI-BLAST search!|PSI-BLAST found no hits/) {
	    $self->_insert_error_msg($translation_md5, 'Not enough sequences found by PSI-BLAST search', 'sift');
	    return;
	}
    }
    
    die "$flat_cmd - exit code: $exit_code: $stderr";
}
    
1;
