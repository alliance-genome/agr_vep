package VepProteinFunction::RunWeka;

use strict;
use warnings;

use File::Copy;
use File::Path qw(make_path remove_tree);

use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;

use base ('VepProteinFunction::BaseVepProteinFunction');

sub run {
    my $self = shift;

    my $translation_md5 = $self->required_param('translation_md5');

    my $feature_file = $self->required_param('feature_file');
    
    # copy stuff to /tmp to avoid lustre slowness
    my ($output_dir, $feature_filename) = $feature_file =~ /(.+)\/([^\/]+)$/;
    $self->param('output_dir', $output_dir);
    my $tmp_dir = "/tmp/mqt_weka_${translation_md5}";
    $self->param('tmp_dir', $tmp_dir);
    make_path($tmp_dir);
    copy($feature_file, $tmp_dir);

    my $input_file = "${tmp_dir}/${feature_filename}";

    # unzip the file if necessary
    if ($input_file =~ /\.gz$/ && -e $input_file) {    
	my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command("gunzip -f $input_file");
	die "Failed to gunzip input file: $input_file: $stderr" unless $exit_code == 0;

	$input_file =~ s/.gz$//;
    }

    chdir $output_dir or die "Failed to chdir to $output_dir";

    my @to_delete;
    for my $conf_prefix (qw(align databases debug options paths programs programs_options)){
	push @to_delete, "${output_dir}/${conf_prefix}.cnf";
    }

    $self->run_weka($translation_md5, $input_file, 'humdiv', $output_dir, $tmp_dir);
    $self->run_weka($translation_md5, $input_file, 'humvar', $output_dir, $tmp_dir) if $self->required_param('mod') eq 'human';
}


sub run_weka {
    my ($self, $translation_md5, $input_file, $dataset, $output_dir, $tmp_dir) = @_;

    my $pph_dir = $self->required_param('pph_dir');
    my $output_file = "${tmp_dir}/pph_${dataset}.txt";
    
    my $model = $self->required_param($dataset . '_model');

    $self->dbc->disconnect_when_inactive(1);

    my $cmd = "$pph_dir/bin/run_weka.pl -l $model $input_file 1> $output_file";
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    
    $self->dbc->disconnect_when_inactive(0);
    
    die "Weka run ($dataset) for $translation_md5 failed - cmd: $flat_cmd: $stderr" unless $exit_code == 0;
    
    open (RESULT, "<$output_file") or die "Failed to open output file: $!";

    
    my ($protein_fasta, $protein_seq) = $self->get_protein_fasta($translation_md5);
    
    my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
            -analysis           => 'polyphen',
            -sub_analysis       => $dataset,
            -peptide_length     => length($protein_seq),
            -translation_md5    => $translation_md5,
	    );
    
    my @fields;
    my $results_available = 0;
    while (<RESULT>) {
	if (/^#/) {
	    s/#//g;
	    @fields = split /\s+/;
	    next;
	}

	die "No header line in result file $output_file?" unless @fields; 

	my @values = split /\t/;

	# trim whitespace
	map { $_ =~ s/^\s+//; $_ =~ s/\s+$// } @values; 
            
	# parse the results into a hash
	my %results = map { $fields[$_] => $values[$_] } (0 .. @fields-1);
            
	my $alt_aa      = $results{o_aa2};
	my $prediction  = $results{prediction};
	my $prob        = $results{pph2_prob};
	my $position    = $results{o_pos};

	next unless $position && $alt_aa;

	$results_available = 1;

	$pred_matrix->add_prediction(
	    $position,
	    $alt_aa,
	    $prediction, 
	    $prob,
	    );
    }
    
    my $pf_dsn = 'dbi:mysql:database=' . $self->required_param('pfdb_name') .
	';host=' . $self->required_param('pfdb_host') .
	';port=' . $self->required_param('pfdb_port');
    my $pf_dbh = DBI->connect($pf_dsn, $self->required_param('pfdb_user'),
			      $self->required_param('pfdb_pass')) or die $DBI::errstr;

    if ($results_available == 1 ){  
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
	my $analysis = $dataset eq 'humdiv' ? 'pph' : 'pph_humvar';
	$matrix_insert_sth->execute($translation_md5, $analysis, $pred_matrix->serialize);
	$self->dbc and $self->dbc->disconnect_if_idle();
    }

    copy($output_file, $output_dir) if $self->param('debug_mode');

    return;
}

sub post_cleanup {
    my $self = shift;
    
    my $tmp_dir = $self->param('tmp_dir');
    my $output_dir = $self->param('output_dir');
    my $pph_root_dir = $self->param('pph_working');

    my $err;
    remove_tree($tmp_dir, {error=> \$err});
    die "remove_tree failed: ".Dumper($err) if $err && @$err;    	
    
    unless ($self->param('debug_mode')){
	chdir $pph_root_dir or die "Failed to chdir to $pph_root_dir";
	remove_tree($output_dir, {error => \$err});
	die "remove_tree failed: ".Dumper($err) if $err && @$err;
    }
}

1;
