package HumanVep::RunPartialVep;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub run {
    my $self = shift;

    my $failed_file = $self->required_param('vep_input_file');
    my $input_file = $self->create_input_file($failed_file);

    my ($chr) = $failed_file =~ /\/chr([^_\/]+)_\d{3}$/;
    my $fasta_file = $self->required_param('rs_fasta_prefix') . $chr .
	$self->required_param('rs_fasta_suffix');
    my $cache = $self->required_param('vep_cache');
    my $cache_dir = $self->required_param('vep_cache_dir');

    # Run VEP on new file
     $self->dbc->disconnect_when_inactive(1);
   
     my $cmd = "vep --cache --offline --$cache -i $input_file --fasta $fasta_file --vcf --sift b --polyphen b --hgvs --hgvsg --dir_cache $cache_dir --output_file $input_file.vep.vcf --force_overwrite --symbol --no_intergenic --distance 0";
 
    $self->param('vep_failure', 1);
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    $self->dbc->disconnect_when_inactive(0);

    $self->join_results($failed_file, $input_file);
	    
    if ($exit_code != 0) {
	die "$flat_cmd - exit code: $exit_code: $stderr";
    }
    else {
	$self->param('vep_failure', 0);
	$self->remove_header() unless $chr eq '1' and $failed_file =~ /_001$/;
	system("rm $failed_file");
    }
   
}

sub create_input_file {
    my ($self, $failed_file) = @_;

    my ($last_pos, $last_id, $last_ref, $last_alt) = $self->last_vep_result_printed();

    my $input_file = $failed_file . '_part';
    for my $file ($input_file, $input_file . '.vep.vcf') {
	system("rm $file") if -e $file;
    }

    my $cmd;
    if ($last_pos == 0) {
       $cmd = "cp $failed_file $input_file";
    }
    else {
	my $grep_cmd = 'grep -n $' . "'" . '\t' . $last_pos . '\t' . $last_id .
	    '\t' . $last_ref . '\t' . $last_alt . '\t' . "' " . $failed_file . ' |';
	open (GREP, $grep_cmd) 
	    or die "Couldn't grep file $failed_file for $last_pos, $last_id, $last_ref, $last_alt";
	my $lines_printed;
	while (<GREP>) {
	    ($lines_printed) = $_ =~ /^(\d+):/;
	}
	close (GREP);
	
	my $input_lines = `wc -l < $failed_file`;
	die "Line count for $failed_file failed: $?" if $?;
	my $lines_remaining = $input_lines - $lines_printed;

	$cmd = "tail -n $lines_remaining $failed_file > $input_file";
    }
    
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    die "Couldn't create input file: $flat_cmd: $exit_code: $stderr" if $exit_code != 0;

    return $input_file;
}

sub join_results {
    my ($self, $failed_file, $input_file) = @_;

    my $failed_results_file = $failed_file . '.vep.vcf';
    my $input_results_file = $input_file . '.vep.vcf';

    system("grep -v '^#' $input_results_file > $input_results_file.copy");
		
    my $cmd = "cat $input_results_file.copy >> $failed_results_file";
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    die "Couldn't concatenate $failed_results_file and $input_results_file: $flat_cmd - $stderr"
	unless $exit_code == 0;
		
    system("rm $input_file*");
    
    return;
}

sub last_vep_result_printed {
    my $self = shift;

    my $last_pos = 0;
    my $last_ref = 'REF';
    my $last_alt = 'ALT';
    my $last_id = '';
    my $last_line = '';

    my $vep_file = $self->required_param('vep_input_file') . '.vep.vcf';
    open(VEP, "<$vep_file") or die "Couldn't open $vep_file";
    open(TMP, ">$vep_file.tmp");
    while (<VEP>) {
	last if $_ !~ /\n$/;
	$last_line = $_;
	print TMP $_;
	next if $_ =~ /^#/;
	my @columns = split("\t", $_);
	$last_pos = $columns[1];
	$last_id = $columns[2];
	$last_ref = $columns[3];
	$last_alt = $columns[4];
    }
    close(VEP);
    close(TMP);

    my @cmds = ("mv $vep_file.tmp $vep_file");
    for my $cmd (@cmds) {
	my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
	die "$flat_cmd: $exit_code: $stderr" unless $exit_code == 0;
    }

    return ($last_pos, $last_id, $last_ref, $last_alt);
}

sub remove_header {
    my $self = shift;

    my $file = $self->required_param('vep_input_file') . '.vep.vcf';
    my @cmds = ("grep -v '^#' $file > $file.tmp", "mv $file.tmp $file");
    for my $cmd (@cmds) {
	my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
	die "Couldn't remove header from $file: $exit_code: $stderr" unless $exit_code == 0;
    }

    return;
}

sub write_output {
    my $self = shift;

    if ($self->param('vep_failure')) {
	$self->dataflow_output_id([{vep_input_file => $self->param('vep_input_file')}], 2);
    }
}
 

1;
