package ModVep::RunPartialVep;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');


sub run {
    my $self = shift;

    my $failed_file = $self->required_param('failed_input_file');
    my $input_file = $self->create_input_file($failed_file);
    my $output_file = $input_file . '.vep.vcf';

    my $fasta_file = $self->required_param('fasta');
    my $gff_file = $self->required_param('gff');

    my $plugin_str = 'ProtFuncAnnot,mod=' . $self->required_param('mod') . ',pass=' .
	$self->required_param('password');


    # Run VEP on new file
        my $cmd = "vep -i $input_file -gff $gff_file -fasta $fasta_file --vcf --hgvs --hgvsg -shift_hgvs=0 --symbol --no_intergenic --distance 0 --output_file $output_file --force_overwrite --plugin $plugin_str";
    if ($self->param('bam')) {
	$cmd .= ' --bam ' . $self->param('bam');
    }

    $self->param('vep_failure', 1);
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    
    $self->dbc->disconnect_when_inactive(0);

    if ($exit_code != 0) {
	$self->param('vep_failure_no', $self->required_param('vep_failure_no') + 1);
	if ($self->param('vep_failure_no') == 3) {
	    die "Too many failures for $failed_file";
	}
	elsif ($exit_code == 2) {
	    $self->join_results($failed_file, $input_file);
	    $self->warning($flat_cmd . ': ' . $stderr);
	}
	else {
	    die "$flat_cmd - exit code: $exit_code: $stderr" if $exit_code != 0;
	}
    }
    else {
	$self->join_results($failed_file, $input_file);
	$self->param('vep_failure', 0);
	$self->remove_header() unless $failed_file =~ /_0+1$/;
	system("rm $failed_file");
    }
}


sub create_input_file {
    my ($self, $failed_file) = @_;

    my ($last_pos, $last_ref, $last_alt) = $self->last_vep_result_printed();

    my $last_vep_success = join('|', $last_pos, $last_ref, $last_alt);
    if (defined $self->param('last_vep_success')) {
	if ($self->param('last_vep_success') ne $last_vep_success) {
	    $self->param('vep_failure_no', 1);
	}
    }
    $self->param('last_vep_success', $last_vep_success);

    my $input_file = $failed_file . '_' . $last_pos . '_' . $self->param('vep_failure_no');
    open (INPUT, '>', $input_file) or die "Could not open $input_file for writing";
    open (FAILED, '<', $failed_file) or die "Could not read $failed_file";

    my $last_line_seen = 0;
    while (<FAILED>) {
	chomp;
	my @columns = split("\t", $_);
	if ($columns[1] == $last_pos and $columns[3] eq $last_ref and $columns[4] eq $last_alt) {
	    $last_line_seen = 1;
	    next;
	}
	next unless $last_line_seen;
	
	print INPUT "$_\n";
    }
    close (INPUT);
    close (FAILED);

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

    my ($last_pos, $last_ref, $last_alt);

    my $vep_file = $self->required_param('failed_input_file') . '.vep.vcf';
    open(VEP, "tail -n 1 $vep_file |") or die "Couldn't open $vep_file";
    while (<VEP>) {
	my @columns = split("\t", $_);
	$last_pos = $columns[1];
	$last_ref = $columns[3];
	$last_alt = $columns[4];
    }
    close(VEP);

    return ($last_pos, $last_ref, $last_alt);
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

    if ($self->param('vep_substr_failure') > 0) {
	$self->dataflow_output_id([{failed_input_file => $self->param('failed_input_file'),
				    vep_failure_no => $self->param('vep_failure_no')}], 2);
    }
}
 

1;
