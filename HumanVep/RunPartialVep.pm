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
 
    $self->param('vep_substr_failure', 1);
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    $self->dbc->disconnect_when_inactive(0);

    if ($exit_code != 0) {
	if ($exit_code == 2 or $stderr =~ /outside\sof\sstring\sat.+line 906,\s<\$__ANONIO__>\sline\s\d+\./) {
	    $self->join_results($failed_file, $input_file);
	    $self->warning($flat_cmd . ': ' . $stderr);
	}
	else {
	    die "$flat_cmd - exit code: $exit_code: $stderr" if $exit_code != 0;
	}
    }
    else {
	$self->join_results($failed_file, $input_file);
	$self->param('vep_substr_failure', 0);
	$self->remove_header() unless $failed_file =~ /_1\/chr1_001$/;
	system("rm $failed_file");
    }
   
}

sub create_input_file {
    my ($self, $failed_file) = @_;

    my ($last_pos, $last_ref, $last_alt) = $self->last_vep_result_printed();

    my $input_file = $failed_file . '_part';

    my $cmd;
    if ($last_pos == 0) {
       $cmd = "cp $failed_file $input_file";
    }
    else {
	open (GREP, 'grep -n $"\t' . $last_pos . '\t" ' . $failed_file . ' |') 
	    or die "Couldn't grep file $failed_file";
	my $lines_printed;
	while (<GREP>) {
	    ($lines_printed) = $_ =~ /^(\d+):/;
	}
	close (GREP);

	my $lines_remaining = $self->required_param('lines_per_input_file') - $lines_printed;
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

    my $vep_file = $self->required_param('vep_input_file') . '.vep.vcf';
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

    if ($self->param('vep_substr_failure')) {
	$self->dataflow_output_id([{vep_input_file => $self->param('vep_input_file')}], 2);
    }
}
 

1;
