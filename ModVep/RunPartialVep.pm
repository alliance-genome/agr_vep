package ModVep::RunPartialVep;

use strict;
use warnings;

use Path::Class;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');


sub run {
    my $self = shift;

    my $failed_file = $self->required_param('vep_input_file');
    my $input_file = $self->create_input_file($failed_file);
    my $output_file = $input_file . '.vep.vcf';

    my $fasta_file = $self->required_param('fasta');
    my $gff_file = $self->required_param('gff');

    my $plugin_str = 'ProtFuncAnnotHTP,mod=' . $self->required_param('mod') . ',pass=' .
	$self->required_param('password');


    # Run VEP on new file
    $self->dbc->disconnect_when_inactive(1);
    my $cmd = "vep -i $input_file -gff $gff_file --format vcf -fasta $fasta_file --vcf --hgvs --hgvsg -shift_hgvs=0 --symbol --numbers --distance 0 --output_file $output_file --force_overwrite --plugin $plugin_str --plugin GenomePos --plugin TranscriptName,gff=$gff_file --check_ref --flag_pick_allele_gene" ;
    if ($self->param('bam')) {
	$cmd .= ' --bam ' . $self->param('bam');
    }

    $self->param('vep_failure', 1);
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    $self->dbc->disconnect_when_inactive(0);

    $self->join_results($failed_file, $input_file);

    if ($exit_code != 0) {
	die "$flat_cmd - exit code: $exit_code: $stderr" if $exit_code != 0;
    }
    else {
	$self->param('vep_failure', 0);
	unlink $failed_file;
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
    
    my @cmds = ("touch $failed_results_file",
	        "grep '^#' $input_results_file > $input_results_file.header",
		"grep -v '^#' $failed_results_file > $failed_results_file.body",
		"grep -v '^#' $input_results_file > $input_results_file.body",
		"cat $input_results_file.header $failed_results_file.body $input_results_file.body > $failed_results_file",
		"rm $input_file* $failed_results_file.*");

    for my $cmd (@cmds) {
	my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
	die "ERROR: $exit_code: $flat_cmd - $stderr" unless $exit_code == 0 or ($exit_code == 1 and $cmd =~ /^grep/);
    }

    return;
}


sub last_vep_result_printed {
    my $self = shift;

    my $last_pos = 0;
    my $last_ref = 'REF';
    my $last_alt = 'ALT';
    my $last_id = '';
    my $last_line = '';

    my ($exit_code, $stderr, $flat_cmd);

    my $vep_file = $self->required_param('vep_input_file') . '.vep.vcf';
    
    if (-e $vep_file) {
	my $vep_fh = file($vep_file)->openr;
	my $tmp_fh = file($vep_file . '.tmp')->openw;
	while (my $line = $vep_fh->getline) {
	    last if $line !~ /\n$/;
	    $tmp_fh->print($line);
	    next if $line =~ /^#/;
	    my @columns = split("\t", $line);
	    $last_pos = $columns[1];
	    $last_id = $columns[2];
	    $last_ref = $columns[3];
	    $last_alt = $columns[4];
	}
	
	($exit_code, $stderr, $flat_cmd) =
	    $self->run_system_command("mv $vep_file.tmp $vep_file");
    }
    else {
	($exit_code, $stderr, $flat_cmd) =
	    $self->run_system_command("touch $vep_file");
    }
    die "$flat_cmd: $exit_code: $stderr" unless $exit_code == 0;


    return ($last_pos, $last_id, $last_ref, $last_alt);
}


sub post_cleanup {
    my $self = shift;

    my $branch_nr = $self->param('vep_failure') ? 3 : 2; 
    $self->dataflow_output_id([{vep_input_file => $self->param('vep_input_file')}], $branch_nr);
}
 

1;
