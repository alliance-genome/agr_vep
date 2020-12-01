package ModVep::RunVep;

use strict;

use File::Copy;
use File::Path qw(make_path remove_tree);
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub run {
    my $self = shift;
    
    my $input_file = $self->required_param('vep_input_file');
    my $output_file = $input_file . '.vep.vcf';
    
    my $fasta_file = $self->required_param('fasta');
    my $gff_file = $self->required_param('gff');

    my $plugin_str = 'ProtFuncAnnot,mod=' . $self->required_param('mod') . ',pass=' .
	$self->required_param('password');

    # Run VEP
    $self->dbc->disconnect_when_inactive(1);
    my $cmd = "vep -i $input_file -gff $gff_file -fasta $fasta_file --vcf --hgvs --hgvsg -shift_hgvs=0 --symbol --distance 0 --output_file $output_file --force_overwrite --plugin $plugin_str";
    
    if ($self->param('bam')) {
	$cmd .= ' --bam ' . $self->param('bam');
    }
	
    
    $self->param('vep_failure', 1);
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    
    $self->dbc->disconnect_when_inactive(0);

    if ($exit_code != 0) {
	if ($exit_code == 2) {
	    $self->warning($flat_cmd . ': ' . $stderr);
	}
	else {
	    die "$flat_cmd - exit code: $exit_code: $stderr";
	}
    }
    else {
	$self->param('vep_failure', 0);
	$self->remove_header() unless $input_file =~ /_0+1$/;
	system("rm $input_file");
    }
    system("rm ${output_file}_summary.html");
}


sub remove_header {
    my $self = shift;

    my $file = $self->required_param('vep_input_file') . '.vep.vcf';
    my $tmp_file = $file . '.tmp';
    open (IN, '<', $file) or die $!;
    open (OUT '>', $tmp_file) or die $!;
    while (<IN>) {
	print OUT $_ unless $_ =~ /^#/;
    }
    close (IN);
    close (OUT);
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command("mv $tmp_file $file");
    die "Couldn't remove header from $file: $exit_code: $stderr" unless $exit_code == 0;

    return;
}


sub post_cleanup {
    my $self = shift;

    if ($self->param('vep_failure')) {
	$self->dataflow_output_id([{failed_input_file => $self->param('vep_input_file')}], 2);
    }
}
    
    
1;
