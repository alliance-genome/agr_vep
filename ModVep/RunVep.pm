package ModVep::RunVep;

use strict;
use warnings;

use Try::Tiny;
use IPC::System::Simple qw(system);

use base ('ModVep::BaseModVep');

sub run {
    my $self = shift;
    
    my $input_file = $self->required_param('vep_input_file');
    my $output_file = $input_file . '.vep.vcf';
    
    my $fasta_file = $self->required_param('fasta');
    my $gff_file = $self->required_param('gff');
    
    $self->convert_vcf_chromosomes($self->required_param('mod'), $input_file,
				   $self->required_param('refseq_chromosomes'));

    my $plugin_str = 'ProtFuncTranscriptNameHTP,mod=' . $self->required_param('mod') . ',pass=' .
	$self->required_param('password');

    # Run VEP
    $self->dbc->disconnect_when_inactive(1);
    my $cmd = "vep -i $input_file -gff $gff_file --format vcf -fasta $fasta_file --vcf --hgvs --hgvsg -shift_hgvs=0 --symbol --numbers --distance 0 --output_file $output_file --force_overwrite --plugin $plugin_str --plugin GenomePos --check_ref --flag_pick_allele_gene --safe";
    
    if ($self->param('bam')) {
	$cmd .= ' --bam ' . $self->param('bam');
    }
	
    $self->param('vep_failure', 0);
    try {
	system($cmd);
    }
    catch {
	$self->param('vep_failure', 1);
	$self->warning($_);
    };
    $self->dbc->disconnect_when_inactive(0);

    unlink $input_file unless $self->param('vep_failure') or $self->required_param('debug');
    unlink "${output_file}_summary.html";
}


sub write_output {
    my $self = shift;

    my $branch_nr = $self->param('vep_failure') ? 3 : 2;
    $self->dataflow_output_id([{vep_input_file => $self->param('vep_input_file')}], $branch_nr);
}
    
    
1;
