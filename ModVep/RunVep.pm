package ModVep::RunVep;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub run {
    my $self = shift;
    
    my $input_file = $self->required_param('vep_input_file');
    my $output_file = $input_file . '.vep.vcf';
    
    my $fasta_file = $self->required_param('fasta');
    my $gff_file = $self->required_param('gff');
    my $map_file = $self->required_param('vep_working') . '/transcript_id_name_map.dmp';

    my $plugin_str = 'ProtFuncAnnotHTP,mod=' . $self->required_param('mod') . ',pass=' .
	$self->required_param('password');

    my $tn_plugin_str = 'TranscriptName,name=' . lc($self->required_param('pipeline_db_name')) . ',host=' . $self->required_param('pipeline_host') .
	',user=' . $self->required_param('pipeline_user') . ',port=' . $self->required_param('pipeline_port') . ',pass=' . $self->required_param('password');


    # Run VEP
    $self->dbc->disconnect_when_inactive(1);
    my $cmd = "vep -i $input_file -gff $gff_file --format vcf -fasta $fasta_file --vcf --hgvs --hgvsg -shift_hgvs=0 --symbol --numbers --distance 0 --output_file $output_file --force_overwrite --plugin $plugin_str --plugin GenomePos --plugin $tn_plugin_str --check_ref --flag_pick_allele_gene --safe";
    
    if ($self->param('bam')) {
	$cmd .= ' --bam ' . $self->param('bam');
    }
	
    $self->param('vep_failure', 1);
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    
    $self->dbc->disconnect_when_inactive(0);

    if ($exit_code != 0) {
	die "$flat_cmd - exit code: $exit_code: $stderr";
    }
    else {
	$self->param('vep_failure', 0);
	unlink $input_file unless $self->required_param('debug');
    }
    system("rm ${output_file}_summary.html");
}


sub post_cleanup {
    my $self = shift;

    my $branch_nr = $self->param('vep_failure') ? 3 : 2;
    $self->dataflow_output_id([{vep_input_file => $self->param('vep_input_file')}], $branch_nr);
}
    
    
1;
