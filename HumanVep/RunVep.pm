package HumanVep::RunVep;

use strict;

use File::Copy;
use File::Path qw(make_path remove_tree);
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub run {
    my $self = shift;
    
    my $input_file = $self->required_param('vep_input_file');
    my ($chr) = $input_file =~ /\/chr([^_\/]+)_\d{3}$/;
    
    my $fasta_file = $self->required_param('rs_fasta_prefix') . $chr .
	$self->required_param('rs_fasta_suffix');
    my $cache = $self->required_param('vep_cache');
    my $cache_dir = $self->required_param('vep_cache_dir');

    # Run VEP
    $self->dbc->disconnect_when_inactive(1);
    my $cmd = "vep --cache --offline --$cache -i $input_file --fasta $fasta_file --vcf --sift b --polyphen b --hgvs --hgvsg --dir_cache $cache_dir --output_file $input_file.vep.vcf --force_overwrite --symbol --no_intergenic --distance 0";
    
    $self->param('vep_failure', 1);
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    $self->dbc->disconnect_when_inactive(0);

    if ($exit_code != 0) {
	if ($exit_code == 2 or $stderr =~ /outside\sof\sstring\sat.+line 906,\s<\$__ANONIO__>\sline\s\d+\./) {
	    $self->warning($flat_cmd . ': ' . $stderr);
	}
	else {
	    die "$flat_cmd - exit code: $exit_code: $stderr" if $exit_code != 0;
	}
    }
    else {
	$self->param('vep_failure', 0);
	$self->remove_header() unless $chr eq '1' and $input_file =~ /_001$/;
	system("rm $input_file");
    }
    system("rm $input_file.vep.vcf_summary.html");
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


sub write_output {
    my $self = shift;

    if ($self->param('vep_failure')) {
	$self->dataflow_output_id([{vep_input_file => $self->param('vep_input_file')}], 2);
    }
}
    
    
1;
