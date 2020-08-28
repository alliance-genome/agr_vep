package VepProteinFunction::RunPolyPhen2;

use strict;
use warnings;

use File::Copy;
use File::Path qw(make_path remove_tree);
use Data::Dumper;

use base ('VepProteinFunction::BaseVepProteinFunction');

sub run {
    my $self = shift;

    my $translation_md5 = $self->required_param('translation_md5');
    my $pph_dir         = $self->required_param('pph_dir');
    my $working_dir     = $self->required_param('pph_working');
    
    my $dir        = substr($translation_md5, 0, 2);
    my $output_dir = "$working_dir/$dir/$translation_md5";

    my @pph_dirs    = qw{alignments blastfiles profiles structures lock source_files .pph polyphen_tmp};
    my @pipe_dirs   = qw{features errors};
    
    chdir $output_dir or die "Failed to chdir to $output_dir";

    for my $conf_prefix (qw(align databases debug options paths programs programs_options)){
	my $conf_file = $self->required_param('pph_conf_dir') . '/' . $conf_prefix . '.cnf';
	copy($conf_file, $output_dir . '/.pph');
    }

    my $subs_file       = "${output_dir}/source_files/subs.txt";
    my $fasta_file      = "${output_dir}/source_files/protein.fa";
    my $aln_file        = "${output_dir}/alignments/${translation_md5}.aln";
    my $output_file     = "${output_dir}/features/features.txt";
    my $errors_file     = "${output_dir}/errors/errors.log";
   
    # get the protein sequence, write PolyPhen-2 input file

    my ($protein_fasta, $protein_seq) = $self->get_protein_fasta($translation_md5);
    $self->create_pph_subs_file($subs_file, $translation_md5, $protein_seq);
    
    # run polyphen
    $self->dbc->disconnect_when_inactive(1);

    # use -A option to disable polyphen's own LSF support (which conflicts with the hive)
    my $cmd = "$pph_dir/bin/run_pph.pl -A -d $output_dir -s $fasta_file $subs_file 1> $output_file";
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);
    
    $self->dbc->disconnect_when_inactive(0);
    
    if ($exit_code != 0){
	die "Failed to run polyphen for $translation_md5: $stderr";
    }
    
    # delete scratch files unless in debug mode
    unless ($self->param('debug_mode')) {
	my $err;
	remove_tree(@pph_dirs, {error => \$err});
	die "remove_tree failed: ".Dumper($err) if $err && @$err;
    }

    ($exit_code, $stderr, $flat_cmd)=  $self->run_system_command("gzip -f $output_file");
    
    if ($exit_code == 0) {
        $output_file .= '.gz';
    }
    else {
	$self->warning("Failed to gzip $output_file: $stderr"); 
    }

    $self->param('feature_file', $output_file);
}

sub create_pph_subs_file{
    my ($self, $subs_file, $translation_md5, $protein_seq) = @_;

    my $codontable_id = $self->required_param('codontable_id');

    open (SUBS, ">$subs_file") or die "Failed to open file for SUBS $subs_file: $!";
    my $pos = 0;
    for my $aa (split(//, $protein_seq)) {
	$pos++;
	next unless $aa =~ /^[GPAVLIMCFYWHKRQNEDST]$/;
	for my $possible_aa(@{$self->possible_aa_changes($codontable_id, $aa)}){
	    next if $possible_aa eq $aa;
	    next unless $possible_aa =~ /^[GPAVLIMCFYWHKRQNEDST]$/;
	    print SUBS join("\t", $translation_md5, $pos, $aa, $possible_aa) . "\n";
	}
    }
    close SUBS;
};

sub write_output {
    my $self = shift;
    
    if (my $feature_file = $self->param('feature_file')) {
        $self->dataflow_output_id( [{
            translation_md5 => $self->param('translation_md5'),
            feature_file    => $feature_file,
				    }], 2);
    }
}


1;
