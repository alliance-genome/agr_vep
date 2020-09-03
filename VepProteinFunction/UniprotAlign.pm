package VepProteinFunction::UniprotAlign;

use strict;
use warnings;

use File::Path qw(make_path);
use Data::Dumper;

use base ('VepProteinFunction::BaseVepProteinFunction');

sub run {
    my $self = shift;

    my $translation_md5 = $self->required_param('translation_md5');
    my $working_dir = $self->required_param('pph_working');

    # Setup folders for this and later PolyPhen2 run
    my $dir = substr($translation_md5, 0, 2);
    my $output_dir = "$working_dir/$dir/$translation_md5";

    unless (-d $output_dir) {
	my $err;
	make_path($output_dir, {error => \$err});
	die "make_path failed: ".Dumper($err) if $err && @$err;
    }
    
    chdir $output_dir or die "Failed to chdir to $output_dir";

    my @pph_dirs    = qw{alignments blastfiles profiles structures lock source_files .pph polyphen_tmp};
    my @pipe_dirs   = qw{features errors};

    my $err;
    make_path(@pph_dirs, @pipe_dirs, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;

    # Retrieve and dump the protein sequence
    my $fasta_file      = "${output_dir}/source_files/protein.fa";
    my ($protein_fasta, $protein_seq) = $self->get_protein_fasta($translation_md5);
    open (FASTA_FILE, ">$fasta_file") or die "Failed to open protein $fasta_file for writing: $!";
    print FASTA_FILE $protein_fasta;
    close FASTA_FILE;

    # Check if sequence is (sub)sequence in Uniprot FASTA
    my ($acc, $uniprot_seq) = $self->grep_uniprot_fasta($protein_seq);

    # If Uniprot sequence found, dump sequence and set as BLAST target
    my $blast_query = $fasta_file;
    my $blast_out = "$output_dir/blastfiles/$translation_md5.blast";
    if (defined $acc){
	$blast_out = "$output_dir/blastfiles/$acc.blast";
	$blast_query = "$output_dir/polyphen_tmp/$acc.fasta";
	open (QUERY_FILE, ">$blast_query") or die "Failed to open $blast_query for writing: $!";
	print QUERY_FILE ">$acc\n$uniprot_seq";
	close QUERY_FILE;
    }
    
    if (-e $blast_out) {
	my $blast_done = 0;
	open (BLAST_TAIL, "tail -n2 $blast_out | head -n1 |");
	while (<BLAST_TAIL>) {
	    chomp;
	    $blast_done = 1 if $_ =~ /<\/BlastOutput>/;
	}
	close BLAST_TAIL;

	return if $blast_done;
    }

    # Carry out BLAST alignment of seq against Uniref100 DB
    $self->dbc->disconnect_when_inactive(1);

    my $uniref_db = $self->required_param('pph_dir') . '/nrdb/uniref100';
    my $blast_cmd = "blastp -seg yes -evalue 1e-3 -num_threads 1 -max_target_seqs 1000 -outfmt 5 -db $uniref_db -query $blast_query -out $blast_out";
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($blast_cmd, {timeout => 36000});

    $self->dbc->disconnect_when_inactive(0);

    if ($exit_code != 0) {
	die "Failed to run BLAST against Uniref100 for $translation_md5 (exit code $exit_code):$stderr";
    }
}


sub grep_uniprot_fasta{
    my ($self, $protein_seq) = @_;

    my $uniprot_fasta = $self->uniprot_fasta_filename;
    open GREP, "grep -s -m 1 -B 1 '$protein_seq' $uniprot_fasta |";
    my $header = <GREP>;
    close GREP, return () unless defined $header;
    my $uniprot_seq = <GREP>;
    close GREP;

    chomp $header;
    chomp $uniprot_seq;
    
    my ($db, $acc, $rest) = split(/\|/, $header);
    
    return ($acc, $uniprot_seq);
}

sub uniprot_fasta_filename{
    my $self = shift;

    my $mod = $self->required_param('mod');
    
    my $common_name;
    if ($mod eq 'FB') {
	$common_name = 'fruitfly';
    }
    elsif ($mod eq 'MGI') {
	$common_name = 'mouse';
    }
    elsif ($mod eq 'RGD') {
	$common_name = 'rat';
    }
    elsif ($mod eq 'SGD') {
	$common_name = 'yeast';
    }
    elsif ($mod eq 'WB') {
	$common_name = 'roundworm';
    }
    elsif ($mod eq 'ZFIN') {
	$common_name = 'zebrafish';
    }
    elsif ($mod eq 'human') {
	$common_name = 'human';
    }
    else {
	die "Unrecognised MOD: $mod";
    }

    my $filename = $self->required_param('pph_dir') . '/uniprot/' . $common_name . '.seq';
    unless (-e $filename) {
	die "Uniprot database does not exist: $filename";
    }

    return $filename;
}

sub write_output {
    my $self = shift;

    $self->dataflow_output_id( [{
	translation_md5 => $self->param('translation_md5'),
	codontable_id   => $self->param('codontable_id'),
				}], 2);
}

1;
