package HumanVep::SplitInput;

use strict;

use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub fetch_input {
    my $self = shift;
    
    my $working_dir = $self->required_param('vep_working');
    make_path($working_dir) unless (-d $working_dir);	    

    my @split_files;
    my @chromosomes = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT);

    for my $chr (@chromosomes) {
	if ($self->param('debug_mode')) {
	    next unless $chr eq 'Y';
	}

	my $output_dir = $working_dir . '/chr' . $chr;
	unless (-d $output_dir) {
	    my $err;
	    make_path($output_dir, {error => \$err});
	    die "make_path failed: ".Dumper($err) if $err && @$err;
	}
    
	my $vcf_file = $self->required_param('ens_var_prefix') . $chr .
	    $self->required_param('ens_var_suffix');
	
	my $nr_lines = $self->required_param('lines_per_input_file');

	my $cmd = "cat $vcf_file | split -d -a 3 -l $nr_lines --numeric-suffixes=1 --verbose - $output_dir/chr${chr}_ | wc -l 1>&2";
	my ($exit_code, $nr_files, $flat_cmd) = $self->run_system_command($cmd);
	die "Splitting $vcf_file failed" unless $exit_code == 0;

	chomp $nr_files;
	for my $n (1 .. $nr_files) {
	    while (length $n < 3) {
		$n = '0' . $n;
	    }
	    
	    my $file = $output_dir .'/chr' . $chr . '_' . $n;
	    push @split_files, $file;
	    
	}
    }
    
    $self->param('vep_output_ids', [map {{vep_input_file => $_}} @split_files]);
}

sub write_output {
    my $self = shift;
    
    $self->dataflow_output_id($self->param('vep_output_ids'), 2);

}

1;
