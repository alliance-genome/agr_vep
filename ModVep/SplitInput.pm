package ModVep::SplitInput;

use strict;

use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub fetch_input {
    my $self = shift;
    
    my $output_dir = $self->required_param('vep_working');
    unless (-d $output_dir) {
	my $err;
	make_path($output_dir, {error => \$err});
	die "make_path failed: ".Dumper($err) if $err && @$err;
    }
    
    my @split_files;

    my $split_file_length = 500000;

    my $header_lines = $self->get_header_lines;
    my $vcf_file = $self->required_param('vcf');
    my $mod = $self->required_param('mod');

    my $cmd = "grep -v '^#' $vcf_file | split -d -a 5 -l 500000 --numeric-suffixes=1 --verbose - $output_dir/${mod}_ | wc -l 1>&2";
    my ($exit_code, $nr_files, $flat_cmd) = $self->run_system_command($cmd);
    die "Splitting $vcf_file failed: $exit_code: $nr_files" unless $exit_code == 0;

    chomp $nr_files;
    for my $n (1 .. $nr_files) {
	while (length $n < 5) {
	    $n = '0' . $n;
	}
	    
	my $file = $output_dir .'/' . $self->required_param('mod') . '_' . $n;
	push @split_files, $file;
	    
	$self->add_header_lines($file, $header_lines, $output_dir)
	    
    }
    
    $self->param('vep_output_ids', [map {{vep_input_file => $_}} @split_files]);
}

sub add_header_lines {
    my ($self, $file, $header_lines, $output_dir) = @_;

    $header_lines =~ s/"/\\"/g;
    my $prepend_cmd = "sed -i '1s|^|${header_lines}|' $file";
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($prepend_cmd);
    die "Prepending headers to $file failed ($exit_code): $flat_cmd - $stderr" unless $exit_code == 0;

    return;
}	    

sub get_header_lines {
    my $self = shift;

    my $header_lines;
    my $file = $self->required_param('vcf');
    open (VCF, "grep '^#' $file |");
    while (<VCF>) {
	chomp;
	$header_lines .= $_ . '\\n';
    }
    close VCF;

    return $header_lines;
}    

sub write_output {
    my $self = shift;
    
    $self->dataflow_output_id($self->param('vep_output_ids'), 2);

}

1;
