package HumanVep::CombineOutput;

use strict;

use File::Path qw(remove_tree);
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');


sub run {
    my $self = shift;

    my $working_dir = $self->required_param('vep_working');
    my $out_file = $self->required_param('out_file');

    my @chromosomes = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
;

    for my $chr (@chromosomes) {
	my $output_dir = $working_dir . '/chr' . $chr;
	next unless -d $output_dir;

	my @chr_files = <$output_dir/*>;
	for my $file (@chr_files) {
	    next unless $file =~ /\.vep.vcf$/;
	    my $cmd = "cat $file >> $out_file";
	    my ($exit_code, $std_err, $flat_cmd) = $self->run_system_command($cmd);
	    die "Rejoining VEP output failed [$exit_code]: $std_err" unless $exit_code == 0;
	}

	my $err;
	remove_tree($output_dir, {error => \$err});
	die "remove_tree failed for chr $chr: " .Dumper($err) if $err && @$err;
    }
}

1;
