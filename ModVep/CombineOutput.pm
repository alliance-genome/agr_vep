package ModVep::CombineOutput;

use strict;

use File::Path qw(remove_tree);
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');


sub run {
    my $self = shift;

    my $working_dir = $self->required_param('vep_working');
    my $out_file = $self->required_param('out_file');

    my @files = <$working_dir/*>;
    for my $file (@files) {
	next unless $file =~ /\.vep\.vcf$/;
	my $cmd = "cat $file >> $out_file";
	my ($exit_code, $std_err, $flat_cmd) = $self->run_system_command($cmd);
	die "Rejoining VEP output failed [$exit_code]: $std_err" unless $exit_code == 0;    
    }

    my $err;
    remove_tree($working_dir, {error => \$err});
    die "remove_tree failed for $working_dir: " .Dumper($err) if $err && @$err;
}

1;
