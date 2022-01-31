package ModVep::CombineChromosomes;

use strict;
use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');
use Path::Class;
use File::Path qw(remove_tree);
use Data::Dumper;

sub run {
    my $self = shift;

    my $mod = $self->required_param('mod');
    my $out_file = $self->required_param('out_file_prefix') . 'all.vcf';
    my $out_fh = file($out_file)->openw;
    if ($mod eq 'RGD') {
	my $first_file = 1;
	my $working_dir = dir($self->required_param('vep_working'));
	for my $chromosome_dir ($working_dir->children()) {
	    next unless $chromosome_dir->is_dir;
	    for my $subdir ($chromosome_dir->children()) {
		for my $file ($subdir->children()) {
		    next unless $file->stringify() =~ /\.vep\.vcf\.processed$/;
		    my $in_fh = file($file)->openr;
		    while ( my $line = $in_fh->getline()) {
			next if !$first_file and $line =~ /^#/;
			chomp $line;
			$out_fh->print($line . "\n");
		    }
		    $first_file = 0;
		}
	    }
	    my $err;
	    remove_tree($chromosome_dir->stringify, {error => \$err});
	    die "remove_tree failed for $chromosome_dir: " .Dumper($err) if $err && @$err;
	}

	my $cmd = "gzip -9 $out_file";
	my ($exit_code, $std_err, $flat_cmd) = $self->run_system_command($cmd);
	die "$cmd failed [$exit_code]: $std_err" unless $exit_code == 0;
    }
}

1;
