package ModVep::CombineChromosomes;

use strict;
use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');
use Path::Class;
use File::Path qw(remove_tree);
use Data::Dumper;

sub run {
    my $self = shift;

    my $mod = $self->required_param('mod');
    my $out_file = $self->required_param('out_file_prefix') . 'vcf';
    my $out_fh = file($out_file)->openw;
    my $first_file = 1;
    my $working_dir = dir($self->required_param('vep_working'));
    for my $chromosome_dir ($working_dir->children()) {
	next unless $chromosome_dir->is_dir;
	my $chr_part_nr = 1;
	while (-d "${chromosome_dir}/${chr_part_nr}") {
	    my @files = <${chromosome_dir}/${chr_part_nr}/*>;
	    for my $file (@files) {
		next unless $file =~ /\.vep\.vcf\.processed$/;
		my $in_fh = file($file)->openr;
		while ( my $line = $in_fh->getline()) {
		    next if !$first_file and $line =~ /^#/;
		    chomp $line;
		    $out_fh->print($line . "\n");
		}
		$first_file = 0;
	    }
	    $chr_part_nr++;
	}
	my $err;
	remove_tree($chromosome_dir->stringify, {error => \$err});
	die "remove_tree failed for $chromosome_dir: " .Dumper($err) if $err && @$err;
    }

    if (-e $out_file . ".gz") {
	system("rm ${out_file}.gz");
    }
    if (-e $out_file . ".gz.tbi") {
	system("rm ${out_file}.gz.tbi");
    }
    for my $cmd( "bgzip -l 9 $out_file", "tabix -p vcf ${out_file}.gz") {
	my ($exit_code, $std_err, $flat_cmd) = $self->run_system_command($cmd);
	die "$cmd failed [$exit_code]: $std_err" unless $exit_code == 0;
    }
}

1;
