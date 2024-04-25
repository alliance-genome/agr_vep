package ModVep::CombineOutput;

use strict;

use File::Path qw(remove_tree);
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');


sub run {
    my $self = shift;

    my $working_dir = $self->required_param('vep_working');
    my $chromosome_nr = $self->required_param('chromosome_nr');
    my $mod = $self->required_param('mod');
    my $chromosome_id = `tail -n 10 ${working_dir}/${chromosome_nr}/1/${mod}_00001.vep.vcf.processed | head -n 1 | cut -f 1`;
    chomp($chromosome_id);
    my $out_file = $self->required_param('out_file_prefix') . $chromosome_id . '.vcf';
    
    my $folder_nr = 1;
    my ($exit_code, $std_err, $flat_cmd);
    while (-d "${working_dir}/${chromosome_nr}/${folder_nr}") {
	my @files = <${working_dir}/${chromosome_nr}/${folder_nr}/*>;
	for my $file (@files) {
	    next unless $file =~ /\.vep\.vcf\.processed$/;
	    my $cmd = "cat $file >> $out_file";
	    ($exit_code, $std_err, $flat_cmd) = $self->run_system_command($cmd);
	    die "Rejoining VEP output failed [$exit_code]: $std_err" unless $exit_code == 0;    
	}
	$folder_nr++;
    }

    if (-e $out_file . ".gz") {
	system("rm ${out_file}.gz");
    }
    if (-e $out_file . ".gz.tbi") {
	system("rm ${out_file}.gz.tbi");
    }
    for my $cmd("bgzip -l 9 ${out_file}", "tabix -p vcf ${out_file}.gz", "vcf-validator ${out_file}.gz") {
	($exit_code, $std_err, $flat_cmd) = $self->run_system_command($cmd);
	die "$cmd failed [$exit_code]: $std_err" unless $exit_code == 0;
    }

}



1;
