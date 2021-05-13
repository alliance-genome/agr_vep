package ModVep::BaseModVep;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub convert_vcf_chromosomes {
    my ($self, $mod, $vcf_file, $chr_map) = @_;

    my $reverse_map = $self->get_reverse_chromosome_map($chr_map, $mod);

    my $out_file = $vcf_file . '.refseq';
    open (IN, '<', $vcf_file);
    open (OUT, '>', $out_file);
    while (<IN>) {
	if ($_ !~ /^#/) {
	    my @columns = split("\t", $_);
	    if (exists $chr_map->{$mod}{$columns[0]}) {
		$columns[0] = $chr_map->{$mod}{$columns[0]};
		print OUT join("\t", @columns);
		next;
	    }
	    else {
		die "Could not map $mod chromosome " . $columns[0] . " to RefSeq ID\n"
		    unless exists $reverse_map->{$columns[0]};
	    }
	}
	print OUT $_;
    }
    close (IN);
    close (OUT);

    system("mv $out_file $vcf_file");

    return;
}


sub get_reverse_chromosome_map {
    my ($self, $chr_map, $mod) = @_;

    my %reverse_map;
    for my $chr (keys %{$chr_map->{$mod}}) {
	$reverse_map{$chr_map->{$mod}{$chr}} = $chr;
    }

    return \%reverse_map;
}


1;
