package ModVep::ProcessOutput;

use strict;
use warnings;

use File::Copy;
use Path::Class;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');


sub run {
    my $self = shift;

    my $input_file = $self->required_param('vep_input_file') . '.vep.vcf';
    my $output_file = $input_file . '.processed';
    unlink $output_file if -e $output_file;

    my ($nr_csq_cols, $allele_ix, $pick_flag_ix, $consequence_ix);
    my $in_fh = file($input_file)->openr;
    my $out_fh = file($output_file)->openw;
    while (my $line = $in_fh->getline) {
	chomp $line;
	if ($line =~ /^#/) {
	    if ($line =~ /^##INFO=<ID=CSQ.+Format:\s(.+)>$/) {
		my @csq_cols = split(/\|/, $1);
		$nr_csq_cols = scalar(@csq_cols);
		($allele_ix) = grep {$csq_cols[$_] eq 'Allele'} (0 .. $nr_csq_cols - 1);
	        ($pick_flag_ix) = grep {$csq_cols[$_] eq 'PICK'} (0 .. $nr_csq_cols - 1);
	        ($consequence_ix) = grep {$csq_cols[$_] eq 'Consequence'} (0 .. $nr_csq_cols - 1);
		$line =~ s/\|PICK\|/\|Gene_level_consequence\|/;
	    }
	    $out_fh->print("$line\n") if $input_file =~ /_0+1$/;
	    next;
	}
	my @columns = split("\t", $line);
	my ($pre_csq, $csq_string, $post_csq) = $columns[7] =~ /^(.+CSQ=)([^;]+)(.*)$/;
	die "Could not parse line for CSQ:\n$line\n" unless defined $pre_csq and defined $csq_string and defined $post_csq;
	my @csq_entries = split(',', $csq_string);
	my %most_severe_consequences;
	my @csq_splits;
	for my $csq_entry(@csq_entries) {
	    my @csq_cols = split(/\|/, $csq_entry, -1);
	    die "Mismatch in nr. of CSQ fields found - expecting ${nr_csq_cols}, found " . scalar @csq_cols . "\n$line\n"
		unless @csq_cols == $nr_csq_cols;

	    push @csq_splits, \@csq_cols;
	    if ($csq_cols[$pick_flag_ix]) {
		$most_severe_consequences{$csq_cols[$allele_ix]} = $csq_cols[$consequence_ix];
	    }
	}
	my @new_csq_entries;
	for my $csq_entry_split(@csq_splits) {
	    $csq_entry_split->[$pick_flag_ix] = $most_severe_consequences{$csq_entry_split->[$allele_ix]};
	    push @new_csq_entries, join('|', @$csq_entry_split);
	}
	$columns[7] = $pre_csq . join(',', @new_csq_entries) . $post_csq;
	$out_fh->print(join("\t", @columns) . "\n");
    }

    move($output_file, $input_file) or die "Overwriting $input_file with $output_file failed";;
}
    
    
1;
