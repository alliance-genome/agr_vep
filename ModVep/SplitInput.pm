package ModVep::SplitInput;

use strict;

use File::Path qw(make_path);
use Data::Dumper;
use DBI;
use Path::Class;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub fetch_input {
    my $self = shift;
    
    my $output_dir = $self->required_param('vep_working');
    system("rm -r $output_dir") if -d $output_dir;
    my $err;
    make_path($output_dir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;


    # Store transcript ID to name map for later use by TranscriptName VEP plugin
    $self->store_transcript_id_map();
    
    my @split_files;
    my $split_file_length = $self->required_param('lines_per_file');
    
    my $vcf_file = $self->required_param('vcf');
    my $mod = $self->required_param('mod');

    my $in_fh = file($vcf_file)->openr or die "Could not open $vcf_file for reading\n";
    my $out_fh = file($output_dir . '/headers.vcf')->openw or die "Could not open ${output_dir}/headers.vcf for writing\n";
    my $line_nr = 0;
    my $chr_nr = 0;
    my $chr_file_nr = 0;
    my $chr_folder_nr = 0;
    my $previous_chr = 'none';
    my $chr_file_suffix;
    while (my $line = $in_fh->getline()) {
	if ($line =~ /^#/) {
	    $out_fh->print($line);
	}
	else {
	    $line_nr++;
	    my ($chr) = $line =~ /^(\S+)\s/;
	    if ($chr ne $previous_chr) {
		$previous_chr = $chr;
		$chr_nr++;
		$chr_file_nr = 0;
		$chr_folder_nr = 0;
		$line_nr = 1;
		make_path("${output_dir}/${chr_nr}", {error => \$err});
		die "make_path failed: ".Dumper($err) if $err && @$err;
	    }
	    if ($line_nr % $self->required_param('lines_per_file') == 1) {
		$chr_file_nr++;
		if ($chr_file_nr % $self->required_param('files_per_folder') == 1) {
		    $chr_folder_nr++;
		    make_path("${output_dir}/${chr_nr}/${chr_folder_nr}", {error => \$err});
		    die "make_path failed: ".Dumper($err) if $err && @$err;
		}
		$chr_file_suffix = $chr_file_nr;
		while (length $chr_file_suffix < 5) {
		    $chr_file_suffix = '0' . $chr_file_suffix;
		}
		my $filename = "${output_dir}/${chr_nr}/${chr_folder_nr}/${mod}_${chr_file_suffix}";
		$out_fh = file($filename)->openw or die "Could not open $filename for writing\n";
		$self->print_header($out_fh) if $chr_file_nr == 1;
		push @split_files, $filename;
	    }
	    $out_fh->print($line);
	}
    }
    
    $self->param('vep_output_ids', [map {{vep_input_file => $_}} @split_files]);

}

sub print_header {
    my ($self, $out_fh) = @_;

    my $in_fh = file($self->required_param('vep_working') . '/headers.vcf')->openr or die "Could not open headers file for reading\n";
    while (my $line = $in_fh->getline()) {
	$out_fh->print($line);
    }

    return;
}


sub store_transcript_id_map {
    my $self = shift;

    my $dsn = 'dbi:mysql:database=agr_pathogenicity_predictions_' . $self->required_param('mod') . 
	';host=' . $ENV{'VEP_DBHOST'} .';port=' . $ENV{'VEP_DBPORT'};
    my $dbh = DBI->connect($dsn, $ENV{'VEP_DBUSER'}, $self->required_param('password')) or die $DBI::errstr;
    my @create_sql = (
	"DROP TABLE IF EXISTS transcript_map;",
	"CREATE TABLE transcript_map (transcript_id varchar(50), transcript_name varchar(50));"
	);
    foreach my $create_sql (@create_sql) {
	$dbh->do($create_sql) or $self->throw("Failed to execute: $create_sql");
    }
  
    my $insert_sql = "INSERT INTO transcript_map VALUES(?, ?)";
    my $sth = $dbh->prepare($insert_sql);

    my %map;
    my $gff_file = $self->required_param('gff');
    open (GFF, "gunzip -c $gff_file|");
    while (<GFF>) {
	next if $_ =~ /^#/;
	my @columns = split("\t", $_);
	next if $columns[2] eq 'exon';
	my %attributes = split(/[=;]/, $columns[8]);
	my $transcript_id = $attributes{'transcript_id'};
	next unless $transcript_id;
	$transcript_id =~ s/\s+$//;
	my $transcript_name = $attributes{'Name'} || $attributes{'name'} || $transcript_id;
	$transcript_name =~ s/\s+$//;
	$sth->execute($transcript_id, $transcript_name);
    }
    close (GFF);

    $dbh->do("CREATE INDEX transcript_idx ON transcript_map (transcript_id)");

    return;
}


sub write_output {
    my $self = shift;
    
    $self->dataflow_output_id($self->param('vep_output_ids'), 2);

}


1;
