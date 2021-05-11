package ModVep::SplitInput;

use strict;

use File::Path qw(make_path);
use Data::Dumper;
use DBI;

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

    my $cmd = "split -d -a 5 -l $split_file_length --numeric-suffixes=1 --verbose $vcf_file $output_dir/${mod}_ | wc -l 1>&2";
    my ($exit_code, $nr_files, $flat_cmd) = $self->run_system_command($cmd);
    die "Splitting $vcf_file failed: $exit_code: $nr_files" unless $exit_code == 0;

    chomp $nr_files;
    my $folder_nr = 0;
    for my $n (1 .. $nr_files) {
	if ($n % $self->required_param('files_per_folder') == 1) {
	    $folder_nr++;
	    make_path("${output_dir}/${folder_nr}", {error => \$err});
	    die "make_path failed: ".Dumper($err) if $err && @$err;
	}
	while (length $n < 5) {
	    $n = '0' . $n;
	}
	    
	my $file_location = "${output_dir}/" . $self->required_param('mod') . "_$n";
	my $file_destination = "${output_dir}/${folder_nr}/" . $self->required_param('mod') . "_$n";

	$self->run_system_command("mv ${file_location} ${file_destination}");
	push @split_files, $file_destination;
    }
    
    $self->param('vep_output_ids', [map {{vep_input_file => $_}} @split_files]);

}


sub store_transcript_id_map {
    my $self = shift;

    my $dsn = 'dbi:mysql:database=agr_pathogenicity_predictions_' . $self->required_param('mod') . 
	';host=' . $ENV{'WORM_DBHOST'} .';port=' . $ENV{'WORM_DBPORT'};
    my $dbh = DBI->connect($dsn, $ENV{'WORM_DBUSER'}, $self->required_param('password')) or die $DBI::errstr;
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
