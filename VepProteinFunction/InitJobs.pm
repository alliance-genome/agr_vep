package VepProteinFunction::InitJobs;

use strict;

use VepProteinFunction::Constants qw(FULL UPDATE NONE);
use Digest::MD5 qw(md5_hex);
use Bio::DB::Fasta;
use Bio::Seq;
use File::Path qw(make_path);

use base qw(VepProteinFunction::BaseVepProteinFunction);

sub fetch_input {
    my $self = shift;
    
    my $sift_run_type = $self->required_param('sift_run_type');
    my $pph_run_type = $self->required_param('pph_run_type');
    
    foreach my $type (qw/sift pph/) {
	if ($self->required_param("$type\_run_type") != NONE) {
	    my $working_dir = $self->required_param("$type\_working");
	    make_path($working_dir) unless (-d $working_dir);	    
	}
    }  

    my $pf_dsn = 'dbi:mysql:database=' . $self->required_param('pfdb_name') .
	';host=' . $self->required_param('pfdb_host') .
	';port=' . $self->required_param('pfdb_port');
    my $pf_dbh = DBI->connect($pf_dsn, $self->required_param('pfdb_user'),
			      $self->required_param('pfdb_pass')) or die $DBI::errstr;

    for my $table ('protein_function_prediction', 'valid_failure', 'analysis_attempt'){
	my $remove_existing_results_query = 
	    "DELETE t.* FROM " . $table . " t, analysis a " .
	    "WHERE t.analysis_id = a.analysis_id AND a.analysis = ?";
    
	my $remove_existing_results_sth = $pf_dbh->prepare($remove_existing_results_query);
	$remove_existing_results_sth->execute('sift') if $sift_run_type == FULL;
	$remove_existing_results_sth->execute('pph') if $pph_run_type == FULL;
    }
    
    my ($translation_seqs, $transcript_id_md5_map, $translation_md5_codontable_ids) =
	$self->get_agr_translation_seqs();
    
    $self->update_pfdb_transcript_md5s($pf_dbh, $transcript_id_md5_map);
     
    my $pep_fasta = $self->required_param('pep_fasta');
    my @dir = split('/', $pep_fasta);
    pop @dir;
    make_path(join('/', @dir));
    
    open my $PEP_FASTA, ">$pep_fasta" or die "Failed to open $pep_fasta for writing\n";
	
    # get rid of any existing index file
    if (-e "$pep_fasta.fai") {
	unlink "$pep_fasta.fai" or die "Failed to delete $pep_fasta index file\n";
    }
	
    for my $md5 (keys %$translation_seqs) {
	my $seq = $translation_seqs->{$md5};
	$seq =~ s/(.{80})/$1\n/g;
	chomp $seq;
	print $PEP_FASTA ">$md5\n$seq\n";
    }
    close $PEP_FASTA;
	
    $self->run_cmd("samtools faidx $pep_fasta");

    my %existing_translation_md5s;
    if ($sift_run_type == UPDATE || $pph_run_type == UPDATE){
	my $existing_md5s_query = qq{
            SELECT a.analysis, t.translation_md5
                FROM protein_function_prediction p
                INNER JOIN translation_md5 t
                    ON p.translation_md5_id = t.translation_md5_id
                INNER JOIN analysis a
                    ON p.analysis_id = a.analysis_id
            UNION
            SELECT a.analysis, t.translation_md5
                FROM valid_failure v
                INNER JOIN translation_md5 t
                    ON v.translation_md5_id = t.translation_md5_id
                INNER JOIN analysis a
                    ON v.analysis_id = a.analysis_id
            UNION
            SELECT a.analysis, t.translation_md5
                FROM analysis_attempt at
                INNER JOIN translation_md5 t
                    ON at.translation_md5_id = t.translation_md5_id
                INNER JOIN analysis a
                    ON at.analysis_id = a.analysis_id
                WHERE at.attempt >= ?
        };
    
	my $existing_sth = $pf_dbh->prepare($existing_md5s_query);
	$existing_sth->execute($self->required_param('max_attempts'));
	while (my $existing_result = $existing_sth->fetchrow_arrayref){
	    $existing_translation_md5s{$existing_result->[0]}{$existing_result->[1]} = 1;
	}
    }
    
    my @sift_md5s;
    if ($sift_run_type == FULL) {
	@sift_md5s = keys %$translation_seqs;
    }
    elsif ($sift_run_type == UPDATE) {
	for my $md5 (keys %$translation_seqs){
	    push @sift_md5s, $md5 unless exists $existing_translation_md5s{'sift'}{$md5};
	}
    }
    
    my @pph_md5s;
    if ($pph_run_type == FULL) {
	@pph_md5s = keys %$translation_seqs;
    }
    elsif ($pph_run_type == UPDATE) {
	for my $md5 (keys %$translation_seqs){
	    push @pph_md5s, $md5 unless exists $existing_translation_md5s{'pph'}{$md5};
	}
    }

    $self->update_pfdb_attempts($pf_dbh, 'sift', \@sift_md5s);
    $self->update_pfdb_attempts($pf_dbh, 'pph', \@pph_md5s);
    
    my %required_md5s = map {$_ => 1} (@sift_md5s, @pph_md5s);
        
    if ($self->param('debug_mode')) {
	my $debug_nr = 100;
	@sift_md5s = @sift_md5s[0 .. $debug_nr - 1];
	@pph_md5s = @pph_md5s[0 .. $debug_nr - 1];
    }

    $self->param('sift_output_ids', [map { {translation_md5 => $_, codontable_id => $translation_md5_codontable_ids->{$_}} } @sift_md5s]);
    $self->param('pph_output_ids', [map { {translation_md5 => $_, codontable_id => $translation_md5_codontable_ids->{$_}} } @pph_md5s]);
}


sub cigar_exon_bases {
    my $cigar = shift;
    
    my $exon_bases = 0;
    while (length $cigar > 0) {
	my ($bases, $operator) = $cigar =~ /^(\d+)(\D)/;
	$cigar = substr($cigar, length($bases) + 1);
	if ($operator =~ /^[MX=D]$/) {
	    $exon_bases += $bases;
	    
	}
    }

    return $exon_bases;
}	


sub construct_cds_from_fasta {
    my ($self, $transcript, $details, $fasta_db, $chromosome_map) = @_;

    my $cds_seq = '';
    my @cds_starts = sort {$a<=>$b} keys %{$details->{$transcript}{'CDS'}};
    for my $cds_start (@cds_starts) {
	my $chromosome = exists $chromosome_map->{$details->{$transcript}{'chromosome'}} ?
	    $chromosome_map->{$details->{$transcript}{'chromosome'}} : $details->{$transcript}{'chromosome'};
	$cds_seq .= $fasta_db->seq($chromosome, $cds_start, $details->{$transcript}{'CDS'}{$cds_start}{'end'});
    }

    $cds_seq = trim_incomplete_codons($cds_seq, $transcript, $details, \@cds_starts);

    return $cds_seq;
}


sub get_agr_translation_seqs {
    my $self = shift;
    
    my $agr_fasta = $self->required_param('agr_fasta');
    if (-e "$agr_fasta.index") {
	unlink "$agr_fasta.index" or die "Failed to delete $agr_fasta.index\n";
    }
    my $agr_fasta_db = Bio::DB::Fasta->new($agr_fasta);

 
    my $alignments = $self->get_alignments();
    my ($transcript_details, $chromosome_map) = $self->parse_gff($alignments);


    my (%translation_seqs, %transcript_id_translation_md5_map,
	%translation_md5_codontable_ids);
    

    for my $transcript (keys %$transcript_details) {

	# GFF parse will have pulled out details for non-coding transcripts too, so skip them
	next unless exists $transcript_details->{$transcript}{'CDS'};
	next unless exists $transcript_details->{$transcript}{'transcript_id'};
	
	my $cds_seq;
	if (exists $alignments->{$transcript_details->{$transcript}{'transcript_id'}}) {
	    $cds_seq = 
		$self->get_cds_from_alignment($transcript, $transcript_details, $alignments);
	}
	else {
	    $cds_seq = $self->construct_cds_from_fasta($transcript, $transcript_details,
						       $agr_fasta_db,
						       $chromosome_map);
	}
	next if !defined $cds_seq;
	
	my $translated_seq = Bio::Seq->new(
	    -seq => $cds_seq,
	    -id => 'transcript',
	    -display_id => 'transcript',
	    -alphabet => 'dna',
	    );
	$translated_seq = $translated_seq->revcom if $transcript_details->{$transcript}{'strand'} eq '-';
	my $translation = $translated_seq->translate(
	    -codontable_id => $transcript_details->{$transcript}{'codontable_id'},
	    )->seq;
	$translation =~ s/\*.*$//;
	my $translation_md5 = md5_hex($translation);

	$translation_seqs{$translation_md5} = $translation;
	$transcript_id_translation_md5_map{$transcript} = $translation_md5;
	$translation_md5_codontable_ids{$translation_md5} =
	    $transcript_details->{$transcript}{'codontable_id'};
    }

    return (\%translation_seqs, \%transcript_id_translation_md5_map,
	\%translation_md5_codontable_ids);
}


sub get_alignments {
    my $self = shift;

    my %alignments;
    my %aligned_transcripts;
    
    my $bam = $self->required_param('agr_bam');
    return \%alignments unless $bam;

    open (BAM, "samtools view $bam |");
    while (<BAM>) {
	my @columns = split("\t", $_);
	my $transcript = $columns[0];

	next if $self->required_param('mod') eq 'human' and $columns[2] !~ /^NC_/; # only want alignments against assembled chromosomes
	if (exists $aligned_transcripts{$transcript}) {
	    warn "Multiple alignments of $transcript to assembled chromosomes - 1st used\n"
		unless $self->param('mod') eq 'human' and (
		    ($columns[2] =~ /23\./ and $aligned_transcripts{$transcript} =~ /24\./)
		    or ($columns[2] =~ /24\./ and $aligned_transcripts{$transcript} =~ /23\./));
	    next;
	}
	$aligned_transcripts{$transcript} = $columns[2];

	my $strand;
	if ($columns[1] == 0) {
	    $strand = '+';
	}
	else {
	    die "Unexpected bitcode" unless $columns[1] == 16;
	    $strand = '-';
	}

	$alignments{$transcript}{'start'} = $columns[3];
	$alignments{$transcript}{'cigar'} = $columns[5];
	$alignments{$transcript}{'strand'} = $strand;
	$alignments{$transcript}{'seq'} = $columns[9];
    }
    close (BAM);

    return \%alignments;
}


sub get_attributes {
    my $attributes_string = shift;

    my %attributes;
    my @key_value_pairs = split(';', $attributes_string);
    for my $kvp (@key_value_pairs){
	next unless $kvp =~ /=/;
	my ($key, $value) = $kvp =~ /^(.+)=(.*)$/;
	$attributes{$key} = $value;
    }

    return \%attributes;
}


sub get_cds_from_alignment {
    my ($self, $transcript, $details, $alignment) = @_;

    my @cds_starts = sort {$a<=>$b} keys %{$details->{$transcript}{'CDS'}};
    my @exon_starts = sort {$a<=>$b} keys %{$details->{$transcript}{'exon'}};

    my $cds_start = $cds_starts[0];
    my $cds_end = $details->{$transcript}{'CDS'}{$cds_starts[-1]}{'end'};

    my $exon_bases_before_cds = 0;
    my $exon_bases_after_cds = 0;

    my $gff_exon_bases = 0;
    for my $exon_start (@exon_starts) {
	my $exon_end = $details->{$transcript}{'exon'}{$exon_start}{'end'};
	$gff_exon_bases += ($exon_end - $exon_start) + 1;
	if ($exon_start < $cds_start) {
	    if ($exon_end < $cds_start) {
		$exon_bases_before_cds += ($exon_end - $exon_start) + 1;
	    }
	    else{
		$exon_bases_before_cds += $cds_start - $exon_start;
	    }
	}
	if ($exon_end > $cds_end) {
	    if ($exon_start > $cds_end) {
		$exon_bases_after_cds += ($exon_end - $exon_start) + 1;
	    }
	    else {
		$exon_bases_after_cds += ($exon_end - $cds_end);
	    }
	}
    }

    my $alignment_seq = $alignment->{$details->{$transcript}{'transcript_id'}}{'seq'};
    my $cigar = $alignment->{$details->{$transcript}{'transcript_id'}}{'cigar'};
    my $cigar_exon_bases = cigar_exon_bases($cigar);
    if ($cigar_exon_bases != $gff_exon_bases) {
	warn "BAM/GFF mismatch for " . $details->{$transcript}{'transcript_id'} . " ${cigar_exon_bases}:${gff_exon_bases}\n";
	return;
    }

    my ($bases, $operator);
    my $bases_to_trim_before = $exon_bases_before_cds;
    my $possible_three_prime_coding_insertion = 0;
    if ($exon_bases_before_cds > 0) {
	my $bases_counted = 0;
	while (length $cigar > 0) {
	    ($bases, $operator) = $cigar =~ /^(\d+)(\D)/;
	    $cigar = substr($cigar, length($bases) + 1);
	    if ($operator =~ /^[MX=]$/) {
		$bases_counted += $bases;
		last if $bases_counted >= $exon_bases_before_cds;
	    }
	    elsif ($operator =~/^[DN]$/) {
		$bases_to_trim_before -= $bases if $operator eq 'D';
		$bases_counted += $bases;
		last if $bases_counted >= $exon_bases_before_cds;
	    }
	    elsif ($operator eq 'I') {
		$bases_to_trim_before += $bases;
	    }
	}
    }
    else {
	if ($cigar =~ /^\d+I/) {
	    my ($inserted_bases) = $cigar =~ /^\d+I/;
	    if ($alignment->{$details->{$transcript}{'transcript_id'}}{'strand'} eq '+') {
		my $start_codon_ix = rindex(substr($alignment_seq, 0, $inserted_bases + 2), 'ATG');
		$bases_to_trim_before = $start_codon_ix unless $start_codon_ix == -1;
	    }
	}
    }
       
    $cigar = $alignment->{$details->{$transcript}{'transcript_id'}}{'cigar'};
    my $bases_to_trim_after = $exon_bases_after_cds;
    if ($exon_bases_after_cds > 0) {
	my $bases_counted = 0;
	while (length $cigar > 0) {
	    ($bases, $operator) = $cigar =~ /\D(\d+)(\D)$/;
	    if (!defined $bases) {
		($bases, $operator) = $cigar =~ /^(\d+)(\D)$/;
	    }
	    $cigar = substr($cigar, 0, -(length($bases) + 1));
	    if ($operator =~ /^[MX=]$/) {
		$bases_counted += $bases;
		last if $bases_counted >= $exon_bases_after_cds;
	    }
	    elsif ($operator =~ /DN/) {
		$bases_to_trim_after -= $bases if $operator eq 'D';
		$bases_counted += $bases;
		last if $bases_counted >= $exon_bases_after_cds;
	    }
	    elsif ($operator eq 'I') {
		$bases_to_trim_after += $bases;
	    }
	}
    }
    else {
	if ($cigar =~ /I$/) {
	    my ($inserted_bases) = $cigar =~ /\D*(\d+)I$/;
	    if ($alignment->{$details->{$transcript}{'transcript_id'}}{'strand'} eq '-') {
		my $start_codon_ix = index(substr($alignment_seq, -($inserted_bases+2)), 'CAT');
		$bases_to_trim_after = length($alignment_seq) - ($start_codon_ix + 3) unless $start_codon_ix == -1;
	    }
	}
    }

    my $cds_seq = substr($alignment_seq, $bases_to_trim_before);
    if ($exon_bases_after_cds > 0){
	$cds_seq = substr($cds_seq, 0, -$bases_to_trim_after);
    }

    return $cds_seq;
}


sub parse_gff {
    my $self = shift;

    my (%details, %chromosome_map);

    open (GFF, '< ' . $self->required_param('agr_gff')) or die "Couldn't open GFF for reading\n";
    while (<GFF>) {
	next if $_ =~ /^#/;
	chomp;
	my @columns = split("\t", $_);

	if ($self->param('mod') eq 'human') {
	    # only want transcripts on assembled chromosomes
	    next unless $columns[0] =~ /^NC_/ or $columns[0] =~ /^\d\d?$/ or $columns[0] =~ /^[XY]$/;
	}

	my $attributes = get_attributes($columns[8]);
	if ($columns[2] eq 'mRNA') {
	    my $id = $attributes->{'ID'} || $attributes->{'id'};
	    my $transcript_id = $attributes->{'Name'} || $attributes->{'name'} ||
		$attributes->{'transcript_id'} || $id;
	    $details{$id}{'transcript_id'} = $transcript_id;
	    
	    $details{$id}{'chromosome'} = $columns[0];
	    $details{$id}{'strand'} = $columns[6];
	    
	    my $codontable_id = 1;
	    if ($columns[0] =~ /^M[tT]/ or $columns[0] =~ /^mito/i or
		$columns[0] =~ /^chrM[tT]/) {
		# Set codon table to invertebrate mitochondrial for FB and WB, 
		# vertebrate mitochondrial for other MODs
		$codontable_id = $self->required_param('mod') =~ /^[FW]B$/ ? 5 : 2;	  
	    }
	    $details{$id}{'codontable_id'} = $codontable_id;
	}
	elsif ($columns[2] eq 'exon' or $columns[2] eq 'CDS') {
	    my $start = $columns[3];
	    for my $parent (split(',', $attributes->{'Parent'})) {
		$details{$parent}{$columns[2]}{$start}{'end'} = $columns[4];
		$details{$parent}{$columns[2]}{$start}{'phase'} = $columns[7];
	    }
	}
	elsif ($columns[2] eq 'region') {
	    $chromosome_map{$columns[0]} = $attributes->{'chromosome'};
	}
	else {
	    # Other entries can be ignored
	}
    }
    close (GFF);
	    
    return (\%details, \%chromosome_map);
}


sub trim_incomplete_codons {
    my ($cds_seq, $transcript, $details, $cds_starts) = @_;

    if ($details->{$transcript}{'strand'} eq '+') {
	my $first_cds_phase = $details->{$transcript}{'CDS'}{$cds_starts->[0]}{'phase'};
	if ($first_cds_phase != 0) {
	    my $bases_to_trim = 3 - $first_cds_phase;
	    $cds_seq = substr($cds_seq, $bases_to_trim);
	}
	my $incomplete_end_codon_bases = (length $cds_seq) % 3;
	$cds_seq = substr($cds_seq, 0, -$incomplete_end_codon_bases) if $incomplete_end_codon_bases != 0;
    }
    else {
	my $first_cds_length = ($details->{$transcript}{'CDS'}{$cds_starts->[-1]}{'end'} - $cds_starts->[-1]) + 1;
	my $first_cds_phase = $details->{$transcript}{'CDS'}{$cds_starts->[-1]}{'phase'};
	if ($first_cds_phase != 0) {
	    my $bases_to_trim = 3 - $first_cds_phase;
	    $cds_seq = substr($cds_seq, 0, -$bases_to_trim);
	}
	my $incomplete_end_codon_bases = (length $cds_seq) % 3;
	$cds_seq = substr($cds_seq, $incomplete_end_codon_bases);
    }
	
    return $cds_seq;
}


sub update_pfdb_attempts {
    my ($self, $pf_dbh, $analysis, $md5s) = @_;

    my $update_attempts_query = qq{
        INSERT INTO analysis_attempt (translation_md5_id, analysis_id, attempt)
	    VALUES (
                (SELECT translation_md5_id FROM translation_md5
                    WHERE translation_md5 = ?),
                (SELECT analysis_id FROM analysis
                    WHERE analysis = ?),
                1
            )
        ON DUPLICATE KEY
            UPDATE attempt = attempt + 1
    };
    my $update_attempts_sth = $pf_dbh->prepare($update_attempts_query);
    for my $md5 (@$md5s) {
	$update_attempts_sth->execute($md5, $analysis);
    }

    return;
}


sub update_pfdb_transcript_md5s {
    my ($self, $pf_dbh, $transcript_id_md5_map) = @_;

    my $update_translation_md5_query = qq{
         INSERT IGNORE INTO translation_md5(translation_md5)
             VALUES(?)
    };
    my $update_transcript_query = qq{
         INSERT IGNORE INTO transcript(transcript_id, translation_md5_id)
             VALUES(
                 ?,
                 (SELECT translation_md5_id FROM translation_md5
                     WHERE translation_md5 = ?)
         )
    };
    my $update_translation_md5_sth = $pf_dbh->prepare($update_translation_md5_query);
    my $update_transcript_sth = $pf_dbh->prepare($update_transcript_query);
    for my $transcript_id (keys %$transcript_id_md5_map){
	$update_translation_md5_sth->execute($transcript_id_md5_map->{$transcript_id});
	$update_transcript_sth->execute($transcript_id, $transcript_id_md5_map->{$transcript_id});
    }

    return;
}

sub write_output {
    my $self = shift;
    
    unless ($self->param('sift_run_type') == NONE) {
        $self->dataflow_output_id($self->param('sift_output_ids'), 3);
    }
    unless ($self->param('pph_run_type') == NONE) {
	$self->dataflow_output_id($self->param('pph_output_ids'), 2);
    }

}

1;
