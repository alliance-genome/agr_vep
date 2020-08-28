#!/usr/bin/env perl
use strict;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

update_polyphen_uniprot.pl

=head1 SYNOPSIS

update_polyphen_uniprot.pl [options]

Where options are:

 -n <organism>    specify common organism name instead of the default human;
                  supported names are: chimp, orangutan, macaque (crab-eating),
                  mouse, rat, dog, cow, zebrafish and fruitfly; some taxonomic
                  groups of species are also supported, based on UniProtKB set
                  of taxonomic divisions. For a full list, please see:
                    ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/
 -s <dat>[,<idm>] specify local source file(s): <dat> with the protein entries in
                  UniProtKB plain text format and optional <idm> with idmapping data;
                  default is to fetch all source files from ftp.uniprot.org
 -k               keep source files; default is to remove source files after
                  they were successfully processed
 -h               print this help

=head1 DESCRIPTION

Downloads Swiss-Prot, TrEMBL, ID mapping, and VARSPLIC data from ftp.uniprot.org
and prepares protein sequence database and other auxillary database files
for use with the PolyPhen-2.  This script has been modified from the original 
script supplied with PolyPhen-2 due to changes in the UniProt download format.

=head1 AUTHOR

 Vasily Ramensky
 Steffen Schmidt
 Ivan Adzhubey
 Mark Quinton-Tulloch

=head1 SUBVERSION

 $LastChangedDate: 2020-08-27 10:41:00 $
 $LastChangedBy: Mark Quinton-Tulloch $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use DBI;
use LWP::Simple;
use Getopt::Std;
use Pod::Usage;
use File::Basename;
use Storable qw(nstore);

# We parse total 26 features of interest
my @fttokens = qw(
  ACT_SITE BINDING CARBOHYD CA_BIND COILED COMPBIAS CROSSLNK DISULFID DNA_BIND
  INIT_MET INTRAMEM LIPID METAL MOD_RES MOTIF NON_STD NP_BIND PEPTIDE PROPEP
  REGION REPEAT SIGNAL SITE TRANSIT TRANSMEM ZN_FING
);

my $UNIPROT_FTP  = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase';
#-----------------------------------------------------------------------------

my %opts;
getopts('hkln:s:', \%opts) or pod2usage();
exists $opts{'h'} and pod2usage({ -verbose=>99 });

# Organism common name: human (default), orangutan, macaque (crab-eating), mouse, rat, dog, cow, zebrafish
# For UniProt organism codes (spname), see: http://www.uniprot.org/docs/speclist
my $CNAME = 'human';

# SPECIES hash indexed by CNAME
my %SPECIES = (
    # Human
    'human'    => {
	'os'     => 'Homo sapiens',
	'spname' => 'HUMAN',
	'taxn'   => 9606,
	'file'   => 'human',
    },
    # Chimpanzee
    'chimp'    => {
	'os'     => 'Pan troglodytes',
	'spname' => 'PANTR',
	'file'   => 'mammals',
    },
    # Orangutan
    'orangutan'    => {
	'os'     => 'Pongo pygmaeus|Pongo abelii',
	'spname' => 'PONPY|PONAB',
	'file'   => 'mammals',
    },
    # Crab-eating macaque
    'macaque'    => {
	'os'     => 'Macaca fascicularis',
	'spname' => 'MACFA',
	'file'   => 'mammals',
    },
    # Mouse
    'mouse'    => {
	'os'     => 'Mus musculus',
	'spname' => 'MOUSE',
	'taxn'   => 10090,
	'file'   => 'rodents',
    },
    # Rat
    'rat'      => {
	'os'     => 'Rattus norvegicus|Rattus rattus',
	'spname' => 'RAT|RATRT',
	'file'   => 'rodents',
    },
    # Dog
    'dog'      => {
	'os'     => 'Canis lupus familiaris',
	'spname' => 'CANLF',
	'file'   => 'mammals',
    },
    # Cow
    'cow'      => {
	'os'     => 'Bos taurus',
	'spname' => 'BOVIN',
	'file'   => 'mammals',
    },
    # Zebrafish
    'zebrafish' => {
	'os'      => 'Danio rerio',
	'spname'  => 'DANRE',
	'taxn'   => 7955,
	'file'    => 'vertebrates',
    },
    # Fruit fly
    'fruitfly' => {
	'os'      => 'Drosophila melanogaster',
	'spname'  => 'DROME',
	'taxn'   => 7227,
	'file'    => 'invertebrates',
    },
  # Roundworm
    'roundworm' => {
	'os' => 'Caenorhabditis elegans',
	'spname' => 'CAEEL',
	'taxn' => 6239,
	'file' => 'invertebrates',
    },
  # Yeast
    'yeast' => {
	'os' => 'Saccharomyces cerevisiae',
	'spname' => 'YEAST',
	'taxn' => 559292,
	'file' => 'fungi',
    },
  
  #
  # Taxonomic divisions
  #
  # bacteria
  'bacteria' => {
    'os'      => '',
    'spname'  => '',
    'file'    => 'bacteria',
  },
);

if (my $cname = $opts{'n'}) {
  die "Illegal or unsupported organism name: $cname\n"
    unless exists $SPECIES{$cname};
  $CNAME = $cname;
  warn "## Do not forget to make the following\n## changes to your options.cnf file\n## (see README for details):\n";
  warn "REFORGANISM  = $SPECIES{$CNAME}{'os'}\n";
  warn "REFORGCOMMON = $CNAME\n";
  warn "MAPGENE      = 0\n";
}

my $AA = '[ACDEFGHIKLMNPQRSTVWY]';

my $seqfile     = "$CNAME.seq";
my $ftfile      = "$CNAME.ft";
my $swfile      = "$CNAME.dat";
my $id2accfile  = "$CNAME.id2acc";
my $snp2accfile = "$CNAME.snp2acc";
my $xdb2accfile = "$CNAME.xdb2acc";
my $dbfile      = "$CNAME.sqlite";

#-----------------------------------------------------------------------------

my @names = (
  'uniprot_sprot_'  . $SPECIES{$CNAME}{'file'} . '.dat',
  'uniprot_trembl_' . $SPECIES{$CNAME}{'file'} . '.dat'
);

if (exists $SPECIES{$CNAME}{'taxn'} && length $SPECIES{$CNAME}{'taxn'}) {
  push @names, $SPECIES{$CNAME}{'spname'} . '_' . $SPECIES{$CNAME}{'taxn'} . '_idmapping_selected.tab';
} else {
  push @names, 'idmapping_selected.tab';
}

$| = 1;

# Load precompiled data into SQLite database and quit
createdb($dbfile), exit(0) if $opts{'l'};

if (defined $opts{'s'}) {

  my @files = split /,/, $opts{'s'};
  foreach my $file (@files) {
    foreach my $ext (qw{ .seq .ft .dat .id2acc .snp2acc .xdb2acc}) {
      die "ERROR: Can't simultaneously read from and write to the same file, please rename input file: $file\n"
        if $file eq "$CNAME$ext";
    }
  }
  @names = @files;

} else {

  my $suf = '.gz';

  foreach my $name (@names) {
    my $infile;
    if ($name =~ /idmapping_selected/) {
      if (exists $SPECIES{$CNAME}{'taxn'} && length $SPECIES{$CNAME}{'taxn'}) {
        $infile = "$UNIPROT_FTP/idmapping/by_organism/${name}${suf}";
      } else {
        $infile = "$UNIPROT_FTP/idmapping/${name}${suf}";
      }
    } else {
      $infile = "$UNIPROT_FTP/taxonomic_divisions/${name}${suf}";
    }
    my ($outfile) = fileparse($infile);
    unlink $outfile if -e $outfile; # otherwise mirror() will fail
    print "Fetching $infile ...\n";
    my $rc = mirror($infile, $outfile);
    die "Error fetching URL $infile: ".status_message($rc)."\n" unless is_success($rc);
    print "Decompressing $outfile ...\n";
    `gunzip -f $outfile`;
    print "Done.\n";
  }

}

open(SEQ,     ">$seqfile")      || die "Can't create file: $seqfile";
open(FT,      ">$ftfile" )      || die "Can't create file: $ftfile";
open(SWALL,   ">$swfile" )      || die "Can't create file: $swfile";
open(ID2ACC,  ">$id2accfile" )  || die "Can't create file: $id2accfile";
open(SNP2ACC, ">$snp2accfile" ) || die "Can't create file: $snp2accfile";

my $dbase_len = 0;

my $total_filtered = 0;
foreach my $file (@names) {

  next if $file =~ /idmapping_selected/;

  print "Processing $file ...\n";

  my $num_filtered = 0;
  my $copy = 0;
  my (@entrylines, @snps, @ftlines, @seqlines, @accs, @gnames);
  my ($id, $seq, $de, $len);

  open(TMP, $file) || die "Can't open $file\n";

  LINE: while (<TMP>) {

    if (/^ID.*\b(\d+) AA/) {
      @seqlines = @snps = @ftlines = @entrylines = @accs = @gnames = ();
      $copy = 0;
      $id = (split)[1];
      $len = $1;
      undef $de;
    } elsif( /^OS   / ) {
      # there exist multiple OS-lines
      $copy++ if !$SPECIES{$CNAME}{'os'} || m/^OS   (?:$SPECIES{$CNAME}{'os'})(?:\s+\(.+\))?\.*$/o;
    } # end if
    elsif ( /^AC / ) {
      my @accs_tmp = split;
      shift @accs_tmp; # remove ^AC
      map { s/;$// } @accs_tmp;
      push @accs, @accs_tmp;
    } # end elsif
    elsif ( /^GN / ) {
	my $gnline = $_;
	chomp $gnline;
	while (<TMP> ){
	    last unless $_ =~ /^GN /;
	    chomp($gnline .= substr($_, 4));
	}
      push @gnames, $1 if $gnline =~ /Name=([^;]+)/;
      if ($gnline =~ /Synonyms=([^;]+)/) {
        my @n = split /,\s*/, $1;
        push @gnames, @n;
      }
    }
    elsif ( /^FT\s+(\w+)\s/ ) {
	push @entrylines, $_;
	my $feature = $1;
	my $ftline  = $_;
	chomp $ftline;
	while (<TMP>) {
	    my $cont = substr($_, 21);
	    last unless substr($_, 2, 19) eq ' 'x19;
	    chomp($ftline .= ' ' . $cont) unless $cont =~ m|^/id=\w+|;
	    push @entrylines, $_;
	}
	$ftline =~ s|\.\s*$||;
	if ($feature eq 'VARIANT') {
	    my ($from, $to);
	    if ($ftline =~ /VARIANT\s+\S+\.\.\S+\s/){
		($from, $to) = $ftline =~ /^FT\s+VARIANT\s+(\S+)\.\.(\S+)/;
	    }
	    else{
		($from) = $ftline =~ /^FT\s+VARIANT\s+(\S+)/;
		$to = $from;
	    }
	    my ($desc) = $ftline =~ /note="\/([^\/"]+)"/;
	    #Validate that the line annotates a single residue variant
	    next LINE if $from =~ /[?<>]/ || $to =~ /[?<>]/;
	    next LINE unless
		$from =~ /^\d+$/ && $to =~ /^\d+$/ &&
		$to - $from == 0;
	    if (defined $desc){
		$desc =~ s/\s$//;
		$desc =~ s/^($AA) -> ($AA)\b\s*//o;
		my ($aa1, $aa2) = ($1, $2);
		$desc =~ s/^\((.+)\)?/$1/;
		# We have an dbSNP rsID or possibly several ones annotated for the same variant
		my $dbcount;
		while ($desc =~ /dbSNP:(rs\d+)/g) {
		    push @snps, [ $1, join('', sort($aa1, $aa2)), $from, $aa1, $aa2, $desc ];
		    $dbcount++;
		}
	    }
	} else {
	    $ftline =~ s/^FT\s+//;
	    $ftline =~ s/\/evidence="([^"]+)"/$1/;
	    if ($ftline =~ /\S\.\.\S/){
		$ftline =~ s/\.\./ /;
	    }
	    else{
		$ftline =~ s/^(\S+\s+)(\S+)/$1$2 $2/;
	    }
	    push @ftlines, [ (split(' ', $ftline, 4))[0..3] ] if grep { $feature eq $_ } @fttokens;
	}
	redo LINE;
    } # end elsif
    elsif ( /^DE\s+(\w.*)\n$/ ) {
      my $w = $1;
      # skip "Contains:" description subsection(s)
      if ($w eq 'Contains:') {
        while (<TMP>) {
          next if /^DE\s{5}/; # skip all continuation lines
          redo LINE;          # restart parsing at the next section line
        }
      } else {
        # replace square brackets in definition text to avoid NCBI makeblastdb
        # parser erroneously recognizing them as (malformed) SeqID modifiers;
        # for details, see: http://www.ncbi.nlm.nih.gov/Sequin/modifiers.html
        $w =~ tr/[]/{}/;
        $de .= " $w";
      }
    } # end elsif
    elsif ( /^ / ) {
      $seq = $_;
      $seq =~ s/\s+//g;
      push @seqlines, $seq;
    } # end elsif
    elsif ( /^\/\// ) {
      if ($copy) {
        print SWALL @entrylines;
        print SWALL "//\n";
        my $primacc = $accs[0];
        my $canonical = $file =~ /uniprot_sprot/ ? 1 : 0;
        $de = ' Canonical;' . $de if $canonical;
        # truncate description if longer than 1000 symbols
        # (this is a hard-coded NCBI makeblastdb limit)
        my $lenstr = " Length: $len";
        my $maxlen = 1000 - length($lenstr);
        $de = substr($de, 0, $maxlen) if length($de) > $maxlen;
        print SEQ ">sp|$primacc|$id$de$lenstr\n";
        print SEQ join('', @seqlines), "\n";
        print FT map { $primacc . "\t" . join("\t", @$_) . "\n" } @ftlines if @ftlines;
        print ID2ACC "$primacc\t$canonical\t$id\n";
        print ID2ACC map { "$primacc\t$canonical\t$_\n" } @accs;
        print ID2ACC map { "$primacc\t$canonical\t$_\n" } @gnames if @gnames;
        print SNP2ACC map { $primacc . "\t" . join("\t", @$_) . "\n" } @snps if @snps;
        $num_filtered++;
        $dbase_len += $len;
      } # end if
    } # end elsif

    push @entrylines, $_;

  } # end while

  close TMP;

  print "Found $num_filtered known $CNAME proteins.\n";
  $total_filtered += $num_filtered;

  unlink $file or die "ERROR: Can't delete file: $file\n" unless $opts{'k'};

} # end foreach

close SNP2ACC; close ID2ACC; close SWALL; close FT; close SEQ;

print "Compressing $CNAME.dat ...\n";
`gzip -f $CNAME.dat`;

# Download VARSPLIC file with the isoform sequences and convert it
# to a greppable form (i.e., with the sequences on single line),
# then append to the main sequence file.
my $gzfile  = 'uniprot_sprot_varsplic.fasta.gz';
my $fafile  = $gzfile; $fafile =~ s/\.gz$//;
my $url  = "$UNIPROT_FTP/complete/$gzfile";
unless ($opts{'s'}) {
  print "Fetching $url ...\n";
  unlink $gzfile if -e $gzfile; # otherwise mirror() will fail
  my $rc = mirror($url, $gzfile);
  die "ERROR: Failed to fetch URL $url: ".status_message($rc)."\n" unless is_success($rc);
  print "Decompressing $gzfile ...\n";
  `gunzip -f $gzfile`;
}
# Skip VARSPLIC file processing if -s option was specified
# but uniprot_sprot_varsplic.fasta file is not present
my $varcount = 0;
unless ($opts{'s'} && ! -e $fafile) {
  print "Processing $fafile ...\n";
  my $orgtag = '_' . $SPECIES{$CNAME}{'spname'} if $SPECIES{$CNAME}{'spname'};
  open(FIN, $fafile) or die "Can't open file: $fafile\n";
  open(FOUT, ">>$seqfile") or die "Can't append to file: $seqfile\n";
  open(ACC, ">>$id2accfile") or die "Can't append to file: $id2accfile\n";
  my $seqline = '';
  my $defline = '';
  while (<FIN>) {
    chomp;
    if (/^>/) {
      if (length $defline) {
        my ($def, $desc) = split /\s+/, $defline, 2;
        if (!defined $orgtag || $def =~ /\w+$orgtag$/o) {
          print FOUT $defline, "\n";
          print FOUT $seqline, "\n" if length $seqline;
          $varcount++;
          $dbase_len += length $seqline;
          my $acc = (split /\|/, $defline, 3)[1];
          if ($acc && $acc =~ /-\d+$/) {
            print ACC "$acc\t0\t$acc\n";
          } else {
            die "Missing/malformed varsplic accession: $acc\n";
          }
        }
      }
      $defline = $_; $seqline = '';
      next;
    }
    $seqline .= $_;
  }
  if (length $defline) {
    my ($def, $desc) = split /\s+/, $defline, 2;
    if (!defined $orgtag || $def =~ /\w+$orgtag$/o) {
      print FOUT $defline, "\n";
      print FOUT $seqline, "\n" if length $seqline;
      $varcount++;
      $dbase_len += length $seqline;
      my $acc = (split /\|/, $defline, 3)[1];
      if ($acc && $acc =~ /-\d+$/) {
        print ACC "$acc\t0\t$acc\n";
      } else {
        die "Missing/malformed varsplic accession: $acc\n";
      }
    }
  }
  close(ACC);
  close(FOUT);
  close(FIN);
  print "Appended $varcount known $CNAME protein isoforms.\n" if $varcount;
  if ($opts{'k'}) {
    unless ($opts{'s'}) {
      print "Compressing $fafile ...\n";
      `gzip -f $fafile`;
    }
  } else {
    unlink $fafile or die "ERROR: Can't delete file: $fafile\n";
  }
}

print "Total number of proteins in '$CNAME' database: ", $total_filtered + $varcount, "\n";
print "Total length of protein sequences in '$CNAME' database: $dbase_len residues\n";

my $MAKEBLASTDB_NAME = 'makeblastdb';
my $mkdb_logfile     = "$CNAME.makeblastdb.log";
my $MAKEBLASTDB_ARGS = "-in $seqfile -dbtype prot -out $CNAME -title $CNAME -parse_seqids -logfile $mkdb_logfile";
my $MAKEBLASTDB_PATH;
print "Running $MAKEBLASTDB_NAME ...\n";
# Attempt to locate makeblastdb executable
my $mkmsg = "Please run the following command manually:\n\t$MAKEBLASTDB_NAME $MAKEBLASTDB_ARGS\n";
if      (-x "$ENV{PPH}/blast/bin/$MAKEBLASTDB_NAME") {
  $MAKEBLASTDB_PATH = "$ENV{PPH}/blast/bin/$MAKEBLASTDB_NAME";
} elsif (-x "$ENV{PPH}/bin/$MAKEBLASTDB_NAME") {
  $MAKEBLASTDB_PATH = "$ENV{PPH}/bin/$MAKEBLASTDB_NAME";
} else {
  my $rc = `which $MAKEBLASTDB_NAME 2>&1`; chomp $rc;
  $MAKEBLASTDB_PATH = $rc unless $rc =~ /no $MAKEBLASTDB_NAME/o;
}
die "Failed to locate $MAKEBLASTDB_NAME executable!\n$mkmsg" unless $MAKEBLASTDB_PATH;

my $cmd = "$MAKEBLASTDB_PATH $MAKEBLASTDB_ARGS";
my $rc  = `$cmd 2>&1`;
if ($?) {
  die "Failed to execute makeblastdb: $!\n$mkmsg";
} elsif ($rc) {
  die "Failed: $rc\n$mkmsg";
} else {
  local $/ = undef;
  open(LOG, "$mkdb_logfile") or
    die "Failed: Can't open $mkdb_logfile\n$mkmsg";
  my $log = <LOG>;
  close(LOG);
  $log =~ s/^\s+//s;
  die "Failed:\n$log$mkmsg" if $log =~ /^Error:/m;
}
print "Done.\n";
# Load compiled files into SQLite database
createdb($dbfile);
print "All done.\n";

#-----------------------------------------------------------------------------
#
#   Tables: id2acc, features, snp2acc, xdb2acc
#
sub createdb {
  my ($dbfile) = @_;
  die "ERROR: createdb: Missing database file name\n" unless length $dbfile;
  my @tables = qw{ id2acc features snp2acc xdb2acc };
  my @ncols  = qw{      3        5       7       3 };
  # Database schema
  my %tables = (
    id2acc =>
      q{CREATE TABLE id2acc (
          acc   VARCHAR(12) NOT NULL,
          sprot BOOLEAN NOT NULL,
          name  VARCHAR(24) NOT NULL
        )
      },
    features =>
      q{CREATE TABLE features (
          acc   VARCHAR(12) NOT NULL,
          name  VARCHAR(10) NOT NULL,
          start INTEGER NOT NULL,
          end   INTEGER NOT NULL,
          desc  TEXT
        )
      },
    snp2acc =>
      q{CREATE TABLE snp2acc (
          acc    VARCHAR(12) NOT NULL,
          rsid   VARCHAR(16) NOT NULL,
          aapair CHAR(2) NOT NULL,
          pos    INTEGER NOT NULL,
          aa1    CHAR(1) NOT NULL,
          aa2    CHAR(1) NOT NULL,
          desc   TEXT
        )
      },
    xdb2acc =>
      q{CREATE TABLE xdb2acc (
          acc   VARCHAR(12) NOT NULL,
          type  VARCHAR(16) NOT NULL COLLATE NOCASE,
          name  VARCHAR(24) NOT NULL COLLATE NOCASE
        )
      }
  );
  my %inserts = (
    id2acc   => q{INSERT INTO id2acc   VALUES (} . join(', ', ('?')x$ncols[0]) . ')',
    features => q{INSERT INTO features VALUES (} . join(', ', ('?')x$ncols[1]) . ')',
    snp2acc  => q{INSERT INTO snp2acc  VALUES (} . join(', ', ('?')x$ncols[2]) . ')',
    xdb2acc  => q{INSERT INTO xdb2acc  VALUES (} . join(', ', ('?')x$ncols[3]) . ')'
  );
  my %indices = (
    id2acc => [
      q{CREATE INDEX x_idacc  ON id2acc (acc)},
      q{CREATE INDEX x_idname ON id2acc (name)},
    ],
    features => [
      q{CREATE INDEX x_ftacc  ON features (acc)},
      q{CREATE INDEX x_ftname ON features (name)},
      q{CREATE INDEX x_ftacnm ON features (acc, name)},
      q{CREATE INDEX x_ftacst ON features (acc, start)},
      q{CREATE INDEX x_ftacen ON features (acc, end)},
    ],
    snp2acc => [
      q{CREATE INDEX x_snpacc ON snp2acc (acc)},
      q{CREATE INDEX x_snprs  ON snp2acc (rsid)},
      q{CREATE INDEX x_snppos ON snp2acc (pos)},
      q{CREATE INDEX x_snpapp ON snp2acc (acc, pos, aapair)},
    ],
    xdb2acc => [
      q{CREATE INDEX x_xdbacc  ON xdb2acc (acc)},
      q{CREATE INDEX x_xdbtype ON xdb2acc (type)},
      q{CREATE INDEX x_xdbname ON xdb2acc (name)},
      q{CREATE INDEX x_xdbtyna ON xdb2acc (type, name)},
    ],
  );
  # idmapping_selected.tab ID type names (23 total, many are optional)
  my @typenames = qw(
    UniProtKB-AC UniProtKB-ID GeneID RefSeq GI PDB GO IPI UniRef100 UniRef90 UniRef50 UniParc
    PIR NCBI-taxon MIM UniGene PubMed EMBL EMBL-CDS Ensembl Ensembl_TRS Ensembl_PRO PubMed-Extra
  );
  print 'Creating SQLite database ...';
  unlink $dbfile or die "\nERROR: Failed to remove old copy of SQLite file: $dbfile\n" if -e $dbfile;
  my $dbargs = { AutoCommit=>1, RaiseError=>1, PrintError=>0 };
  my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs);
  # Set some pragmas
  $dbh->do('PRAGMA legacy_file_format = OFF');
  $dbh->do('PRAGMA foreign_keys = ON');
  # Create tables and prepare insert statements
  my %sth;
  foreach my $table_name (@tables) {
    $dbh->do($tables{$table_name});
    $sth{$table_name} = $dbh->prepare($inserts{$table_name});
  }
  print " done.\n";
  # Populate tables
  my @datafiles = ($id2accfile, $ftfile, $snp2accfile);
  push @datafiles, $names[-1] if @names > 1 && $names[-1] =~ /idmapping_selected/; # idmapping_selected.tab
  print "Populating database (this will take a while):\n";
  my $batch_size;
  my $count;
  for (my $i=0; $i<@datafiles; $i++) {
    my $data_file  = $datafiles[$i];
    my $table_name = $tables[$i];
    next unless length $table_name;
    print 'Loading data into ', sprintf('%-8s ', $table_name);
    open(FIN, $data_file) or die "\nERROR: Failed to open file: $data_file\n";
    $batch_size = $table_name =~ /^xdb/ ? 50_000 : 5_000;
    $count = 0;
    LINE: while (<FIN>) {
      next LINE if /^#/ || /^\s*$/;
      chomp;
      my @r;
      # idmapping_selected.tab requires special processing
      if ($table_name =~ /^xdb/) {
        my $orgtag = '_' . $SPECIES{$CNAME}{'spname'} if $SPECIES{$CNAME}{'spname'};
        # For format description see:
        #   ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
        my @a = split /\t/, $_, 23;
        next LINE if defined $orgtag && $a[1] !~ /^\w+$orgtag$/o;  # skip all but CNAME proteins
        my $acc = $a[0];
        # Load the following ID types:
        # GeneID RefSeq GI PDB GO IPI UniParc PIR UniGene EMBL EMBL-CDS Ensembl Ensembl_TRS Ensembl_PRO
        foreach my $n (2..7, 11, 12, 15, 17..21) {
          my @r = split /;\s*/, $a[$n];
          next unless @r;
          map { $_ = length() ? $_ : undef } @r;
          foreach my $id (@r) {
            unless ($count % $batch_size) {
              if ($count) { print '.'; $dbh->commit; }
              $dbh->begin_work;
            }
            $id =~ s/^(ENSG[0-9]+)\.\s*\[[A-Z0-9-]+\]\s*/\1/;
            eval { $sth{$table_name}->execute($acc, $typenames[$n], $id); };
            if ($@) {
              warn $@;
              $dbh->rollback;
              print "\n" if $count;
              exit 1;
            }
            $count++;
          }
        }
      # all other tables are loaded as one data line per row
      } else {
        @r = split /\t/, $_, $ncols[$i];
        map { $_ = length() ? $_ : undef } @r;
        unless ($count % $batch_size) {
          if ($count) { print '.'; $dbh->commit; }
          $dbh->begin_work;
        }
        eval { $sth{$table_name}->execute(@r); };
        if ($@) {
          warn $@;
          $dbh->rollback;
          print "\n" if $count;
          exit 1;
        }
        $count++;
      }
    }
    close(FIN);
    print "\n";
    $dbh->commit if $count;
    # Index table
    print 'Indexing ', sprintf('%-8s', $table_name), ' 'x10;
    foreach my $index (@{ $indices{$table_name} }) {
      print '.';
      $dbh->do($index);
    }
    print "\n";
    print 'Optimizing ', sprintf('%-8s', $table_name), ' 'x8, '...';
    $dbh->do("ANALYZE $table_name");
    print "\n";
    if ($table_name =~ /^xdb/) {
	if ($opts{'k'}) {
	    print "Compressing $data_file ...\n";
	    `gzip -f $data_file`;
	} else {
	    unlink $data_file or die "ERROR: Can't delete file: $data_file\n";
	}
    } else {
	unlink $data_file or die "ERROR: Can't delete file: $data_file\n" unless $opts{'k'};
    }
  }
  $dbh->disconnect;
  print "Finished.\n";
}
