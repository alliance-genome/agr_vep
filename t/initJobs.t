use strict;
use warnings;

use Test::More tests => 7;

my $test_fasta = 't/test.fa'; 
my $test_gff = 't/test.gff'; 
my $mod = 'MGI';

my $forward_protein_seq  = 'MTKAADTKHHHSMIQRLLILLSCGYTKLLAQSPTLCCSFDVNNTFNDNVTSGLWNYEVQGEVKTVPFILNRNNKCHVTSDFENRLNATEICEKQLHSLQGQVYHFQDVLLQMRGENNTIREPLTLQSIVCGWYADERFMGSWKVCLNGSKIFHGDIKRWLHIYSGTNWTEEILEKIKNLNDFLNRTSQGEFKNKFKEYNLHCKENQEPTALSTTADVGRPSSRACTSNPSVLLIMLSCFLLYVF';
my $reverse_protein_seq1 = 'MGKSTNTPSSTWLGGETDTPISPPRCSHQFLPQPASDPSRRHQFYCSFLLLILMVSGHCSHACLQLAPVIEGNARNAGNALLVMALLLLILESCSAGTYALNCKLTVKYRTLQGLCSVNGKTFLDFGDENHEGNATMLCPALYQSLTDISEVMWSLQSGNDALNVTTRSQYYQGEFIDGFWDINTDEQHSIYVYPLNKTWRESHSDNSSAMEQWKNKNLEKDIRNVLMVDFSCCLNKSSPHFREMPR';
my $reverse_protein_seq2 = 'MGKSTNTPSSTWLGGETDTPISPPRCSHQFLPQPASDPSRRHQFYCSFLLLILMVSGHCHACLQLAPVIEGNARNAGNALLVMALLLLILESCSAGTYALNCKLTVKYRTLQGLCSVNGKTFLDFGDENHEGNATMLCPALYQSLTDISEVMWSLQSGNDALNVTTRSQYYQGEFIDGFWDINTDEQHSIYVYPLNKTWRESHSDNSSAMEQWKNKNLEKDIRNVLMVDFSCCLNKSSPHFREMPTLPTTAAHVDQPRSMACKSSPFDGLIMILLIYIL';
my $mt_seq               = 'VFFINILTLLVPILIAMAFLTLVERKILGYMQLRKGPNIVGPYGILQPFADAMKLFMKEPMRPLTTSMSLFIIAPTLSLTLALSLWVPLPMPHPLINLNLGILFILATSSLSVYSILWSGWASNSKYSLFGALRAVAQTISYEVTMAIILLSVLLMNGSYSLQTLITTQEHMWLLLPAWPMAMMWFISTLAETNRAPFDLTEGESELVSGFNVEYAAGPFALFFMAEYTNIILMNALTTIIFLGPLYYINLPELYSTNFMMEALLLSSTFLWIRASYPRFRYDQLMHLLWKNFLPLTLALCMWHISLPIFTAGVPPYM';


my $forward_id  = 'MGI_C57BL6J_3588256_transcript_2';
my $reverse_id1 = 'MGI_C57BL6J_3774845_transcript_7';
my $reverse_id2 = 'MGI_C57BL6J_3774845_transcript_6';
my $mt_id       = 'MGI_C57BL6J_101787_transcript_1';


my $self = bless {}, 'VepProteinFunction::InitJobs';
$self->param('agr_fasta', $test_fasta);
$self->param('agr_gff', $test_gff);
$self->param('debug_mode', 0);
$self->param('mod', $mod);

BEGIN { use ok 'VepProteinFunction::InitJobs' }

my ($translation_seqs, $transcript_id_translation_md5_map, $translation_md5_codontable_id) = $self->get_agr_translation_seqs;

for my $id(keys %$transcript_id_translation_md5_map){
    print( "Transcript: $id; MD5: " . $transcript_id_translation_md5_map->{$id} . "\n");
}

my $forward_md5      = $transcript_id_translation_md5_map->{$forward_id};
my $reverse_md5_1    = $transcript_id_translation_md5_map->{$reverse_id1};
my $reverse_md5_2    = $transcript_id_translation_md5_map->{$reverse_id2};
my $mt_md5           = $transcript_id_translation_md5_map->{$mt_id};
my $codontable_id    = $translation_md5_codontable_id->{$forward_md5};
my $vm_codontable_id = $translation_md5_codontable_id->{$mt_md5};

is($translation_seqs->{$forward_md5},   $forward_protein_seq,  "create forward strand translation seq");
is($translation_seqs->{$reverse_md5_1}, $reverse_protein_seq1, "create reverse strand translation seq");
is($translation_seqs->{$reverse_md5_2}, $reverse_protein_seq2, "create shared CDS translation seq");
is($translation_seqs->{$mt_md5},        $mt_seq,               "create vertebrate mitochondrial seq");

is($codontable_id,    1, "retrieve codontable ID for standard genetic code");
is($vm_codontable_id, 2, "retrieve condontable ID for vertebrate mitochondrial genetic code");

done_testing();