package ModVep::ModVep_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');


sub default_options {
    my ($self) = @_;

     return {

        # NB: You can find some documentation on this pipeline on confluence here:
        #
        # http://www.ebi.ac.uk/seqdb/confluence/display/EV/Protein+function+pipeline

        # Pipeline wide settings

        # If the debug_mode flag is set to 1 then we will only run the pipeline on chromosome Y.
	# When set to 0 the full pipeline will be run
        
        hive_force_init                => 1,
        hive_use_param_stack           => 0,
        hive_use_triggers              => 0,
        hive_auto_rebalance_semaphores => 0, 
        hive_no_init                   => 0,
        hive_default_max_retry_count   => 0,
        hive_debug_init                => 1,
        # the location of your ensembl checkout, the hive looks here for SQL files etc.

        pipeline_name           => $self->o('mod') . '_vep',
        pipeline_dir            => $self->o('pipeline_base_dir') . '/' . $self->o('pipeline_name'),  
        
        
        # connection details for the hive's own database
	pipeline_db_name => 'agr_htp_' . $self->o('pipeline_name') . '_ehive', 
        pipeline_db => { ## to change
            -host   => $self->o('pipeline_host'), 
            -port   => $self->o('pipeline_port'),
            -user   => $self->o('pipeline_user'),
            -pass   => $self->o('password'),            
            -dbname => $self->o('pipeline_db_name'),
            -driver => 'mysql',
        },

        # configuration for the various resource options used in the pipeline        
	default_lsf_options  => '-q' . $self->o('lsf_queue') . ' -R"select[mem>3000] rusage[mem=3000]" -M3000',    
        highmem_lsf_options  => '-q' . $self->o('lsf_queue') . ' -R"select[mem>8000] rusage[mem=8000]" -M8000',  
        moremem_lsf_options  => '-q' . $self->o('lsf_queue') . ' -R"select[mem>32000] rusage[mem=32000]" -M32000',  

	vep_working            => $self->o('pipeline_dir') . '/' . $self->o('mod') . '_working',
	out_file_prefix        => $self->o('pipeline_dir') . '/' . $self->o('mod') . '.vep.',
        
	lines_per_file => 25000,
	files_per_folder => 200,

	standard_max_workers    => 275,
	highmem_max_workers     => 75,
	hive_max_workers        => 350,

	debug => 1,

	refseq_chromosomes => {
	    'FB'   => {
		'2L' => 'NT_033779.5',
		'2R' => 'NT_033778.4',
		'3L' => 'NT_037436.4',
		'3R' => 'NT_033777.3',
		'4'  => 'NC_004353.4',
		'X'  => 'NC_004354.4',
		'Y'  => 'NC_024512.1',
		'mitochondrion_genome' => 'NC_024511.2',
		'Unmapped_Scaffold_8_D1580_D1567' => 'NW_007931083.1',
		'211000022278279' => 'NW_007931104.1',
		'211000022278436' => 'NW_001845431.1',
		'211000022278449' => 'NW_001845819.1',
		'211000022278760' => 'NW_001846712.1',
		'211000022279165' => 'NW_001846812.1',
		'211000022279188' => 'NW_001845284.1',
		'211000022279264' => 'NW_001847227.1',
		'211000022279392' => 'NW_001846198.1',
		'211000022279681' => 'NW_001845031.1',
		'211000022280328' => 'NW_001844935.1',
		'211000022280341' => 'NW_001846187.1',
		'211000022280347' => 'NW_001845870.1',
		'211000022280481' => 'NW_001845220.1',
		'211000022280494' => 'NW_001845164.1',
		'211000022280703' => 'NW_001845199.1',
		'rDNA' => 'NW_007931121.1',
	    },
	    'MGI'  => {
		# GRCm38    
		#'1'  => 'NC_000067.6',
		#'2'  => 'NC_000068.7',
		#'3'  => 'NC_000069.6',
		#'4'  => 'NC_000070.6',
		#'5'  => 'NC_000071.6',
		#'6'  => 'NC_000072.6',
		#'7'  => 'NC_000073.6',
		#'8'  => 'NC_000074.6',
		#'9'  => 'NC_000075.6',
		#'10' => 'NC_000076.6',
		#'11' => 'NC_000077.6',
		#'12' => 'NC_000078.6',
		#'13' => 'NC_000079.6',
		#'14' => 'NC_000080.6',
		#'15' => 'NC_000081.6',
		#'16' => 'NC_000082.6',
		#'17' => 'NC_000083.6',
		#'18' => 'NC_000084.6',
		#'19' => 'NC_000085.6',
		#'X'  => 'NC_000086.7',
		#'Y'  => 'NC_000087.7',
		#'MT' => 'NC_005089.1',
		# GRCm39    
		'1'  => 'NC_000067.7',
		'2'  => 'NC_000068.8',
		'3'  => 'NC_000069.7',
		'4'  => 'NC_000070.7',
		'5'  => 'NC_000071.7',
		'6'  => 'NC_000072.7',
		'7'  => 'NC_000073.7',
		'8'  => 'NC_000074.7',
		'9'  => 'NC_000075.7',
		'10' => 'NC_000076.7',
		'11' => 'NC_000077.7',
		'12' => 'NC_000078.7',
		'13' => 'NC_000079.7',
		'14' => 'NC_000080.7',
		'15' => 'NC_000081.7',
		'16' => 'NC_000082.7',
		'17' => 'NC_000083.7',
		'18' => 'NC_000084.7',
		'19' => 'NC_000085.7',
		'X'  => 'NC_000086.8',
		'Y'  => 'NC_000087.8',
		'MT' => 'NC_005089.1',
	    },
	    'RGD'  => {
		# mRatBN7.2
		'1'  => 'NC_051336.1',
		'2'  => 'NC_051337.1',
		'3'  => 'NC_051338.1',
		'4'  => 'NC_051339.1',
		'5'  => 'NC_051340.1',
		'6'  => 'NC_051341.1',
		'7'  => 'NC_051342.1',
		'8'  => 'NC_051343.1',
		'9'  => 'NC_051344.1',
		'10' => 'NC_051345.1',
		'11' => 'NC_051346.1',
		'12' => 'NC_051347.1',
		'13' => 'NC_051348.1',
		'14' => 'NC_051349.1',
		'15' => 'NC_051350.1',
		'16' => 'NC_051351.1',
		'17' => 'NC_051352.1',
		'18' => 'NC_051353.1',
		'19' => 'NC_051354.1',
		'20' => 'NC_051355.1',
		'X'  => 'NC_051356.1',
		'Y'  => 'NC_051357.1',
		'MT' => 'NC_001665.2',
		# Rnor60
		# '1'  => 'NC_005100.4',
		# '2'  => 'NC_005101.4',
		# '3'  => 'NC_005102.4',
		# '4'  => 'NC_005103.4',
		# '5'  => 'NC_005104.4',
		# '6'  => 'NC_005105.4',
		# '7'  => 'NC_005106.4',
		# '8'  => 'NC_005107.4',
		# '9'  => 'NC_005108.4',
		# '10' => 'NC_005109.4',
		# '11' => 'NC_005110.4',
		# '12' => 'NC_005111.4',
		# '13' => 'NC_005112.4',
		# '14' => 'NC_005113.4',
		# '15' => 'NC_005114.4',
		# '16' => 'NC_005115.4',
		# '17' => 'NC_005116.4',
		# '18' => 'NC_005117.4',
		# '19' => 'NC_005118.4',
		# '20' => 'NC_005119.4',
		# 'X'  => 'NC_005120.4',
		# 'Y'  => 'NC_024475.1',
		# 'MT' => 'NC_001665.2',    
	    },
	    'SGD'  => {
		'chrI'    => 'NC_001133.9',
		'chrII'   => 'NC_001134.8',
		'chrIII'  => 'NC_001135.5',
		'chrIV'   => 'NC_001136.10',
		'chrV'    => 'NC_001137.3',
		'chrVI'   => 'NC_001138.5',
		'chrVII'  => 'NC_001139.9',
		'chrVIII' => 'NC_001140.6',
		'chrIX'   => 'NC_001141.2',
		'chrX'    => 'NC_001142.9',
		'chrXI'   => 'NC_001143.9',
		'chrXII'  => 'NC_001144.5',
		'chrXIII' => 'NC_001145.3',
		'chrXIV'  => 'NC_001146.8',
		'chrXV'   => 'NC_001147.6',
		'chrXVI'  => 'NC_001148.4',
		'chrmt'   => 'NC_001224.1',
	    },
	    'WB'   => {
		'I'     => 'NC_003279.8',
		'II'    => 'NC_003280.10',
		'III'   => 'NC_003281.10',
		'IV'    => 'NC_003282.8',
		'V'     => 'NC_003283.11',
		'X'     => 'NC_003284.9',
		'MtDNA' => 'NC_001328.1',
	    },
	    'ZFIN' => {
		'1'  => 'NC_007112.7',
		'2'  => 'NC_007113.7',
		'3'  => 'NC_007114.7',
		'4'  => 'NC_007115.7',
		'5'  => 'NC_007116.7',
		'6'  => 'NC_007117.7',
		'7'  => 'NC_007118.7',
		'8'  => 'NC_007119.7',
		'9'  => 'NC_007120.7',
		'10' => 'NC_007121.7',
		'11' => 'NC_007122.7',
		'12' => 'NC_007123.7',
		'13' => 'NC_007124.7',
		'14' => 'NC_007125.7',
		'15' => 'NC_007126.7',
		'16' => 'NC_007127.7',
		'17' => 'NC_007128.7',
		'18' => 'NC_007129.7',
		'19' => 'NC_007130.7',
		'20' => 'NC_007131.7',
		'21' => 'NC_007132.7',
		'22' => 'NC_007133.7',
		'23' => 'NC_007134.7',
		'24' => 'NC_007135.7',
		'25' => 'NC_007136.7',
		'MT' => 'NC_002333.2',
	    },
	    'HUMAN' => {	
		'1'  => 'NC_000001.11',
		'2'  => 'NC_000002.12',
		'3'  => 'NC_000003.12',
		'4'  => 'NC_000004.12',
		'5'  => 'NC_000005.10',
		'6'  => 'NC_000006.12',
		'7'  => 'NC_000007.14',
		'8'  => 'NC_000008.11',
		'9'  => 'NC_000009.12',
		'10' => 'NC_000010.11',
		'11' => 'NC_000011.10',
		'12' => 'NC_000012.12',
		'13' => 'NC_000013.11',
		'14' => 'NC_000014.9',
		'15' => 'NC_000015.10',
		'16' => 'NC_000016.10',
		'17' => 'NC_000017.11',
		'18' => 'NC_000018.10',
		'19' => 'NC_000019.10',
		'20' => 'NC_000020.11',
		'21' => 'NC_000021.9',
		'22' => 'NC_000022.11',
		'X'  => 'NC_000023.11',
		'Y'  => 'NC_000024.10',
		'MT' => 'NC_012920.1'
	    },
	}
    };
}

sub pipeline_create_commands {
  my ($self) = @_;
  return [
    @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation
  ];
}

sub resource_classes {
    my ($self) = @_;
    return {
          'default'  => { 'LSF' => $self->o('default_lsf_options')  },
          'highmem'  => { 'LSF' => $self->o('highmem_lsf_options')  },
          'moremem'  => { 'LSF' => $self->o('moremem_lsf_options')  },
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    
    my @common_params = (
	vep_working        => $self->o('vep_working'),
	mod                => $self->o('mod'),
	vcf                => $self->o('vcf'),
	debug              => $self->o('debug'),
	refseq_chromosomes => $self->o('refseq_chromosomes'),
    );

    return [
	{   -logic_name => 'split_input',
	    -module     => 'ModVep::SplitInput',
	    -parameters => {
		debug_mode       => $self->o('debug_mode'),
		lines_per_file   => $self->o('lines_per_file'),
		files_per_folder => $self->o('files_per_folder'),
		gff              => $self->o('gff'),
		password         => $self->o('password'),
		@common_params,
	    },
	    -input_ids => [{}],
	    -rc_name   => 'default',
	    -max_retry_count => 0,
	    -flow_into => {
		'2->A' => ['run_vep'],
		'A->1' => ['output_factory'],
	    },
	},
        
        {   -logic_name     => 'run_vep',
            -module         => 'ModVep::RunVep',
            -parameters     => {
		vep_dir          => $self->o('vep_dir'),
		password         => $self->o('password'),
		gff              => $self->o('gff'),
		fasta            => $self->o('fasta'),
		bam              => $self->o('bam'),
                @common_params,
            },
            -failed_job_tolerance => 5,
            -max_retry_count => 0,
            -input_ids      => [],
            -analysis_capacity  => $self->o('standard_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name        => 'default',
            -flow_into      => {
	      3 => ['run_partial_vep'],
	      2 => ['process_output'],
	      -1 => ['run_partial_vep'],
            },
        },

	{   -logic_name     => 'run_partial_vep',
	    -module         => 'ModVep::RunPartialVep',
	    -parameters     => {
		vep_dir          => $self->o('vep_dir'),
		password         => $self->o('password'),
		gff              => $self->o('gff'),
		fasta            => $self->o('fasta'),
		bam              => $self->o('bam'),
                @common_params,
	    },
	    -failed_job_tolerance => 100,
	    -max_retry_count => 0,
	    -input_ids       => [],
	    -analysis_capacity => $self->o('highmem_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name => 'highmem',
	    -flow_into => {
		3 => ['run_partial_vep_more_mem'],
		2 => ['process_output'],
		-1 => ['run_partial_vep_more_mem'],
	    },
	},

	{   -logic_name     => 'run_partial_vep_more_mem',
	    -module         => 'ModVep::RunPartialVep',
	    -parameters     => {
		vep_dir          => $self->o('vep_dir'),
		password         => $self->o('password'),
		gff              => $self->o('gff'),
		fasta            => $self->o('fasta'),
		bam              => $self->o('bam'),
                @common_params,
	    },
	    -failed_job_tolerance => 0,
	    -max_retry_count => 0,
	    -input_ids       => [],
	    -analysis_capacity => $self->o('highmem_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name => 'moremem',
	    -flow_into => {
		3 => ['run_partial_vep_more_mem'],
		2 => ['process_output'],
		-1 => ['run_partial_vep_more_mem'],
	    },
	},

	{   -logic_name      => 'process_output',
            -module          => 'ModVep::ProcessOutput',
            -parameters      => {
		@common_params,
	    },
	    -input_ids       => [],
	    -rc_name         => 'default',
	    -max_retry_count => 0,
	    -failed_job_tolerance => 10,
	    -analysis_capacity => $self->o('standard_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -flow_into => {
		2 => ['process_output_highmem'],
		-1 => ['process_output_highmem'],
	    },
	},

	{   -logic_name      => 'process_output_highmem',
            -module          => 'ModVep::ProcessOutput',
            -parameters      => {
		@common_params,
	    },
	    -input_ids       => [],
	    -failed_job_tolerance => 0,
	    -rc_name         => 'moremem',
	    -analysis_capacity => $self->o('highmem_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -max_retry_count => 1,
	},

	{   -logic_name     => 'output_factory',
	    -module         => 'ModVep::OutputFactory',
	    -parameters     => {
		@common_params,
	    },
	    -flow_into      => {
		'2->A' => ['combine_output'],
		'A->1' => ['combine_chromosomes']    
             },
	},
	
	{   -logic_name     => 'combine_output',
	    -module         => 'ModVep::CombineOutput',
	    -parameters     => {
		out_file_prefix  => $self->o('out_file_prefix'),
		@common_params,
	    },
	    -input_ids      => [],
	    -rc_name        => 'highmem',
	},

	{   -logic_name     => 'combine_chromosomes',
	    -module         => 'ModVep::CombineChromosomes',
	    -parameters     => {
		out_file_prefix => $self->o('out_file_prefix'),
		@common_params,
	    },
	    -input_ids      => [],
	    -rc_name        => 'highmem'
	}	
    ];
}

1;

