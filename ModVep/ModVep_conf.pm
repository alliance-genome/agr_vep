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

        pipeline_db => { ## to change
            -host   => $self->o('pipeline_host'), 
            -port   => $self->o('pipeline_port'),
            -user   => $self->o('pipeline_user'),
            -pass   => $self->o('password'),            
            -dbname => 'agr_htp_'.$self->o('pipeline_name').'_ehive',
            -driver => 'mysql',
        },

        # configuration for the various resource options used in the pipeline        
	default_lsf_options  => '-q' . $self->o('lsf_queue') . ' -R"select[mem>3000] rusage[mem=3000]" -M3000',    
        highmem_lsf_options  => '-q' . $self->o('lsf_queue') . ' -R"select[mem>8000] rusage[mem=8000]" -M8000',  
        moremem_lsf_options  => '-q' . $self->o('lsf_queue') . ' -R"select[mem>32000] rusage[mem=32000]" -M32000',  

	vep_working            => $self->o('pipeline_dir') . '/' . $self->o('mod') . '_working',
	out_file               => $self->o('pipeline_dir') . '/' . $self->o('mod') . '.vep.vcf',
        
	standard_max_workers    => 300,
	highmem_max_workers     => 25,
	hive_max_workers        => 325,
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
	vep_working         => $self->o('vep_working'),
	mod                 => $self->o('mod'),
	vcf                 => $self->o('vcf'),
    );

    return [
	{   -logic_name => 'split_input',
	    -module     => 'ModVep::SplitInput',
	    -parameters => {
		debug_mode     => $self->o('debug_mode'),
		@common_params,
	    },
	    -input_ids => [{}],
	    -rc_name   => 'default',
	    -max_retry_count => 0,
	    -flow_into => {
		'2->A' => ['run_vep'],
		'A->1' => ['combine_output'],
	    },
	},
        
        {   -logic_name     => 'run_vep',
            -module         => 'ModVep::RunVep',
            -parameters     => {
		vep_dir         => $self->o('vep_dir'),
		password        => $self->o('password'),
		gff             => $self->o('gff'),
		fasta           => $self->o('fasta'),
		bam             => $self->o('bam'),
                @common_params,
            },
            -failed_job_tolerance => 5,
            -max_retry_count => 0,
            -input_ids      => [],
            -analysis_capacity  => $self->o('standard_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name        => 'default',
            -flow_into      => {
	      2 => ['run_partial_vep'],
            },
        },

	{   -logic_name     => 'run_partial_vep',
	    -module         => 'ModVep::RunPartialVep',
	    -parameters     => {
		vep_dir         => $self->o('vep_dir'),
		password        => $self->o('password'),
		gff             => $self->o('gff'),
		fasta           => $self->o('fasta'),
		bam             => $self->o('bam'),
                @common_params,
	    },
	    -failed_job_tolerance => 10,
	    -max_retry_count => 0,
	    -input_ids       => [],
	    -analysis_capacity => $self->o('highmem_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name => 'highmem',
	    -flow_into => {
		2 => ['run_partial_vep_more_mem'],
	    },
	},

	{   -logic_name     => 'run_partial_vep_more_mem',
	    -module         => 'ModVep::RunPartialVep',
	    -parameters     => {
		vep_dir         => $self->o('vep_dir'),
		password        => $self->o('password'),
		gff             => $self->o('gff'),
		fasta           => $self->o('fasta'),
		bam             => $self->o('bam'),
                @common_params,
	    },
	    -failed_job_tolerance => 10,
	    -max_retry_count => 0,
	    -input_ids       => [],
	    -analysis_capacity => $self->o('highmem_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name => 'moremem',
	    -flow_into => {
		2 => ['run_partial_vep_more_mem'],
	    },
	},

	{   -logic_name     => 'combine_output',
	    -module         => 'ModVep::CombineOutput',
	    -parameters     => {
		out_file       => $self->o('out_file'),
		@common_params,
	    },
	    -input_ids      => [],
	    -rc_name        => 'highmem',
	},	
		
    ];
}

1;

