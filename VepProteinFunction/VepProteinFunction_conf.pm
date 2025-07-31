package VepProteinFunction::VepProteinFunction_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use VepProteinFunction::Constants qw(FULL UPDATE NONE);

sub default_options {
    my ($self) = @_;

     return {

        # NB: You can find some documentation on this pipeline on confluence here:
        #
        # http://www.ebi.ac.uk/seqdb/confluence/display/EV/Protein+function+pipeline

        # Pipeline wide settings

        # If the debug_mode flag is set to 1 then we will only run the pipeline on a single transcript.
	# When set to 0 the full pipeline will be run
        
        hive_force_init                => 1,
        hive_use_param_stack           => 0,
        hive_use_triggers              => 0,
        hive_auto_rebalance_semaphores => 0, 
        hive_no_init                   => 0,
        hive_default_max_retry_count   => 0,
        hive_debug_init                => 1,
        debug_mode                     => 0,
	
        pipeline_name           => 'agr_pathogenicity_predictions_' . lc($self->o('mod')),
        pipeline_dir            => $self->o('pipeline_base_dir') . '/' . $self->o('pipeline_name'),  
        
        # directory used for the hive's own output files

        output_dir => $self->o('pipeline_dir').'/hive_output',

	# location of temporary dir, preferably with fast read/write access
	tmp_root_dir => $self->o('tmp_root_dir'),

	# peptide sequences for all unique translations for this MOD will be dumped to this file

        pep_fasta => $self->o('pipeline_dir') . '/' . $self->o('mod').'.pep.fa',
        
        # connection details for the hive's own database

        pipeline_db => { ## to change
            -host   => $self->o('pipeline_host'), 
            -port   => $self->o('pipeline_port'),
            -user   => $self->o('pipeline_user'),
            -pass   => $self->o('password'),            
            -dbname => $self->o('pipeline_name').'_ehive',
            -driver => 'mysql',
        },

	# connection details for protein function results database
	pfdb_host => $self->o('pipeline_host'), 
	pfdb_port => $self->o('pipeline_port'), 
	pfdb_user => $self->o('pipeline_user'), 
	pfdb_pass => $self->o('password'), 
	pfdb_name => 'agr_pathogenicity_predictions_' . $self->o('mod'), 

        # configuration for the various resource options used in the pipeline
	tinymem_2min_slurm_options => ' --partition=production --time=0:02:00 --mem=500m',
	default_5min_slurm_options => ' --partition=production --time=0:05:00 --mem=2000m',
	medmem_5min_slurm_options => ' --partition=production --time=0:05:00 --mem=8000m',
	highmem_12hr_slurm_options => ' --partition=production --time=12:00:00 --mem=13500m',
	supermem_2hr_slurm_options => ' --partition=production --time=2:00:00 --mem=20000m',
	supermem_12hr_slurm_options => ' --partition=production --time=12:00:00 --mem=20000m',
	highmem_48hr_slurm_options => ' --partition=production --time=2-00:00:00 --mem=13500m',
	
        sift_working            => $self->o('pipeline_dir').'/sift_working',
	pph_working             => $self->o('pipeline_dir').'/pph_working',
        
        # weka classifier models

	humdiv_model            => $self->o('pph_dir').'/models/HumDiv.UniRef100.NBd.f11.model',
	humvar_model            => $self->o('pph_dir').'/models/HumVar.UniRef100.NBd.f11.model',

        # the run type for sift (& polyphen) can be one of FULL to run predictions for
        # all translations regardless of whether we already have predictions in the
        # database, NONE to exclude this analysis, or UPDATE to run predictions for any
        # new or changed translations in the database. 

        sift_run_type           => UPDATE,
	pph_run_type            => UPDATE,
        
	# max. no. of attempts at SIFT/PolyPhen run before ignoring for update runs
	max_attempts            => 3,

	blast_max_workers       => 310,
	sift_max_workers        => 190,
	pph_max_workers         => 60,
	weka_max_workers        => 20,
	supermem_max_workers    => 20,
	hive_max_workers        => 600,
    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
	@{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation
	$self->db_cmd('CREATE TABLE IF NOT EXISTS failure_reason (
        translation_md5 char(32) NOT NULL,
        error_msg varchar(500) NOT NULL,
        analysis char(32) NOT NULL,
        PRIMARY KEY (translation_md5),
        UNIQUE KEY md5_error_analysis  (translation_md5, error_msg, analysis)
        ) ENGINE=InnoDB DEFAULT CHARSET=latin1;'),
	];
}

sub resource_classes {
    my ($self) = @_;
    return {
	'tinymem_2min'  => { 'SLURM' => $self->o('tinymem_2min_slurm_options')  },
	'default_5min'  => { 'SLURM' => $self->o('default_5min_slurm_options')  },
	'medmem_5min'   => { 'SLURM' => $self->o('medmem_5min_slurm_options')   },
	'highmem_12hr'  => { 'SLURM' => $self->o('highmem_12hr_slurm_options')  },
	'highmem_48hr'  => { 'SLURM' => $self->o('highmem_48hr_slurm_options')  },
	'supermem_2hr'  => { 'SLURM' => $self->o('supermem_2hr_slurm_options')  },
	'supermem_12hr' => { 'SLURM' => $self->o('supermem_12hr_slurm_options') },
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    
    my @common_params = (
        mod                 => $self->o('mod'),
        debug_mode          => $self->o('debug_mode'),
	pfdb_name           => $self->o('pfdb_name'),
	pfdb_host           => $self->o('pfdb_host'),
	pfdb_user           => $self->o('pfdb_user'),
	pfdb_port           => $self->o('pfdb_port'),
	pfdb_pass           => $self->o('pfdb_pass'),
	pep_fasta           => $self->o('pep_fasta'),
    );

    return [	
        {   -logic_name => 'init_jobs',
            -module     => 'VepProteinFunction::InitJobs',
            -parameters => {
		sift_working    => $self->o('sift_working'),
		pph_working     => $self->o('pph_working'),
                sift_run_type   => $self->o('sift_run_type'),
		pph_run_type    => $self->o('pph_run_type'),
		max_attempts    => $self->o('max_attempts'),
		agr_fasta       => $self->o('agr_fasta'),
		agr_gff         => $self->o('agr_gff'),
		agr_bam         => $self->o('agr_bam'),
		@common_params,
            },
	    -input_ids => [{}],
            -rc_name    => 'supermem_2hr',
            -max_retry_count => 0,
            -flow_into  => {
                2 => [ 'uniprot_align' ],
		3 => [ 'run_sift' ],
            },
        },

        {   -logic_name     => 'run_sift',
            -module         => 'VepProteinFunction::RunSift',
            -parameters     => {
		sift_run_type   => $self->o('sift_run_type'),
                sift_dir        => $self->o('sift_dir'),
                sift_working    => $self->o('sift_working'),
                blastdb         => $self->o('blastdb'),
		ncbi_dir        => $self->o('ncbi_dir'),
                @common_params,
            },
            -failed_job_tolerance => 100,
            -max_retry_count => 0,
            -analysis_capacity  => $self->o('sift_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name        => 'highmem_48hr',
            -flow_into      => {
              -1 => ['run_sift_supermem'],
            },
        },

        {   -logic_name     => 'run_sift_supermem',
            -module         => 'VepProteinFunction::RunSift',
            -parameters     => {
		sift_run_type   => $self->o('sift_run_type'),
                sift_dir        => $self->o('sift_dir'),
                sift_working    => $self->o('sift_working'),
                blastdb         => $self->o('blastdb'),
		ncbi_dir        => $self->o('ncbi_dir'),
                @common_params,
            },
            -input_ids      => [],
	    -analysis_capacity  => $self->o('supermem_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
            -rc_name        => 'supermem_12hr',
            -failed_job_tolerance => 100,
        },

	{   -logic_name     => 'uniprot_align',
	    -module         => 'VepProteinFunction::UniprotAlign',
	    -parameters     => {
		pph_dir      => $self->o('pph_dir'),
		pph_blast_db => $self->o('pph_blast_db'),
		pph_working  => $self->o('pph_working'),
		uniprot_dbs  => $self->o('uniprot_dbs'),
		@common_params,
	    },
	    -max_retry_count => 1,
	    -analysis_capacity => $self->o('blast_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name       => 'highmem_12hr',
	    -flow_into     => {
		2   => [ 'run_pph' ],
	    },
	    -failed_job_tolerance => 100,
	},

	{   -logic_name     => 'run_pph',
	    -module         => 'VepProteinFunction::RunPolyPhen2',
	    -parameters     => {
		pph_dir      => $self->o('pph_dir'),
		pph_conf_dir => $self->o('pph_conf_dir'),
		pph_working  => $self->o('pph_working'),
		@common_params,
	    },
	    -max_retry_count => 0,
	    -analysis_capacity  => $self->o('pph_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name        => 'default_5min',
	    -flow_into      => {
		2   => [ 'run_weka' ],
		-1  => [ 'run_pph_medmem' ],
	    },
            -failed_job_tolerance => 100,
	},

	{   -logic_name     => 'run_pph_medmem',
	    -module         => 'VepProteinFunction::RunPolyPhen2',
	    -parameters     => {
		pph_dir      => $self->o('pph_dir'),
		pph_conf_dir => $self->o('pph_conf_dir'),
		pph_working  => $self->o('pph_working'),
		@common_params,
	    },
	    -max_retry_count => 0,
	    -analysis_capacity  => $self->o('pph_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name        => 'medmem_5min',
	    -flow_into      => {
		2   => [ 'run_weka' ],
	    },
            -failed_job_tolerance => 100,
	},

	{   -logic_name     => 'run_weka',
            -module         => 'VepProteinFunction::RunWeka',
            -parameters     => { 
                pph_dir         => $self->o('pph_dir'),
		pph_working  => $self->o('pph_working'),
                humdiv_model    => $self->o('humdiv_model'),
		humvar_model    => $self->o('humvar_model'),
                @common_params,
            },
            -max_retry_count => 0,
            -analysis_capacity  => $self->o('weka_max_workers'),
	    -hive_capacity => $self->o('hive_max_workers'),
	    -rc_name        => 'tinymem_2min',
            -flow_into      => {},
            -failed_job_tolerance => 5,
        },
		
    ];
}

1;

