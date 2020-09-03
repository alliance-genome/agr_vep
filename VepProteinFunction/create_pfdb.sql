CREATE DATABASE mqt_vep_protein_function_WB;

USE mqt_vep_protein_function_WB;

CREATE TABLE translation_md5 (
       translation_md5_id	SERIAL		PRIMARY KEY,
       translation_md5		CHAR(32)	NOT NULL	UNIQUE
);

CREATE TABLE transcript (
       transcript_id	VARCHAR(50)	PRIMARY KEY,
       translation_md5_id		BIGINT	NOT NULL	REFERENCES translation_md5(translation_md5_id)
);

CREATE TABLE analysis (
       analysis_id		SERIAL		PRIMARY KEY,
       analysis			VARCHAR(10)	NOT NULL
);

CREATE TABLE protein_function_prediction (
       translation_md5_id	BIGINT	 	 NOT NULL	REFERENCES translation_md5(translation_md5_id),
       analysis_id		SMALLINT	 NOT NULL	REFERENCES analysis(analysis_id),
       prediction_matrix	BLOB	 	 NOT NULL,
       CONSTRAINT protein_function_prediction_pk
       		  PRIMARY KEY (translation_md5_id, analysis_id)
);

CREATE TABLE valid_failure (
       translation_md5_id	BIGINT		NOT NULL	REFERENCES translation_md5(translation_md5_id),
       analysis_id		SMALLINT	NOT NULL	REFERENCES analysis(analysis_id),
       CONSTRAINT valid_failure_pk
       		  PRIMARY KEY (translation_md5_id, analysis_id)
);

CREATE TABLE analysis_attempt (
       translation_md5_id	BIGINT		NOT NULL	REFERENCES translation_md5(translation_md5_id),
       analysis_id		SMALLINT	NOT NULL	REFERENCES analysis(analysis_id),
       attempt			SMALLINT	NOT NULL,
       CONSTRAINT analysis_attempts_pk
       		  PRIMARY KEY (translation_md5_id, analysis_id)
);

INSERT INTO analysis (analysis)
VALUES
	('sift'),
	('pph');
