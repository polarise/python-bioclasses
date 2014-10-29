# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import logging

class Run( object ):
	def __init__( self, row_dict ):
		"""
		assay_type	assembly	assemblyname	biome	bioproject	biosample	center_name	collection_date	consent	env_package	estimated_size	feature	g1k_analysis_group	g1k_pop_code
		geo_loc_name	host	host_taxid	investigation_type	isol_growth_condt	lat_lon	library_name	material	mbases	mbytes	num_replicons	platform	project_name	project_type
		ref_biomaterial	run	sample_name	sequencing_method	sop	source	source_mat_id	sra_sample	sra_study	strain	
		"""
		self.sra_study = row_dict['sra_study']
		self.sra_sample = row_dict['sra_sample']
		self.run = row_dict['run']
		self.mbases = row_dict['mbases']
		self.mbytes = row_dict['mbytes']
		#self.bioproject = row_dict['bioproject']
		#self.biosample = row_dict['biosample']
		#self.sample_name = row_dict['sample_name']
		#self.library_name = row_dict['library_name']
		#self.assay_type = row_dict['assay_type']
		#self.center_name = row_dict['center_name']
		#self.consent = row_dict['consent']
		#self.assembly = row_dict['assembly']
		#self.biome = row_dict['biome']
		#self.estimated_size = row_dict['estimated_size']
		#self.feature = row_dict['feature']
		#self.host = row_dict['host']
		#self.host_taxid = row_dict['host_taxid']
		#self.material = row_dict['material']
		#self.platform = row_dict['platform']
		#self.project_name = row_dict['project_name']
		#self.project_type = row_dict['project_type']
		#self.sequencing_method = row_dict['sequencing_method']
		#self.sop = row_dict['sop']
		#self.source_mat_id = row_dict['source_mat_id']
		#self.strain = row_dict['strain']

class SraRunTable( object ):
	def __init__( self, fn, ignore="" ):	
		# meta
		self.filename = fn
		self.sra_study = None
		self.runs = dict()
		self.samples = list()
		self.no_of_runs = 0
		self.no_of_samples = 0
		self.ignore_runs = ignore.split( "," )
		
	def read( self ):
		field_names = list()
		named_study = False
		with open( self.filename ) as f:
			row_count = 0
			for row in f:
				if row_count == 0:
					field_names = row.strip( "\n" ).split( "\t" )
					row_count += 1
					continue
				else:
					l = row.strip( "\n" ).split( "\t" )
					row_dict = dict( zip( field_names, l ))
					R = Run( row_dict )
					if R.run in self.ignore_runs:
						logging.info( "Ignoring run %s..." % R.run )
					else:
						logging.info( "Adding run %s..." % R.run )
						self.runs[R.run] = R
						self.no_of_runs += 1
						if R.sra_sample not in self.samples:
							self.samples.append( R.sra_sample )
							self.no_of_samples += 1
					# get the study name
					if not named_study:
						self.sra_study = R.sra_study
						named_study = True
	
	def __repr__( self ):
		out_str = "\t".join([ "Run", "Sample", "Study", "Mbases", "MBytes" ])
		out_str += "\n"
		for run,R in self.runs.iteritems():
			out_str += "\t".join([ R.run, R.sra_sample, R.sra_study, R.mbases, \
				R.mbytes ])
			out_str += "\n"
		return out_str
		
