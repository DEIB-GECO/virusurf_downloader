from sqlalchemy import create_engine  
from sqlalchemy import Column, ForeignKey
from sqlalchemy import String, Integer, Boolean, Float, Date
from sqlalchemy.ext.declarative import declarative_base  
from sqlalchemy.orm import sessionmaker

# https://www.compose.com/articles/using-postgresql-through-sqlalchemy/



base = declarative_base()


class ExperimentType(base):  
    __tablename__ = 'experiment_type'
    
    experiment_type_id = Column(Integer, primary_key=True)
                            
    sequencing_technology = Column(String)
    assembly_method = Column(String)
    coverage = Column(String)
                            

class SequenceProject(base):  
    __tablename__ = 'sequence_project'
    
    sequence_project_id = Column(Integer, primary_key=True)
                            
    sequencing_lab = Column(String)
    submission_date = Column(Date)
    bioproject_id = Column(String)
    database_source = Column(String)
    
                                    
class Virus(base):
    __tablename__ = 'virus'
    #     
    virus_taxonomy_id = Column(Integer, primary_key=True, autoincrement=False)
    
    family = Column(String)
    sub_family = Column(String)
    genus  = Column(String)
    species_name = Column(String)
    species_taxon_id = Column(String)
    genbank_acronym = Column(String)
    equivalent_list = Column(String)
    molecule_type = Column(String)
    is_single_stranded = Column(String)
    is_positive_stranded = Column(String)
    
    
class HostSample(base):
    __tablename__ = 'host_sample'
    
    host_sample_id = Column(Integer, primary_key=True)
    
    species = Column(String)
    species_taxon_id = Column(Integer)
    collection_date = Column(Integer)
    isolation_source  = Column(String)
    originating_lab = Column(String)
    country = Column(String)
    region = Column(String)
    geo_group = Column(String)

    


class Sequence(base):  
    __tablename__ = 'sequence'
    
    sequence_id = Column(Integer, primary_key=True)
    experiment_type_id = Column(Integer, ForeignKey(ExperimentType.experiment_type_id), nullable=False)
    sequence_project_id = Column(Integer, ForeignKey(SequenceProject.sequence_project_id), nullable=False)

    accession_id = Column(String, unique=True, nullable=False)
    alternative_accession_id = Column(String, unique=True, nullable=False)
    strain_name = Column(String)
    is_reference = Column(Boolean, nullable=False)
    is_complete = Column(Boolean)
    nucleotide_sequence = Column(String)
    strand = Column(String)
    length = Column(Integer)
    gc_percentage = Column(Float)


class Annotation(base):  
    __tablename__ = 'annotation'
 
    annotation_id = Column(Integer, primary_key=True)
    sequence_id = Column(Integer, ForeignKey(Sequence.sequence_id), nullable=False)

    feature_type = Column(String, nullable=False)
    start = Column(Integer)
    stop = Column(Integer)
    gene_name = Column(String)
    product = Column(String)
    external_reference = Column(String)
    
    
class Variant(base):  
    __tablename__ = 'variant'
 
    variant_id = Column(Integer, primary_key=True)
    sequence_id = Column(Integer, ForeignKey(Sequence.sequence_id), nullable=False)
    
    start = Column(Integer)
#     start2 = Column(Integer)
    length = Column(Integer, nullable=False)
    original_sequence = Column(String, nullable=False)
    alternative_sequence = Column(String, nullable=False)
    variant_type = Column(String, nullable=False)

    def get_list(self):
        return [self.start, self.length, self.original_sequence, self.alt_sequence, self.variant_type]
    def get_list_columns():
        return ['start', 'length', 'original_sequence', 'alt_sequence',  'variant_type']
    