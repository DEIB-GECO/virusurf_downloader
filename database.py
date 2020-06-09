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
                            

class SequencingProject(base):  
    __tablename__ = 'sequencing_project'
    
    sequencing_project_id = Column(Integer, primary_key=True)
                            
    sequencing_lab = Column(String)
    submission_date = Column(Date)
    database_source = Column(String)
    bioproject_id = Column(String)
    
                                    
class Virus(base):
    __tablename__ = 'virus'
    #  
    virus_id = Column(Integer, primary_key=True)
    
    taxon_id = Column(Integer)
    taxon_name = Column(String)
    
    family = Column(String)
    sub_family = Column(String)
    genus  = Column(String)
    species = Column(String)
    equivalent_list = Column(String)
    molecule_type = Column(String)
    is_single_stranded = Column(Boolean)
    is_positive_stranded = Column(Boolean)
    
    
class HostSample(base):
    __tablename__ = 'host_sample'
    
    host_sample_id = Column(Integer, primary_key=True)
    
    host_taxon_id = Column(Integer)
    host_taxon_name = Column(String)
    
    collection_date = Column(String)
    isolation_source  = Column(String)
    originating_lab = Column(String)
    country = Column(String)
    region = Column(String)
    geo_group = Column(String)
    
#     extra 
    age = Column(Integer)
    gender = Column(String)
    

    


class Sequence(base):  
    __tablename__ = 'sequence'
    
    sequence_id = Column(Integer, primary_key=True)
    # FKs     
    experiment_type_id = Column(Integer, ForeignKey(ExperimentType.experiment_type_id), nullable=False)
    virus_id = Column(Integer, ForeignKey(Virus.virus_id), nullable=False)
    host_sample_id = Column(Integer, ForeignKey(HostSample.host_sample_id), nullable=False)
    sequencing_project_id = Column(Integer, ForeignKey(SequencingProject.sequencing_project_id), nullable=False)
    
    accession_id = Column(String, unique=True, nullable=False)
    alternative_accession_id = Column(String, unique=True, nullable=True)
    strain_name = Column(String)
    is_reference = Column(Boolean, nullable=False)
    is_complete = Column(Boolean)
    nucleotide_sequence = Column(String)
    strand = Column(String)
    length = Column(Integer)
    gc_percentage = Column(Float)
    linage = Column(String)
    clade = Column(String)


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
    
    sequence_original = Column(String, nullable=False)
    sequence_alternative = Column(String, nullable=False)
    start_original = Column(Integer)
    start_alternative = Column(Integer)
    variant_length = Column(Integer, nullable=False)
    variant_type = Column(String, nullable=False)

    def get_list(self):
        return [self.start_original, self.variant_length, self.sequence_original, self.sequence_alternative, self.variant_type]
    def get_list_columns():
        return ['start', 'length', 'sequence_original', 'alt_sequence',  'variant_type']
    