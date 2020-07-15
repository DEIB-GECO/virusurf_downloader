class VirusSource:
    def __init__(self):
        pass

    def taxon_id(self):
        raise NotImplementedError

    def taxon_name(self):
        raise NotImplementedError

    def family(self):
        raise NotImplementedError

    def sub_family(self):
        raise NotImplementedError

    def genus(self):
        raise NotImplementedError

    def species(self):
        raise NotImplementedError

    def equivalent_names(self):
        raise NotImplementedError

    def molecule_type(self):
        raise NotImplementedError

    def is_single_stranded(self):
        raise NotImplementedError

    def is_positive_stranded(self):
        raise NotImplementedError

    def virus_samples(self):
        raise NotImplementedError
