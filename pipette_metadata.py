import yaml
import yaml_local



class PipetteMetadata(object):
    """Class for handling pipette.yml files.
    """
    def __init__(self, yml_file=None):
        self.yml_file = None
        self.meta = None
        if yml_file is not None:
            self.load_yml(yml_file)

    def load_yml(self, yml_file):
        self.meta = yaml.load(open(yml_file, 'rb'))
