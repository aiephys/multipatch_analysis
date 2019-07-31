import os
import yaml
from . import yaml_local


class PipetteMetadata(object):
    """Class for handling pipette.yml files.
    """
    def __init__(self, site_dir=None):
        self.site_dir = site_dir
        self.pipettes = None
        if site_dir is None:
            self.yml_file = None
        else:
            yml_file = os.path.join(site_dir, 'pipettes.yml')
            if os.path.isfile(yml_file):
                self.load_yml(yml_file)

    def load_yml(self, yml_file):
        self.yml_file = yml_file

        if hasattr(yaml, 'FullLoader'):
            # pyyaml new API
            self.pipettes = yaml.load(open(yml_file, 'rb'), Loader=yaml.FullLoader)
        else:
            # pyyaml old API
            self.pipettes = yaml.load(open(yml_file, 'rb'))

