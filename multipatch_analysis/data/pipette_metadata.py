import os
from collections import OrderedDict
import yaml
from .. import yaml_local


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
            pipettes = yaml.load(open(yml_file, 'rb'), Loader=yaml.FullLoader)
        else:
            # pyyaml old API
            pipettes = yaml.load(open(yml_file, 'rb'))

        # convert integer pipette IDs to str, to match data models downstream
        self.pipettes = OrderedDict()
        for k,v in pipettes.items():
            self.pipettes[str(k)] = v
            if v.get('synapse_to') is not None:
                v['synapse_to'] = [str(x) for x in v['synapse_to']]
            if v.get('gap_to') is not None:
                v['gap_to'] = [str(x) for x in v['gap_to']]
            