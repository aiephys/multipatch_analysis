"""
Add custom behavior to YAML serializer:

* Load and save OrderedDict as a regular dict, but with key order preserved.

"""
import sys
import yaml
import collections

if sys.version[0] > '2':
    unicode = str

# Credit: https://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts
_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def unicode_representer(dumper, data):
    return dumper.represent_str(str(data))

def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))

yaml.add_representer(collections.OrderedDict, dict_representer)
yaml.add_representer(unicode, unicode_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)

