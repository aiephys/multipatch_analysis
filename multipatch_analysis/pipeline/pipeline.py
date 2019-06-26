from collections import OrderedDict
from pyqtgraph import toposort
from .pipeline_module import PipelineModule, DatabasePipelineModule


class Pipeline(object):
    """A sequence of pipeline modules and methods for invoking them.
    """
    name = None
    module_classes = []
    
    @classmethod
    def all_pipelines(cls):
        pipelines = OrderedDict()
        for pcls in cls.__subclasses__():
            if pcls.name is not None:
                pipelines[pcls.name] = pcls
            pipelines.update(pcls.all_pipelines())
        return pipelines
    
    def __init__(self, **kwds):
        self.kwds = kwds
        self.modules = [mcls(self) for mcls in self.module_classes]

        excluded = [PipelineModule, DatabasePipelineModule]
        deps = {c:c.upstream_modules() for c in self.modules if type(c) not in excluded}
        self._sorted_modules = OrderedDict([(mod.name, mod) for mod in toposort(deps)])

    def sorted_modules(self):
        """Return an ordered dictionary mapping {name: module} of all modules in this pipeline,
        topologically sorted by dependencies (least dependent to most dependent).
        """
        return self._sorted_modules
    
    def get_module(self, module_name):
        return self.sorted_modules()[module_name]
        
    def update(self, modules=None, job_ids=None):
        if modules is None:
            modules = self.modules

        
    def drop(self, modules=None, job_ids=None):
        if modules is None:
            modules = self.modules
        for module in modules:
            if job_ids is None:
                print("Dropping and reinitializing module %s" % module.name)
                module.drop_all(reinitialize=True)
            else:
                print("Dropping %d jobs in module %s" % (len(job_ids), module.name))
                module.drop_jobs(job_ids=job_ids)
        
    def rebuild(self, modules=None, job_ids=None):
        if modules is None:
            modules = self.modules
        
    def __getstsate__(self):
        return self.kwds
        
