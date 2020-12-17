import re
from collections import OrderedDict
from ..util import toposort
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
            modules = self.sorted_modules().keys()
        else:
            # if modules were specified, then generate the complete sorted list of dependent modules 
            # that need to be dropped as well
            mods = set()
            for module in modules:
                mods.add(module)
                mods = mods | set(module.all_downstream_modules())
            modules = [m for m in self.sorted_modules().values() if m in mods]

        if job_ids is None:
            for module in reversed(modules):
                print("Dropping module %s" % module.name)
                module.drop_all()
            for module in modules:
                print("Reinitializing module %s" % module.name)
                module.initialize()
        else:
            for module in reversed(modules):
                print("Dropping %d jobs in module %s" % (len(job_ids), module.name))
                module.drop_jobs(job_ids=job_ids)
        
    def rebuild(self, modules=None, job_ids=None):
        if modules is None:
            modules = self.modules
        
    def __getstsate__(self):
        return self.kwds
        
    def report(self, modules=None, job_ids=None):
        """Return a printable string report describing the state of the pipeline plus error messages
        
        If job IDs are specified, then return information about the status of each job.
        """
        if modules is None:
            modules = self.sorted_modules()
        job_ids = job_ids or []
        mod_name_len = max([len(mod.name) for mod in modules])
            
        # only report the first error encountered for each job ID
        report = []
        failed_job_ids = set()
        totals = {}
        job_status = {}
        for module in modules:
            module_errors = []
            totals[module.name] = [0, 0]
            job_status[module.name] = module.job_status()
            for job_id, (success, error, meta) in job_status[module.name].items():
                if success is False and job_id not in failed_job_ids:
                    module_errors.append((job_id, error, meta))
                    failed_job_ids.add(job_id)
                totals[module.name][0 if success else 1] += 1

            report.append("\n=====  %s  =====\n" % module.name)
            
            # sort by error 
            module_errors.sort(key=lambda e: e[1])

            # decide exactly how to describe each error
            err_strs = []
            for job_id, job_error, meta in module_errors:
                # Attempt to summarize the job error string into a single line
                err_lines = job_error.strip().split('\n')
                
                # strip off all traceback lines
                tb_lines = []
                while len(err_lines) >= 2:
                    if err_lines[0].startswith('Traceback '):
                        err_lines = err_lines[1:]
                    if re.match(r'  File \".*\", line \d+, .*', err_lines[0]) is not None:
                        # found a traceback stack line; might be followed by a code line
                        if err_lines[1].startswith('    '):
                            tb_lines.append(err_lines[:2])
                            err_lines = err_lines[2:]
                        else:
                            tb_lines.append(err_lines[:1])
                            err_lines = err_lines[1:]
                    else:
                        break
                if len(err_lines) == 0:
                    err_msg = "[no error message]"
                else:
                    err_msg = err_lines[0]  # assume the first line after the traceback contains the most useful message
                    
                err_strs.append((
                    str(job_id),
                    "" if meta is None else meta.get('source', ''),
                    err_msg,
                ))
                
            # format into a nice table
            if len(err_strs) > 0:
                col_widths = [max([len(err[i]) for err in err_strs]) for i in range(3)]
                fmt = "{:<%ds}  {:<%ds}  {:s}\n" % tuple(col_widths[:2])
                for cols in err_strs:
                    report.append(fmt.format(*cols))
                
        report.append("\n=====  Pipeline summary  =====\n")
        fmt = "%%%ds   %%4d pass   %%4d fail\n" % mod_name_len
        for module in modules:
            tot_success, tot_error = totals[module.name]
            report.append(fmt % (module.name, tot_success, tot_error))
        
        
        for jid in job_ids:
            report.append("\n----- job: %s -----\n" % jid)
            for module in modules:
                js = job_status[module.name].get(jid, None)
                if js is None:
                    status = '-'
                    error = ''
                else:
                    success, error, meta = js
                    status = 'ok' if success else 'fail'
                    
                report_entry = "  {:20s} :  {:15s} : {:5s} : ".format(module.name, jid, status)
                
                error = error or ''
                indent = ' ' * len(report_entry)
                for i,err_line in enumerate(error.split('\n')):
                    if i > 0:
                        err_line = indent + err_line
                    report_entry = report_entry + err_line + '\n'
                    
                report.append(report_entry)
        
        return ''.join(report)
