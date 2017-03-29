import copy
import sys
import logging

from stammp.utils import execute
from collections import OrderedDict

import stammp.utils.config_validation as cv

logger = logging.getLogger()

class CmdPipelineModule:

    def __init__(self, pipeline, keep_all=False, cfg_req=[]):
        cfg_req.extend([
            ('keep_all', cv.Annot(cv.boolean, default=False, warn_if_missing=False)),
            ('module_info', cv.Annot(str, default='', warn_if_missing=False)),
            ('skip', cv.Annot(cv.boolean, default=False, warn_if_missing=False)),
        ])
        self._cfg_req = OrderedDict()
        for key, value in cfg_req:
            self._cfg_req[key] = value
        self._keep_all = keep_all
        self._tmp_files = []
        self._intermed_files = []
        self._cmds = []
        self._pipeline = pipeline

    def execute(self):
        for cmd in self._cmds:
            yield execute(cmd)

    def prepare(self, cfg):
        self._keep_all = cfg['keep_all']
        self._module_info = cfg['module_info']

    def cleanup(self, keep_intermed=False):
        if not self._keep_all:
            for rm_file in self._tmp_files:
                execute('rm -rf %s' % rm_file, exit=False)
            if not keep_intermed:
                for rm_file in self._intermed_files:
                    execute('rm -rf %s' % rm_file, exit=False)

    @property
    def config_req(self):
        return self._cfg_req

    @property
    def module_info(self):
        return self._module_info


class Pipeline:

    def __init__(self, initial_files, general_cfg, cfg_path):
        self._cur_files = copy.copy(initial_files)
        self._general_cfg = general_cfg
        self._jobs = []
        self._used_files = []
        self._current = 0
        self._cfg_path = cfg_path

    def schedule(self, module):
        self._jobs.append(module)

    def cleanup(self):
        n_jobs = len(self._jobs)
        for i, job in enumerate(self._jobs):
            if i < n_jobs - 1:
                job.cleanup()
            else:
                job.cleanup(keep_intermed=True)

    def __iter__(self):
        self._current = 0
        return self

    def __next__(self):
        if self._current < len(self._jobs):
            self._current += 1
            return self._jobs[self._current - 1]
        else:
            raise StopIteration

    def get_config(self, cfg_name):
        return self._general_cfg.get(cfg_name)

    def has_curfile(self, fmt):
        filepath = self._cur_files.get(fmt)
        return filepath is not None

    def get_curfile(self, fmt):
        filepath = self._cur_files.get(fmt)
        if filepath is None:
            logger.error('No module created a file of type %r', fmt)
            sys.exit(1)
        self._used_files.append(filepath)
        return filepath

    def upd_curfile(self, fmt, filepath):
        self._cur_files[fmt] = filepath

    @property
    def cfg_path(self):
        return self._cfg_path
