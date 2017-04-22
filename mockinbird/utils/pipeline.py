import copy
import sys
import logging

from mockinbird.utils import execute
from collections import OrderedDict

import mockinbird.utils.config_validation as cv

logger = logging.getLogger()


class CmdPipelineModule:
    """Parent module for all modules that run shell commands

    This class stores a list of commands in the ``_cmds`` variable and executes them on
    demand. It defines reasonable default behavior for all necessary pipeline methods.

    Defining your own module requires following steps:

    - overwrite the constructor and pass the configuration requirements of the module
    - overwrite the prepare method and
        - add file paths to temporary files to ``_tmp_files``
        - add intermediate file paths to ``_intermed_files``
        - add shell commands to ``_cmds``
    - use ``_pipeline`` to access the paths to previously generated files and
      register files as module output
    """

    def __init__(self, pipeline, keep_all=False, cfg_req=[]):
        """Initializes a new module with empty list of queued commands

        Args:
            pipeline: the pipeline the module is queued in
            keep_all (:obj:`boolean`): do not remove any temporary files
            cfg_req (:obj:`list`): list of configuration options
        """
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
        """Execute all queued commands"""
        for cmd in self._cmds:
            yield execute(cmd)

    def prepare(self, cfg):
        """Prepares the module by queuing commands to execute

        This is the main method that queues commands and lists files that are to be cleaned up.
        Subclasses should call the parent method before providing their own implementation.

        Args:
            cfg (:obj:`dict`): dictionary of the configuration options
        """
        self._keep_all = cfg['keep_all']
        self._module_info = cfg['module_info']

    def cleanup(self, keep_intermed=False):
        """Cleans up temporary and intermediate files

        Args:
            keep_intermed (:obj:`boolean`): if set to True, intermediate files are not removed
        """
        if not self._keep_all:
            for rm_file in self._tmp_files:
                logger.debug('%s is marked as a temporary file, cleaning up', rm_file)
                execute('rm -rf %s' % rm_file, exit=False)
            if not keep_intermed:
                for rm_file in self._intermed_files:
                    logger.debug('%s is marked as an intermediate file, cleaning up', rm_file)
                    execute('rm -rf %s' % rm_file, exit=False)

    @property
    def config_req(self):
        """The list of config requirements"""
        return self._cfg_req

    @property
    def module_info(self):
        """The module info string"""
        return self._module_info


class Pipeline:
    """The Pipeline queues modules and provides means for modules to communicate"""

    def __init__(self, initial_files, general_cfg, cfg_path):
        """
        Construct a new pipeline

        Args:
            initial_files (:obj:`dict`): dictionary mapping containing initial
                                         (file type, file path) key-value pairs
            general_cfg (:obj:`dict`): dictionary of dictionaries of globally available
                                       config options
            cfg_path (:obj:`str`): file path to the config file
        """

        self._cur_files = copy.copy(initial_files)
        self._general_cfg = general_cfg
        self._jobs = []
        self._used_files = []
        self._current = 0
        self._cfg_path = cfg_path

    def schedule(self, module):
        """Append module to the list of schedule modules"""
        self._jobs.append(module)

    def cleanup(self):
        """Clean up files of all modules

        This method invokes the ``cleanup`` method of all queued modules.
        Intermediate files of the last module are not cleaned up
        """
        n_jobs = len(self._jobs)
        for i, job in enumerate(self._jobs):
            if i < n_jobs - 1:
                job.cleanup()
            else:
                job.cleanup(keep_intermed=True)

    def __iter__(self):
        """Iterate over all queued modules"""
        self._current = 0
        return self

    def __next__(self):
        if self._current < len(self._jobs):
            self._current += 1
            return self._jobs[self._current - 1]
        else:
            raise StopIteration

    def get_config(self, cfg_name):
        """Get config section from global config dictionary

        Args:
            cfg_name (:obj:`str`): name of the config section to retrieve
        """
        return self._general_cfg.get(cfg_name)

    def has_curfile(self, fmt):
        """Check if a file of format ``fmt`` is already available in the pipeline"""
        filepath = self._cur_files.get(fmt)
        return filepath is not None

    def get_curfile(self, fmt):
        """Get the path to the most recently created file of format ``fmt``

        Args:
            fmt (:obj:`str`): file format of file to retrieve

        Exits the program with exit-code ``1`` if the file was not yet queued
        """
        # TODO exiting here is really ugly. Better raise an exception
        filepath = self._cur_files.get(fmt)
        if filepath is None:
            logger.error('No module created a file of type %r', fmt)
            sys.exit(1)
        self._used_files.append(filepath)
        return filepath

    def upd_curfile(self, fmt, filepath):
        """register ``filepath`` as the most recently created file of format ``fmt``"""
        self._cur_files[fmt] = filepath

    @property
    def cfg_path(self):
        """path to the config file"""
        return self._cfg_path
