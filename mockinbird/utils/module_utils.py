import sys
import logging
from pydoc import locate

import yaml
from jinja2 import Environment, PackageLoader

from mockinbird.utils import config_validation as cv
logger = logging.getLogger()


def parse_yaml(cfg_file):
    env = Environment(loader=PackageLoader('mockinbird'))
    with open(cfg_file) as infile:
        data = env.from_string(infile.read())
        try:
            config = yaml.load(data.render())
        except yaml.scanner.ScannerError as e:
            logger.error(e)
            sys.exit(1)
    return config


def queue_pipeline(config, pipeline, def_lookup_path):
    for module_descr in config['pipeline']:
        if isinstance(module_descr, dict):
            if len(module_descr) != 1:
                logger.error('There is an error in the pipeline description syntax')
                logger.error('The pipeline section defines one list of modules. Please make'
                             ' sure you did not assign additional attributes.')
                sys.exit(1)
            # sorry for this ugly dictionary unpacking hack
            for module_str, data in module_descr.items():
                pass
            if data is None:
                data = {}
        else:
            module_str = module_descr
            data = {}
        logger.debug('Parsing configuration for module %r', module_str)

        if data.get('skip', False):
            logger.info('Skipping module %r', module_str)
            continue

        # one of our modules
        if '.' not in module_str:
            module_str = def_lookup_path + '.' + module_str
        module = locate(module_str)
        if not module:
            logger.error('module %r could not be loaded', module_str)
            sys.exit(1)
        module_name = module.__name__
        module = module(pipeline)
        config_tmpl = module.config_req
        try:
            cfg = cv.validate_section(data, config_tmpl)
        except cv.ConfigError:
            logger.error('error in config of module: %r. Exiting.', module_name)
            sys.exit(1)
        if cfg['skip']:
            mod_info = cfg['module_info']
            add_str = ' (%s)' % mod_info if mod_info else ''
            logger.debug('Skipping module %r%s', module.__class__.__name__, add_str)
        else:
            module.prepare(cfg)
            pipeline.schedule(module)


def run_pipeline(pipeline):
    for job in pipeline:
        add_str = ' (%s)' % job.module_info if job.module_info else ''
        logger.info('Executing module %r%s', job.__class__.__name__, add_str)
        res = job.execute()
        if res:
            for stdout, stderr in res:
                if stdout.strip() != '':
                    msg = 'Additional Output:\n\n'
                    logger.info(msg + stdout)
