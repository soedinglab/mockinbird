from setuptools import setup
import os
setup(name='stammp',
      version='0.1',
      description='Analyzing PAR-CLIP data',
      url='https://bitbucket.org/soedinglab',
      author='PT',
      author_email='torkler@genzentrum.lmu.de',
      license='MIT',
      entry_points={'console_scripts': [
        'stammp-bsfinder = stammp.scripts.bsfinder:run',
        'stammp-preprocess = stammp.scripts.preprocess:run',
        'stammp-normalize = stammp.scripts.normalize:run',
        'stammp-convert2quantile = stammp.scripts.convert2quantile:run',
        'stammp-makeCenterBothEnds = stammp.plots.makeCenterBothEnds:run',
        'stammp-makeKmerPerPosition = stammp.plots.makeKmerPerPosition:run',
        'stammp-remove5primeAdapter = stammp.scripts.utils.remove5primeAdapter:run',
        'stammp-removePCRduplicates = stammp.scripts.utils.removePCRduplicates:run',
        'stammp-makeNegSets = stammp.scripts.makeNegSets:run',
        'stammp-xxmotif = stammp.scripts.xxmotif:run',
        'stammp-makeKmerLogOdds = stammp.plots.makeKmerLogOdds:run',
        'stammp-makeJaccard = stammp.plots.makeJaccard:run']
      },
      packages=['stammp','stammp.obj','stammp.scripts','stammp.plots', 'stammp.scripts.utils'],
      include_package_data=True,
      zip_safe=False)

