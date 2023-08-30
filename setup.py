from setuptools import setup, find_packages
from sieveai import __program__, __version__, __subversion__

REQUIREMENTS = [i.strip() for i in open("sieveai.requirements.txt").readlines()]

LONG_DESCRIPTION = "".join(open("README.md").readlines())

setup(
  name=__program__,
  version=f"{__version__}.{__subversion__}",
  packages=find_packages(),
  description='An automated drug discovery pipeline',
  author='Vishal Kumar Sahu',
  author_email='mail@vishalkumarsahu.in',
  url='https://www.github.com/VishalKumarSahu/SieveAI',
  install_requires=REQUIREMENTS,
  entry_points={
    'console_scripts': ['sieveai=sieveai:dock', 'sieveai-rescore=sieveai:rescore'],
  },
  long_description=LONG_DESCRIPTION,
  long_description_content_type='text/markdown',
  classifiers=[
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Topic :: Scientific/Engineering :: Molecular Docking',
  ],
)
