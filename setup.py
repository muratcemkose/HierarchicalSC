"""
Created on Mon Feb  19 14:22:31 2019

@author: Murat Cem Köse
"""

from setuptools import setup, find_packages

setup(name='HierarchicalSC',
      version='0.2',
      url='https://github.com/muratcemkose/HierarchicalSC',
      license='KUL',
      author='Murat Cem Köse',
      author_email='muratcem.kose@gmail.com',
      description='Single cell annotation library',
      packages=find_packages(),
      long_description=open('README.md').read(),
      zip_safe=False)
