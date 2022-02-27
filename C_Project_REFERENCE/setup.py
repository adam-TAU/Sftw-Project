# -*- coding: utf-8 -*-

from setuptools import Extension, setup

setup(name='spkmeans',
     version='0.1.0',
     description='Python wrapper for our kmeans code',
     ext_modules=[
         Extension(
             'spkmeans',  
             sources = ['spkmeans.c','spkmeansmodule.c'],
            )
        ]
    )