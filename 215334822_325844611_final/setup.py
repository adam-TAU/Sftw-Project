from setuptools import setup, find_packages, Extension
import sysconfig

setup(
        name='spkmeans',
        version='0.0.1',
        description='TBD',
        packages = find_packages(),
        license = 'GPL-2',
        ext_modules=[
            Extension(
                'spkmeans',
                ['spkmeansmodule.c', 'spkmeans.c', 'spkmeans_goals.c', 'matrix.c', 'graph.c', 'eigen.c'],
                depends = ['spkmeans.h', 'spkmeans_goals.h', 'matrix.h', 'graph.h', 'eigen.h'],
		),
            ]
        )
