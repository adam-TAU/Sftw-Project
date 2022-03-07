from setuptools import setup, find_packages, Extension


setup(
	name='spkmeans',
	version='0.0.1',
	description='TBD',
	packages = find_packages(),
    license = 'GPL-2',
    ext_modules=[
        Extension(
            'spkmeans',
            ['spkmeansmodule.c', 'spkmeans.c', 'spkmeans.h', 'matrix.c', 'matrix.h', 'graph.c', 'graph.h', 'eigen.c', 'eigen.h'],
        )
    ]
)
