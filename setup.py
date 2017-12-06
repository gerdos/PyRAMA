from setuptools import setup

setup(
    name='pyrama',
    version='2.0.1',
    scripts=['bin/pyrama'],
    packages=['rama_lib'],
    package_data={'rama_lib': ['data/*.data']},
    url='https://github.com/gerdos/PyRAMA',
    license='',
    author='gerdos',
    author_email='gerdos@caesar.elte.hu',
    description='Ramachandran plot generator',
    install_requires=['biopython', 'numpy', 'matplotlib']
)
