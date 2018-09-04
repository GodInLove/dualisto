from setuptools import setup, find_packages
from version import version
self_version = version

try:
    LONG_DESCRIPTION = open("README.rst", "rb").read().decode("utf-8")
except IOError:
    LONG_DESCRIPTION = "an independent demo of KNOWN operon predict method"

setup(
    name='dualisto',
    version=self_version,
    keywords=['dual-rna-seq','dual rna seq'],
    packages=find_packages(),
    url='https://github.com/GodInLove/dualisto',
    license='GPLv3',
    author='yaodongliu',
    author_email='yd.liu.scu@gmail.com',
    description="an independent demo of dual rna seq data analysis",
    long_description=LONG_DESCRIPTION,
    install_requires=[
        'biopython>=1.71',

    ],
    entry_points={
        "console_scripts": ['dualisto = dualisto.dualisto:main']
    },
)
