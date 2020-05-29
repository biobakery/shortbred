import sys

try:
    import setuptools
except ImportError:
    sys.exit("Please install setuptools.")
    
VERSION = "0.9.5"
AUTHOR = "Shortbred Development Team"
AUTHOR_EMAIL = "shortbred-users@googlegroups.com"
MAINTAINER = "Lauren McIver"
MAINTAINER_EMAIL = "lauren.j.mciver@gmail.com"

setuptools.setup(
    name="shortbred",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    version=VERSION,
    license="MIT",
    description="ShortBRED: a system for profiling protein families of interest at very high specificity in shotgun meta’omic sequencing data",
    long_description="ShortBRED is a system for profiling protein families of interest at very high specificity in shotgun meta’omic sequencing data. ShortBRED-Identify collapses proteins of interest into families, and then screens these families (as represented by consensus sequences) against 1) each other and 2) a comprehensive protein reference. ShortBRED-Identify then identifies short, distinguishing peptide sequences (“markers”) for each protein family. This process is performed once for a given set of proteins of interest to produce a reusable marker set (see below for examples of pre-computed markers). ShortBRED-Quantify then screens a metagenome or metatranscriptome against a given marker set to profile the presence/absence and relative abundance of the associated proteins.",
    url="http://huttenhower.sph.harvard.edu/shortbred",
    keywords=['microbial','microbiome','bioinformatics','microbiology','metagenomic','shortbred'],
    platforms=['Linux','MacOS'],
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    install_requires=["biopython"],
    packages=setuptools.find_packages(),
    scripts=['shortbred_identify.py','shortbred_quantify.py'],
    zip_safe = False
 )
