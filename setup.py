from setuptools import setup, find_packages
setup(
    name="wormtools",
    version="0.1",
    package_dir={'':'src'},
    packages=find_packages('src'),
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    scripts=["scripts/search_seq_motif.py"],
    author="David King",
    author_email="David.King@colostate.edu",
    description="Python tools for bioinformatics on c. elegans.",
    keywords="python bioinformatics c.elegans",
    url="https://github.com/meekrob/wormtools"
)
