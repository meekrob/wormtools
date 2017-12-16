from setuptools import setup, find_packages
setup(
    name="wormtools",
    version="0.1",
    packages=find_packages(),
    scripts=["search_seq_motif.py"],
    author="David King",
    author_email="davidcking.inbox@gmail.com",
    description="Basic and advanced tools for bioinformatics on c. elegans.",
    keywords="python bioinformatics c.elegans",
    url="https://github.com/meekrob/wormtools"
)
