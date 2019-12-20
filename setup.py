from glob import glob
from setuptools import setup
from setuptools import find_packages

if __name__== "__main__":
    setup(name="ggr-project",
          version="1.1.0",
          description="GGR analysis",
          author="Daniel Kim",
          author_email='danielskim@stanford.edu',
          url="https://github.com/vervacity/ggr-project",
          license="MIT",
          install_requires=["pandas>=0.23"],
          packages=find_packages(),
          package_data={"":["data/*.json", "data/*.txt"]},
          scripts=["bin/ggr"] + glob("R/*/*.R") + glob("R/*.R")
    )
