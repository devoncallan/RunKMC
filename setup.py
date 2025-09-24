# setup.py
from setuptools import setup
from setuptools.dist import Distribution


class BinaryDistribution(Distribution):
    """Force wheel to be platform-specific, not pure Python"""

    def has_ext_modules(self):
        return True


setup(
    distclass=BinaryDistribution,
)
