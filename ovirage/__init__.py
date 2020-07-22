import sys
from os.path import dirname, abspath, join
from ngs_utils import logger
import numpy as np
from ovirage import polyidus


def package_path():
    return dirname(abspath(__file__))

