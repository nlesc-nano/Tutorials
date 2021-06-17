# -*- coding: utf-8 -*-

import logging

from .__version__ import __version__ as __version__
from ._path_dict import PathMapping

logging.getLogger(__name__).addHandler(logging.NullHandler())
del logging

__author__ = "Juliette Zito"
__email__ = 'juliette.zito@hotmail.fr'

__all__ = ["PathMapping"]
