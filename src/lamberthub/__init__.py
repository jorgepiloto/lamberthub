""" A collection of Lambert's problem solvers """

from lamberthub.ecc_solvers.avanzini import avanzini2008
from lamberthub.universal_solvers.arora import arora2013
from lamberthub.universal_solvers.gooding import gooding1990
from lamberthub.universal_solvers.izzo import izzo2015

__version__ = "0.1.dev0"

ALL_SOLVERS = [gooding1990, avanzini2008, arora2013, izzo2015]
""" A list holding all lamberthub available solvers """

ZERO_REV_SOLVERS = [gooding1990, avanzini2008, arora2013, izzo2015]
""" A list holding all direct-revolution lamberthub solvers """

MULTI_REV_SOLVERS = [gooding1990, izzo2015]
""" A list holding all multi-revolution lamberthub solvers """

ROBUST_SOLVERS = [gooding1990, izzo2015]
""" A list holding all robust lamberthub solvers """
