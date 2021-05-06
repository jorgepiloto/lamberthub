""" A collection of Lambert's problem solvers """

from lamberthub.universal_solvers.gooding import gooding1990

__version__ = "0.1.dev0"

ALL_SOLVERS = [gooding1990]
""" A list holding all lamberthub available solvers """

ZERO_REV_SOLVERS = [gooding1990]
""" A list holding all direct-revolution lamberthub solvers """

MULTI_REV_SOLVERS = [gooding1990]
""" A list holding all multi-revolution lamberthub solvers """

ROBUST_SOLVERS = [gooding1990]
""" A list holding all robust lamberthub solvers """
