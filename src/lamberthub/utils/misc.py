""" Miscellany utilities """


def get_solver_name(solver):
    """
    Returns solver's name.

    Parameters
    ----------
    solver: function
        Solver function.

    Returns
    -------
    name: str
        String representation for the solver.

    """

    # Get the function name and its length
    raw_name = str(solver.__name__.capitalize())
    len_name = len(raw_name)

    # Break into author name and implementation year
    author_name = raw_name[0 : (len_name - 4)]
    year_name = raw_name[(len_name - 4) :]

    # Assemble and return the solver's name
    name = author_name + " " + year_name
    return name
