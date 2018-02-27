# This function returns the parameters that will be needed to run
# the program.
def set_parameters():

    # Parameters the user needs to specify:

    # Path to the data file
    path_to_data = ''

    # Error flag.
    # If there are errors on the velocity data this should be set to True
    # If there aren't (e.g. in simulated data) this should be False.
    error_flag = True

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Parameters the user may wish to specify, or can leave at the
    # default values

    # Range of dr to colour code stars by in figure 2.
    dr_start = 0.
    dr_end = 1.

    # The width of the bins dr pairs are sorted into.
    bin_width = 0.1

    '''
    The default definition of dv is the magnitude of the
    velocity difference, i.e the velocity difference between
    two stars is calculated in the same way you'd calculate
    the distance between two stars. The result is always positive.

    This method will be used if dv_default is set to True.

    An alternative way to define dv is the rate at which the
    distance between the stars is changing, so dr/dt. If the
    distance between the stars is decreasing (so moving towards
    one another) this value is negative. If they are moving
    further apart the distance between the strs increases, so
    the value returned is positive.

    This definition is helpful for differentiating between
    expanding and collapsing clusters/regions. It will be used
    if dv_default is set to False.
    '''
    dv_default = True

    '''
    Observational errors artificially increase the mean dv of all
    pairs measured. This has little impact on analysis of velocity
    structure, nevertheless it is possible to correct for this effect.
    If this flag is set to True then the true mean dv will be estimated,
    and the plot produced will use this corrected value.
    If this flag is set to False no correction will be done.

    Doing the correction may take some time (a few minutes for 1000 stars).
    The correction assumes that the uncertainty on every velocity
    measurement is the same.
    '''
    correct_inflation = False

    return path_to_data, error_flag, dv_default, correct_inflation, bin_width, dr_start, dr_end
