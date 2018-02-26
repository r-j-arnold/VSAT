# Version 1.0, Becky Arnold, January 2018 running on python version 2.7.13
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math
import fort_subroutines
import github_correct_dv


'''
This program takes the positions and velocities of stars
in a cluster and produces a plot of the cluster's velocity
structure as per the method in Arnold et al. (2018) (submitted).

To help users analyse and investigate the velocity structure
a 2nd plot is produced. This plot is the positions of the stars
with the stars colour coded according to their contribution to
different part of the cluster's velocity structure. This is described
in more detail in

SECTION OF PAPER!

The user can choose which part of the cluster's velocity structure to
colour code by by modifying the parameters dr_start and dr_end.

This program is not dependant on the units of the data inputted,
although pc and km/s although are assumed and used to label figure
axes.

The user must fill in the read_data function to import their
data as different people format their data files differently.
Other than that the user only needs to set the error_flag (below).
It should be set to true if there are observational errors and false
if there are not (as may be the case with data from simulations).

The code is split into sections with ~~~~~~~~~~~~~~ breaks
The sections are as follows:

1. Read data function
2. Read in the data
3. Set parameters
4. Calculate dr and dv for every pair of stars
5. Sort pairs into dr bins and calculate the average dv in each bin with errors
6. Tidy up
7. Get rid of any empty bins
8. Adjust for inflation of average dv due to uncertainties
9. Make plots of the data- 1st plot: velocity structure
10. 2nd plot- positions of the stars with the stars colour coded
'''


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1. READ IN DATA FUNCTION

# The user must write this function
def read_data():

    '''
    This is the only function the user needs to write.
    It reads in the position (r), velocity (v) and
    (if there are any) velocity errors (verr) for all stars.
    It is necessary for the user to write this function because
    people format their data files differently.
    The output of this function should be lists, each containing
    lists of the form:

    r = [[x positions of stars], [y positions of stars], [z positions of stars]]
    v = [[vx velocities of stars], [vy velocities of stars], [vz velocities of stars]]
    and (if errors)
    verr = [[measurement errors on vx velocities], [measurement errors on vy velocities], [measurement errors on vz velocities]]

    If the data is observational rather than simulated it is unlikely the full
    3D information will be available.
    If data is absent them omit the relevant lists.
    For example there may only be 2D positional data
    and 1D radial velocity data available. In that case the output
    lists from this function would be:
    r = [[x positions of stars], [y positions of stars]]
    v = [[radial velocities of stars]]
    verr = [[measurement errors on radial velocities]]
    If there are no measurement errors then either don't return verr
    or return it as zero.
    '''


    return r, v

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2. READ IN THE DATA

# Read in the position (r), velocity (v), and,
# if there is any, velocity error (verr) data.
# This is the only function the user must write.
r, v = read_data()

# The user must set the error flag.
# If there are errors on the velocity data this should be set to true
# If there aren't (e.g. in simulated data) this should be false.
error_flag = False


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 3. SET PARAMETERS

# Parameters the user may wish to specify, or can leave at the default values

# For the colour coded plot
dr_start = 0.
dr_end = 1.

'''
The default definition of dv is the magnitude of the
velocity difference, i.e the velocity difference between
two stars is calculated in the same way you'd calculate
the distance between two stars. The result is always positive.

This method will be used if dv_default is set to true

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
Observational errors artificially increase the average dv
measured. Average dv has little impact on analysis of velocity
structure. Nevertheless it is possible to correct for this effect.
If this flag is set to true then the true average dv will be estimated,
and the plot produced will use this corrected value.
Doing the correction may take some time (a few minutes for 1000 stars).
The correction assumes that the uncertainty on every velocity
measurement is the same. The
'''
correct_inflation = False

# The width of the bins dr pairs are sorted into.
bin_width = 0.1







# ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

# No further modification by the user should be necessary beyond this point
# (although users are welcome to modify anything to better suit their needs
# or make improvements).

# ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !








# Count the number of stars in the data
n_stars = len(r[0])

# For this number of stars how many possible pairs are there?
n_pairs = ((n_stars*n_stars) - n_stars)/2.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 4. CALCULATE DR AND DV FOR EVERY PAIR OF STARS


'''
For every possible pair of stars calculate the distance between them, dr
and the difference in their velocities, dv.
If there are errors on the velocities then propagate them to get the
errors on the dvs

A fortran subroutine has been written to do this step because
fortran is significantly faster than python at number crunching.
The time saving due to this is small for a single snapshot
however for a many snapshots (e.g. if analysing every timestep of a simulation)
then it becomes significant.
The package f2py was used to make the fortran subroutine compatible with python.
If for any reason the user wishes to make changes to this subroutine then
they will need to build their altered version with f2py

This subroutine returns an array containing the dr and corresponding dv
for every pair of stars, and the error on dv (if there are velocity errors).
If there are no velocity errors then that column is all zeros.
It also returns the largest dr in the cluster so
we know the biggest dr bin that will be needed.
'''


if error_flag:
    dr_dv, max_dr = fort_subroutines.calc_dr_dv(n_stars, n_pairs, bin_width, len(r), r, len(v), v, error_flag, dv_default, verr)
else:
    dr_dv, max_dr = fort_subroutines.calc_dr_dv(n_stars, n_pairs, bin_width, len(r), r, len(v), v, error_flag, dv_default)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5. SORT PAIRS INTO DR BINS AND CALCULATE THE AVERAGE DV IN EACH BIN

# Define the bin edges, these will mostly be used for
# plotting, and define the number of bins
edges = list(np.arange(0, max_dr, bin_width))
n_bins = len(edges)

'''
As in the above section a fortran subroutine is used for
this step. Again, this is done because fortran is faster
than python at number crunching, and if the user wishes
to make changes to this subroutine then they will need
to build their altered version with f2py.

This subroutine returns the mean dv in each bin, held
in the array called mean_dv. It also returns the errors
on the average dv in each bin, called error. It returns
the number of pairs in each bin, n_in_bins, which will
be used in the next step to tidy up the figure. Finally
it returns a list detailing how many times each star
appears in each dv bin. E.g. if a certain star is part
of many small dr pairs but few large dr ones it will
appear often in low dr bins. This list is called
count_stars_bins, and is used to colour code plots to
better understand how different stars contribute to
velocity structure.
'''

if error_flag:
    mean_dv, error, n_in_bins, count_stars_bins = fort_subroutines.sort_errors(n_stars, n_pairs, n_bins, dr_dv, bin_width)
else:
    mean_dv, error, n_in_bins, count_stars_bins = fort_subroutines.sort_no_errors(n_stars, n_pairs, n_bins, dr_dv, bin_width)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 6. TIDY UP


'''
These changes are just to make the figure look tidier.

High dr bins are sparsely populated because only stars
on the extreme outskirts are far enough apart to fall
into these high dr bins. By definition there are few
stars on the outskirts of a cluster, resulting in there
being few pairs in these bins. If a bin contains few
pairs its average dv is unreliable and often noisy.
Therefore we get rid of bins at high dr. Define high
dr as the dr beyond which the average number of pairs
per bin is < 10. Get rid of the bins above that cutoff
'''


# Go to larger and larger dr until the cutoff is reached
for i in range(n_bins):

    # Get the average number of pairs per bin for bins with dr >= this one
    av_numb_in_bins = sum(n_in_bins[i:])/float(n_bins - i)

    # It the average number of pairs per bin is < 10 then the bins
    # beyond this point are likely noisy and unreliable
    # Therefore want to cut any bins beyond this cutoff point
    if av_numb_in_bins < 10:

        cutoff_bin = i

        break

# Get rid of the bins above the cutoff_point
mean_dv = mean_dv[:cutoff_bin]
error = error[:cutoff_bin]
n_in_bins = n_in_bins[:cutoff_bin]
edges = edges[:cutoff_bin]

# The number of bins has changed now the outermost have been cut. Recalculate.
n_bins = len(mean_dv)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 7. GET RID OF ANY EMPTY BINS

# Figure out which bins to remove
# Reversing the order just makes it easier to delete them without
# screwing up the indecies
remove_list = [i for i in range(n_bins) if n_in_bins[i] == 0.]
remove_list.reverse()

# Remove the empty bins
for empty in remove_list:

    mean_dv = np.delete(mean_dv, empty)
    error = np.delete(error, empty)
    n_in_bins = np.delete(n_in_bins, empty)
    edges = np.delete(edges, empty)

# The number of bins has changed now the empty ones have been
# deleted. Recalculate.
n_bins = len(edges)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 8. ADJUST FOR INFLATION OF AVERAGE DV DUE TO UNCERTAINTIES

# Observational errors artificially increase average dv measured.
# Correct this inflation of dv if there are errors and the user
# has set the correct_inflation flag to true. Note that this
# function assumes the uncertainty on each velocity measurement
# is the same.

if error_flag and correct_inflation:

    mean_dv = github_correct_dv.correct(n_stars, v, verr[0][0], mean_dv)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 9. MAKE PLOTS OF THE DATA- 1ST PLOT: VELOCITY STRUCTURE

# Cosmetic tweaks to the font style for the figure
mpl.rcParams['font.size'] = '18'
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Plot the velocity structure
plt.figure(1)
plt.errorbar(edges, mean_dv, yerr=error, capsize=3.)

# Add axis labels, etc
plt.xlabel('$\Delta r$ (pc)')
plt.ylabel('$\Delta v~ $' + '(km s' + r'$^{-1}$)')
plt.title('Velocity structure')

# Neaten up the display
plt.xlim([0, edges[-1]])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 10. 2ND PLOT: POSITIONS OF THE STARS WITH THE STARS COLOUR CODED

# Can only do the plot if there's at least 2D spatial information
if len(r) > 1:

    '''
    This is a tool to help investigate and understand different
    parts of the velocity structure shown in Fig. 1.

    Say Fig 1 shows a spike or other feature in dv in a certain
    dr range. The user may want to investigate which stars
    contributed most heavily in that dr range. This section plots
    the positions of the stars (assumes x and y position data is
    available) colour coded according to how many times they appear
    in certain dr bins. The user can define which dr bins here
    as dr_start and dr_end. Because the stars have already
    been sorted into bins it is not possible to look at resolution
    below the bin width.

    Example: the user sees a spike in Fig. 1 between 2.23 and 2.74,
    and their bin width is 0.1 (parsecs assumed in this example,
    but irrelevant to the code). They set dr_start = 2.23 and
    dr_end = 2.74. It is not possible to take partial bins now the
    pairs have already been sorted and averaged. The closest match
    to the input values is found. With bin resolution of 0.1 that
    is 2.2 and 2.7, so colour coding will be according to pairs
    with dr between 2.2 and 2.7.
    '''

    # Get the range of bins the user's dr_star and dr_end translates to
    start_bin = min(range(len(edges)), key=lambda i: abs(edges[i]-dr_start))
    end_bin = min(range(len(edges)), key=lambda i: abs(edges[i]-dr_end))
    bin_range = range(start_bin, end_bin)

    # Check that the start and end points make sense
    if end_bin < start_bin:
        sys.exit('\n\nFig. 2 \n\nThe end bin is before or the same as the start bin: decrease dr_end')
    if start_bin < 0:
        sys.exit('\n\nFig. 2 \n\nYour dr_start is < 0. Its not possible for stars to have a negative distance between them.')

    # Count up how many times each star appears in any of the
    # bins in this range
    count_stars_range = [sum([count_stars_bins[a_bin][star] for a_bin in bin_range]) for star in range(n_stars)]

    # Make a new figure
    plt.figure(2)

    # Set a colourmap
    cm = plt.cm.get_cmap('rainbow')

    # Plot the colour coded stars
    plt.scatter(r[0], r[1], c=count_stars_range, alpha=0.8, cmap=cm, s=3.)

    # Add a colourbar
    plt.colorbar()

    # Add axis labels, etc
    plt.xlabel('x (pc)')
    plt.ylabel('y (pc)')
    plt.title('Counts in $\Delta r$ bins between ' + '%.1f' % dr_start + ' and ' + '%.1f' % dr_end + ' pc')

# Display the plots
plt.show()

