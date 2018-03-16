import numpy as np
# Import the subroutines written in fortran
try:
    import fort_subroutines
except ImportError as failed_import:
    print('\n\n', failed_import)
    print ('Likely cause: running the code with python 3 or later.')
    print('Suggest running with python 2.7.\n\n')
    import sys
    sys.exit()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
For every possible pair of stars calculate the distance between them, dr
and the difference in their velocities, dv.
If there are errors on the velocities then propagate them to get the
errors on the dvs.

This function returns the number of pairs (n_pairs), and array with the
drs, corresponding dvs and, if there are errors, the errors on the dvs (dr_dv),
and the largest dr (max_dr).

Much of the work of this function is done by a fortran subroutine as
fortran is significantly faster than python at number crunching.
The time saving due to this is small for a single snapshot
however for a many snapshots (e.g. if analysing every timestep of a simulation)
then it becomes significant.

The package f2py was used to make the fortran subroutine compatible with python.
If for any reason the user wishes to make changes to the subroutine then
they will need to build their altered version with f2py.
'''
def calc_dr_dv(n_stars, r, v, verr, error_flag, dv_default, bin_width):

    # Get the number of possible pairs of stars
    n_pairs = ((n_stars*n_stars) - n_stars) / 2

    '''
    This subroutine returns an array containing the dr and corresponding dv
    for every pair of stars, and the error on dv (if there are velocity errors).
    If there are no velocity errors then that column is all zeros.
    It also returns the largest dr in the cluster so we know the biggest dr bin that
    will be needed.
    '''
    if error_flag:
        dr_dv, max_dr = fort_subroutines.calc_dr_dv(n_stars, n_pairs, bin_width, len(r), r, len(v), v, error_flag, dv_default, verr)

        # Remove any pairs with exacly 0 dr, dv, or error.
        # These pairs will be very rare, and they screw up the error calculation
        remove_list = []
        for pair in range(n_pairs):
            
            if (abs(dr_dv[pair][0]) < 0.00001) or (abs(dr_dv[pair][1]) < 0.00001) or (abs(dr_dv[pair][2]) < 0.00001):
                print 'duplict'
                remove_list.append(pair)
                
    else:
        dr_dv, max_dr = fort_subroutines.calc_dr_dv(n_stars, n_pairs, bin_width, len(r), r, len(v), v, error_flag, dv_default)

        # Remove any pairs with exacly dr or dv.
        # These pairs will be very rare, and they screw up the error calculation
        remove_list = []
        for pair in range(n_pairs):
            
            if (abs(dr_dv[pair][0]) < 0.00001) or (abs(dr_dv[pair][1]) < 0.00001):
                
                remove_list.append(pair)

    #Remove the zero pairs
    remove_list.reverse()
    for remove in remove_list:

        dr_dv = np.delete(dr_dv, (remove), axis=0)

    n_pairs = len(dr_dv)

    return n_pairs, dr_dv, max_dr

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
This function sorts the pairs of stars into dr bins. It then calculates the
mean dv in each bin.

It returns the edges of the bins (edges), the number of
bins (n_bins) the mean dv in each bin (mean_dv), the uncertainty on each mean
(error), the number of pairs in each bin (n_in_bins), and a count of how many
times each star appears in each bin (count_stars_bins).

As in the above section much of the work of this function is done by a fortran
subroutine as fortran is significantly faster than python at number crunching.
If the user wishes to make changes to this subroutine then they will need to
build their altered version with f2py.
'''
def sort_into_bins(n_stars, n_pairs, dr_dv, max_dr, bin_width, error_flag):

    # Define the bin edges, and define the number of bins
    edges = list(np.arange(0, max_dr, bin_width))
    n_bins = len(edges)

    # Sort the pairs of stars into dr bins then calculate the mean dv in each.
    if error_flag:
        mean_dv, error, n_in_bins, count_stars_bins = fort_subroutines.sort_errors(n_stars, n_pairs, n_bins, dr_dv, bin_width)
    else:
        mean_dv, error, n_in_bins, count_stars_bins = fort_subroutines.sort_no_errors(n_stars, n_pairs, n_bins, dr_dv, bin_width)

    return edges, n_bins, mean_dv, error, n_in_bins, count_stars_bins

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
Tidy up function, just makes the resulting plots neater.

High dr bins contain few pairs because only stars on the extreme
outskirts are far enough apart to populate them, and by definition
there are few stars on the outskirts of a cluster. If a bin contains
few pairs its mean dv is unreliable and often noisy. Therefore we get
rid of bins at high dr.

Also get rid of any empty bins.
'''
def tidy_up(n_bins, edges, mean_dv, error, n_in_bins, count_stars_bins):

    # Get rid of high dv bins.
    # Define high dr as the dr beyond which the average number of pairs per
    # bin is < 10. Find that cutoff.
    for i in range(n_bins):

        # Get the average number of pairs per bin for bins with dr >= this one
        av_numb_in_bins = sum(n_in_bins[i:])/float(n_bins - i)

        # It the average number of pairs per bin is < 10 then the bins
        # beyond this point are likely noisy and unreliable
        # Therefore want to cut any bins beyond this cutoff point
        if av_numb_in_bins < 10:

            cutoff_bin = i

            break

    # Get rid of the bins above the cutoff bin
    edges = edges[:cutoff_bin]
    mean_dv = mean_dv[:cutoff_bin]
    error = error[:cutoff_bin]
    n_in_bins = n_in_bins[:cutoff_bin]
    count_stars_bins = count_stars_bins[:cutoff_bin]

    # The number of bins has changed now the outermost have been cut.
    # Recalculate.
    n_bins = len(mean_dv)

    # - - - - - - - - - - - - - - - - -

    # Get rid of any empty bins.

    # Figure out which bins to remove
    # Reversing the order just makes it easier to delete them without
    # screwing up the indecies
    remove_list = [i for i in range(n_bins) if n_in_bins[i] == 0.]
    remove_list.reverse()

    # Remove the empty bins
    for empty in remove_list:

        edges = np.delete(edges, empty)
        mean_dv = np.delete(mean_dv, empty)
        error = np.delete(error, empty)
        n_in_bins = np.delete(n_in_bins, empty)
        count_stars_bins = np.delete(count_stars_bins, empty)

    # The number of bins has changed now the empty ones have been
    # deleted. Recalculate.
    n_bins = len(edges)

    return n_bins, edges, mean_dv, error, n_in_bins, count_stars_bins

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
This function calculates dr and dv for every possible pair of stars.
It then sorts the pairs into dr bins and calculates the mean dv in each
bin with errors. Finally it tidies up, which just produces neater plots.

This function returns the number of bins (n_bins), the edges of the
bins (edges), the mean dv in each bin (mean_dv), the uncertainty on each
mean (error), the number of pairs in each bin (n_in_bins), and a count
of how many times each star appears in each bin (count_stars_bins).
'''
def calc_drdv_and_sort(n_stars, r, v, verr, error_flag, dv_default, bin_width):

    n_pairs, dr_dv, max_dr = calc_dr_dv(n_stars, r, v, verr, error_flag, dv_default, bin_width)

    edges, n_bins, mean_dv, error, n_in_bins, count_stars_bins = sort_into_bins(n_stars, n_pairs, dr_dv, max_dr, bin_width, error_flag)

    n_bins, edges, mean_dv, error, n_in_bins, count_stars_bins = tidy_up(n_bins, edges, mean_dv, error, n_in_bins, count_stars_bins)

    return n_bins, edges, mean_dv, error, n_in_bins, count_stars_bins
