# Version 1.01, Becky Arnold, February 2018 running on python version 2.7.13
import matplotlib.pyplot as plt
import parameters
import read_data
import calc_drdv_and_sort
import correct_dv
import plotting

'''
Set the parameters that will be needed to run the program.
In this function the user should set the path_to_data variable
to the path to the data file. They should also set error_flag to
True if there are errors on the velocity data and False otherwise.
'''
path_to_data, error_flag, dv_default, correct_inflation, bin_width, dr_start, dr_end = parameters.set_parameters()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


'''
Read in the data. The user must fill in this function in read_data.py.
n_stars is the number of stars, r is the position data of the stars,
v is the velocity data, and verr is the error on the velocity data.
If there are no errors on the velocities verr should be 0.
'''
n_stars, r, v, verr = read_data.read_data(path_to_data)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
Calculate dr and dv for every possible pair of stars.
Then sort the pairs into dr bins and calculate the mean dv in each
bin with errors.

This function returns the number of bins (n_bins), the edges of the
bins (edges), the mean dv in each bin (mean_dv), the uncertainty on each
mean (error), the number of pairs in each bin (n_in_bins), and a count
of how many times each star appears in each bin (count_stars_bins).
'''
n_bins, edges, mean_dv, error, n_in_bins, count_stars_bins = calc_drdv_and_sort.calc_drdv_and_sort(n_stars, r, v, verr, error_flag, dv_default, bin_width)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


'''
Observational errors artificially increase mean dv measured.
Correct this inflation of dv if there are errors and the user
has set the correct_inflation flag to True. Note that this
function assumes the uncertainty on each velocity measurement
is the same.
'''
if error_flag and correct_inflation:

    mean_dv = correct_dv.correct(n_stars, v, verr[0][0], mean_dv)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot results

# Plot the velocity structure (dr vs mean_dv)
plt.figure(1)
plotting.v_struct_plot(edges, mean_dv, error)

# Project the stars in 2D and colour code them according to how
# often they appear in bins between dr_star and dr_end (which are
# set in parameters.py).
plt.figure(2)
plotting.col_code_plot(n_stars, r, edges, count_stars_bins, dr_start, dr_end)

# Display the results
plt.show()
