import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

# Cosmetic tweaks to the font style for the figure
mpl.rcParams['font.size'] = '18'
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Plot the velocity structure (dr vs mean_dv)
def v_struct_plot(edges, mean_dv, error):

    # Plot the velocity structure
    plt.errorbar(edges, mean_dv, yerr=error, capsize=3.)

    # Add axis labels, etc
    plt.xlabel('$\Delta r$ (pc)')
    plt.ylabel('$\Delta v~ $' + '(km s' + r'$^{-1}$)')
    plt.title('Velocity structure')

    # Neaten up the display
    plt.xlim([0, edges[-1]])

    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
This is a tool to help investigate and understand different
parts of the velocity structure shown v_struct_plot().

It projects the stars in 2D and colour codes them according to how
often they appear in bins between dr_star and dr_end (which are
set in parameters.py).

Say Fig 1 shows a spike or other feature in dv in a certain
dr range. The user may want to investigate which stars
contributed most heavily in that dr range. They can do this by
modifying dr_start and dr_end to contain the feature of interest.

It is not possible to look at resolution below the bin width.
'''
def col_code_plot(n_stars, r, edges, count_stars_bins, dr_start, dr_end):

    # Can only do the plot if there's at least 2D spatial information
    if len(r) > 1:

        # Get the range of bins the user's dr_start and dr_end translates to
        start_bin = min(range(len(edges)), key=lambda i: abs(edges[i]-dr_start))
        end_bin = min(range(len(edges)), key=lambda i: abs(edges[i]-dr_end))
        bin_range = range(start_bin, end_bin)

        # Check that the start and end points make sense
        if end_bin < start_bin:
            sys.exit('\n\nFig. 2 \n\nThe end bin is before or the same as the start bin: decrease dr_start or increase dr_end')
        if start_bin < 0:
            sys.exit('\n\nFig. 2 \n\nYour dr_start is < 0. Its not possible for stars to have a negative distance between them.')

        # Count up how many times each star appears in any of the
        # bins in this range
        count_stars_range = [sum([count_stars_bins[a_bin][star] for a_bin in bin_range]) for star in range(n_stars)]

        # Set a colourmap
        cm = plt.cm.get_cmap('viridis')

        # Plot the colour coded stars
        plt.scatter(r[0], r[1], c=count_stars_range, alpha=0.8, cmap=cm, s=3.)

        # Add a colourbar
        plt.colorbar()
        cbar.set_label('Counts', rotation=270, labelpad=20)

        # Add axis labels, etc
        plt.xlabel('x (pc)')
        plt.ylabel('y (pc)')
        plt.title('Counts in $\Delta r$ bins between ' + '%.1f' % dr_start + ' and ' + '%.1f' % dr_end + ' pc')

    return
