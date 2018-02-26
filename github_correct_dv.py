import numpy as np
from scipy.stats.kde import gaussian_kde
import fort_subroutines

'''
Observational errors artificially increase the observed mean
dv because they broaden the observed velocity distribution.
This function estimates the true mean dv and adjusts the
observed mean dv values to reflect it.

It does this by assuming the true velocity distribution is
the same shape as the observed distribution, but narrower.
The width of the distribution is adjusted, velocities are
drawn and observational uncertainties applied, resulting in
simulated observed velocity distributions for each width. The
width which best reproduces the observed velocity distribution
is the best fitting model. A large number of velocities are
then drawn from the best fitting model and the mean dv is
calculated. This is the estimated value of the true dv.

This function assumes the uncertainty on each velocity
measurement is the same. It takes the number of stars, the
velocities, the uncertainty on each velocity measurement
and the mean dv values observed in each bin. It returns the
corrected mean dv list.
'''
def correct(n_stars, v, uncertainty, mean_dv):

    # Get the average dv of the observed velocities
    obs_av_dv = fort_subroutines.calc_av_dv(n_stars, v, len(v))

    # For each component of the observed velocities model
    # the true velocity disrtibution. Then draw a large number
    # of velocities from it. These will be used to estimate
    # the average dv of the true distribution. Create a list
    # to hold these drawn velocities
    v_true_sim = []

    # Each component of the velocity (x, y, z) is modelled
    # separately
    for v_component in v:

        # Smooth the observed velocity data into a probability
        # density function (pdf) using a gaussian kernel.
        # The pdf is probability as a function of velocity, hold
        # the velocities (pdf_v) and their corresponding probabilities
        # (pdf_prob) in lists
        kde = gaussian_kde(v_component)
        pdf_v = np.linspace(min(v_component) - 0.2, max(v_component) + 0.2, 100)
        pdf_prob = [i*95 for i in kde(pdf_v)]

        # Model the true velocity distribution by scaling the pdf to
        # be narrower by different amounts, drawing data from these
        # model distributions, applying the observational uncertainties
        # and identifying which best reproduces the observed velocity
        # distribution. Return the scaled velocity range of the best model
        pdf_v_opt = optimum_model(n_stars, v_component, uncertainty, pdf_v, pdf_prob)

        # Generate data drawn from this model velocity distribution
        v_true_sim.append(draw_velocities(pdf_v_opt, pdf_prob, 10000))

    # Find the average dv of the velocities drawn from the model
    # distributions. This is the estimated value of average dv of the
    # true velocities.
    est_true_av_dv = fort_subroutines.calc_av_dv(10000, v_true_sim, len(v_true_sim))

    # Correct the observed mean_dv list to remove inflation
    # due to observational uncertainties.
    correction_factor = obs_av_dv - est_true_av_dv
    mean_dv = [val - correction_factor for val in mean_dv]

    return mean_dv

'''
This function models the true velocity distribution and
returns the best fitting model.

It does this by assuming the true velocity distribution is
the same shape as the observed distribution, but narrower.
The width of the distribution is adjusted, velocities are
drawn and observational uncertainties applied, resulting in
simulated observed velocity distributions for each width. The
width which best reproduces the observed velocity distribution
is the best fitting model.

This function takes the number of stars, the velocities, the
uncertainty on each velocity measurement, and the pdf of the
observed velocities which is held as a list of velocities and
a corresponding list of probabilities. It returns the adjusted
velocity list which has been re-scaled to match the best fitting
model.
'''
def optimum_model(n_stars, v_component, uncertainty, pdf_v, pdf_prob):

    # Figure out the observed width of the velocity component's distribution
    vdist_obs_width = np.std(v_component)

    # Uncertainties broaden distributions so the width of the true distribution
    # will be smaller than the observed width, but likely greater than
    # observed width - uncertainty per measurement. Calculate the minimum
    # likely width, add a little padding to be safe, and make sure the padding
    # is not so much it makes it negative. If it does just set it to a very
    # small value (0.05).
    min_width = max(0.05,  (vdist_obs_width - uncertainty - 0.25))

    # Create a list of 100 possible widths of the true distribution to test
    # between the observed (inflated by uncertainties) width and the
    # minimum likely width.
    widths = np.arange(min_width, vdist_obs_width, (vdist_obs_width - min_width)/100.)

    # Test how well each width reproduces the observed velocity distribution
    # by calculating the area between cumulative distribution function (cdf) of
    # the observed velocities and simulated observations of velocities from each
    # model. If they're a good match the area should be low. Create a list to
    # hold the areas.
    area_list = []

    # Sort the observed velocities to aid making a cdf out of them.
    v_component.sort()

    # Go though each width
    for width in widths:

        # Scale the velocity distribution to the given width.
        ratio_adjust = width / vdist_obs_width
        pdf_v_width = [v_pos*ratio_adjust for v_pos in pdf_v]

        # Do 100 tests for each width, create a list to hold the area
        # between the observed and simulated cdfs for each test.
        test_areas = []
        for test in range(100):

            # Draw velocities from the distribution given this width,
            # and apply observational uncertainties
            v_width_data = draw_velocities(pdf_v_width, pdf_prob, n_stars)
            v_width_data = [np.random.normal(loc=vel, scale=uncertainty) for vel in v_width_data]

            # Sort the simulated velocities (needed to make cdf)
            v_width_data.sort()

            # Calculate the area between the observed and simulated cdfs
            area = 0.
            for datapoint in range(len(v_component)):
                area += abs(v_component[datapoint] - v_width_data[datapoint])
            area = area / len(v_component)

            test_areas.append(area)

        # Add the average area for this width to the area list
        area_list.append(np.mean(test_areas))

    # Find which width minimises the area between the cdfs (so best
    # reproduces the observations)
    optimum_width = widths[area_list.index(min(area_list))]

    # Adjust the observed pdf to fit the optimum width
    ratio_adjust = optimum_width / vdist_obs_width
    pdf_v = [v_pos*ratio_adjust for v_pos in pdf_v]

    return pdf_v


'''
This function draws velocities from an inputted distributions' pdf.
Does this by changing the pdf into a cdf which has a fixed range 0 - 1.
It then draws random numbers between 0 and 1 and uses the cdf to transform
them into velocities.
'''
def draw_velocities(pdf_v_width, pdf_prob, n_stars):

    # Change the pdf to a cdf
    frac = 1. / sum(pdf_prob)
    cdf_v = [sum(pdf_prob[0:i + 1])*frac for i in range(len(pdf_prob))]

    # Create an empty list to hold the drawn velocities
    drawn_vels = []

    # Draw the desired number of velocities
    for draw in range(n_stars):

        # Pick a random number between 0 and 1
        rand = np.random.random()

        # The cdf has a finite number of points between 0 and 1, so find the
        # points either side of the random number and their corresponding
        # velocities. We will then interpolate between these.
        for i in range(len(cdf_v)):

            if cdf_v[i] > rand:

                val_above_in_cumu = cdf_v[i]
                correspond_above_v = pdf_v_width[i]

                val_below_in_cumu = cdf_v[i - 1]
                correspond_below_v = pdf_v_width[i - 1]

                break

        # Interpolate between the two velocities to produce the drawn velocity.
        frac_ad = (rand - val_below_in_cumu)/(val_above_in_cumu - val_below_in_cumu)
        additional = frac_ad*(correspond_above_v - correspond_below_v)
        drawn_v = correspond_below_v + additional

        # Add the velocity to the list of drawn velocities
        drawn_vels.append(drawn_v)

    return drawn_vels
