import numpy as np


# Test 1 generates data with no velocity structure.
# The velocity structure shown in figure 1 should be roughly flat
# with increasing noise as dr increases.
def test_1():

    n_stars = 1000
    r = [[np.random.normal() for i in range(n_stars)] for j in range(3)]
    v = [[np.random.normal() for i in range(n_stars)] for j in range(3)]
    verr = [[np.random.uniform(0.5, 1.) for i in range(n_stars)] for j in range(3)]

    return n_stars, r, v, verr


# Test 2 generates data with linear velocity structure (distance between
# stars is linearly proportional to their velocity distance). The
# velocity structure shown in figure 1 should be linearly increasing.
def test_2():

    n_stars = 1000
    r = [[np.random.normal() for i in range(n_stars)] for j in range(3)]
    v = [[r[j][i] for i in range(n_stars)] for j in range(3)]
    verr = [[np.random.uniform(0.5, 1.) for i in range(n_stars)] for j in range(3)]

    return n_stars, r, v, verr

# Test 2 generates data with linear velocity structure (distance between
# stars is linearly proportional to their velocity distance). The
# velocity structure shown in figure 1 should be linearly increasing.
def test_3():

    n_stars = 1000
    r = [[np.random.normal() for i in range(n_stars)] for j in range(3)]
    v = [[r[j][i] for i in range(n_stars)] for j in range(3)]
    for coord in range(3):
        v[coord][-5] = v[coord][10]
        v[coord][-15] = v[coord][12]
        v[coord][-25] = v[coord][14]

    verr = [[np.random.uniform(0.5, 1.) for i in range(n_stars)] for j in range(3)]

    return n_stars, r, v, verr


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


'''
Read data function

This is the only function the user needs to write.
It must return the number of stars (n_stars), the position data (r),
velocity data (v) and (if there are any) velocity errors (verr) for
all stars. It is necessary for the user to write this function
because people format their data files differently.

The output of this function should be lists, each containing
lists of the form:

r = [[x positions of stars], [y positions of stars], [z positions of stars]]
v = [[vx velocities of stars], [vy velocities of stars], [vz velocities of stars]]
and (if errors)
verr = [[errors on vx velocities], [errors on vy velocities], [errors on vz velocities]]

If the full 3D information is not available then omit the relevant
lists. For example there may only be 2D positional data and 1D
radial velocity data available. In that case the output lists from
this function should be:

r = [[x positions of stars], [y positions of stars]]
v = [[radial velocities of stars]]
verr = [[measurement errors on radial velocities]]
If there are no measurement errors then return verr as 0.
'''


def read_data(path_to_data):

    # To test uncomment one and run main.py.
    # Ensure error_flag is set to True in parameters.py.
    # n_stars, r, v, verr = test_1()
    n_stars, r, v, verr = test_3()

    return n_stars, r, v, verr
