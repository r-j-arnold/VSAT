# Velocity Structure Analysis Tool (VSAT)

Description
---

VSAT is a tool for analysing velocity structure in star clusters/associations 
following the method described in Arnold et al. (2018) 
(in prep). In short, for every possible pair of stars the distance
between them (dr) and their velocity difference (dv) is calculated.
Pairs are then sorted into dr bins. In each bin the mean dv of
the pairs it contains is calculated. These mean dvs are plotted
as a function of dr. This is the velocity structure. 

To help investigate the velocity structure a second plot is produced. It
shows a 2D projection of the stars. The stars are colour coded according to
how often they appear in user-specified dr bins. Understanding which stars
contribute most heavily to different parts of the velocity structure is
often helpful for interpreting it. This is described in more detail in
Arnold et al. (2018) (in prep).

This program is not dependent on the units of the data inputted,
although pc and km/s although are assumed and used to label figure
axes.

Structure of the program
---

1. Set parameters and flags
2. Read in the data
3. Calculate dr and dv for every pair of stars
4. Sort pairs into dr bins and calculate the mean dv in each bin with errors
5. Tidy up and get rid of any empty bins
6. Correct the inflation of mean dv due to uncertainties (if desired)
7. Make plots of the data- 1st plot: velocity structure
8. 2nd plot: positions of the stars with the stars colour coded

How to run
---

* In parameters.py set the path_to_data variable to the path to the data file.
* In parameters.py set the error_flag to True or False depending on whether
or not the velocity data has observational errors.
* In read_data.py fill in the read_data function to read in the data in the necessary format, which is described in comments.
* Run main.py from the command line like any other python script.

The user can set dr_start and dr_end in parameters.py to specify which
dr bins to colour code the second plot by.

This program was written using python version 2.7.13, and may 
encounter errors if run with later versions.

License information
---

VSAT has an MIT license, reflecting that it is free and users are
welcome to modify it to better suit their needs or make improvements,
but there is no warranty.

Contact information
---

For further information please contact Becky Arnold at
r.j.arnold.uk@gmail.com. 
