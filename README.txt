Velocity Structure Analysis Tool (VSAT)

VSAT is a tool for analysing velocity structure in star clusters/
associations following the method described in Arnold et al.(2018)
(submitted). In short, for every possible pair of stars the distance
between them (dr) and their velocity difference (dv) is calculated.
Pairs are then sorted into dr bins. In each bin the average dv of
the pairs it contains is calculated. These average dvs are plotted
as a function of dr. 

This program was written using python version 2.7.13.

Users are welcome to modify anything to better suit their needs or make
improvements.

To use this program the user must fill in the read_data function
in main.py to read in their data in the necessary format, which is
described in comments. The user must also set the error_flag to
True or False depending on whether or not the velocity data has
observational errors. After that the program can be run from the
command line like any other python script.

For further information please contact Becky Arnold at
r.j.arnold.uk@gmail.com. 
