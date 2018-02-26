# Velocity Structure Analysis Tool (VSAT)

Description
---

VSAT is a tool for analysing velocity structure in star clusters/associations 
following the method described in Arnold et al. (2018) 
(in prep). In short, for every possible pair of stars the distance
between them (dr) and their velocity difference (dv) is calculated.
Pairs are then sorted into dr bins. In each bin the average dv of
the pairs it contains is calculated. These average dvs are plotted
as a function of dr. 

How to run
---

To use this program the user must fill in the read_data function
in main.py to read in their data in the necessary format, which is
described in comments. The user must also set the error_flag to
True or False depending on whether or not the velocity data has
observational errors. After that the program can be run from the
command line like any other python script.

This program was written using python version 2.7.13, and may 
encounter errors if run with later versions.

License information
---

This code has  MIT license, reflecting that users are welcome to modify 
it to better suit their needs or make improvements. 

Contact information
---

For further information please contact Becky Arnold at
r.j.arnold.uk@gmail.com. 
