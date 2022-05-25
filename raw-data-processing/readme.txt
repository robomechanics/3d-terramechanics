This folder contains the data collected on the terramechanics testbed. Folders are specified by the gear ratio of the motor used and travel speed of the tests, and within each of these folders is a folder for each day of testing. 

Each day of testing has a bag file with ~10 seconds of force/torque data with the wheel lifted off the ground, called zero.bag. 

Each actual test takes the form XX_XX_XX_x.{txt,bag}, with the numbers being Vx, Vry, beta, and trial number.

To plot raw data: open plot_all_raw_data.m
In Matlab, open either all_smooth_data2.mat or all_grouser_data2.mat

Within plot_all_raw_data.m, scroll down to the plot you want to generate and Run This Section. Don't run the whole script, it will reprocess the data and take forever.
To generate individual plots of each trial, change the value of "disp_plots" from 0 to 1