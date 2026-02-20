# UFV_MD_CoMDipCalc
A program that takes a .gro file for H20 and calculates the correlation function of the system. 


# How to use:
Compile the H20DipCorr.cpp file on your machine and then run the resulting executable with the arguments <input-filename> <out-filename> <shell-size-radius>. The input-filename is the name of the gro file you want the program to use, out-filename is the name you want the program to name the output file containing the correlation function data, and shell-size-radius is the radial width you want for each shell, which will separate the correlation function into shells. A shell width of 1.5 will result in all the particles who's distance from the center of mass of the system is less than 1.5 to be grouped, the next shell will have particles who's distance is between 1.5 and 3, the next shell after that will be for particles between 3 and 4.5, and this goes on until all particles for the system are included in a shell.

Example: If you have a gro file called 'dc.gro' and you want the program to output to a file called 'Output.txt' and you want each shell to have a radial width of 1.5, and you compile the program to an executable called 'H20DipCorr' you would run something similar to this (the specifics of how you compile or execute programs will of course depend on your OS):
H20DipCorr dc.gro Output.txt 1.5


# Understanding the output:
The output, if successful, the first few lines should look something like this:
"Timesteps: 40
Initial Time: 0
Final Time: 78360
Number of Shells: 4
Shell Width: 1.5
Shell: 0
0 2.58408
40 0.846912
80 0.540074
120 0.385863
160 0.305215
200 0.267207
240 0.195219
..."
The first line describes the size of the time steps between timestamps; then the next lines give the initial time and final timestamp of the correlation function(s). The 'Number of Shells' gives the total number of shells the system has (in this case, 4, resulting in 4 different correlation functions), and the 'Shell Width' describes the radial width of each shell. After these initial lines, the next line is 'Shell: 0,' indicating that the lines that follow are the correlation function values for the 1st shell. Each subsequent line gives a timestamp and the corresponding value of the correlation function at that timestamp, separated by whitespace. Here, the first line after 'Shell: 0' is '0 2.58408', indicating that the correlation function's value at t = 0 is 2.58408. The line '120 0.385863' says that at t = 120, the value is 0.385863. These lines will continue until the final timestamp/value for the correlation function. If you have more than 1 shell for your system, the line following the last correlation function value for 'Shell:0' will be 'Shell 1:' which will then be followed by subsequent lines that describe the correlation function for that particular shell. This process repeats until you have gone through all the shells in your system.

# A Note on Scaling:
This program does not scale any information given by the gro file. This was done for simplicity, as a constant factor can easily scale the analysis of the correlation function. For example, if the gro file lists a time as 't=  40.00000', we would treat this value as seconds. Additionally, for a line like '8963SOL      H35851  -1.799  -3.270   1.633 -0.1893 -0.6977  1.1569', we would take -1.799, 3.270, 1.633 as if they were in meters. What this equates to is that if the timesteps are physically represented picoseconds (10^-12) in the gro file, you would have to scale time values in the output file of this program by a factor of  10^-12. If the positions were in nanometers (10^-9), you would have to scale the correlation function by a factor of 10^(2*(-9 - 19)). Here, the -19 comes from the fact that we calculated the dipole moments using a charge value of 0.6732053036584803, whereas physically it would be 0.6732053036584803 x 10^-19. As the correlation function results from the dot product of the net dipole moments between timestamps, and each dipole moment in this example should be scaled by 10^(-9 -19), we have a total scaling of 10^(2*(-9 -19)), as previously stated.