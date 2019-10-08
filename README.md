# APT-Grain-Boundary-Analysis

Script for calculating Gibbsian interfacial excess values and compositions of interfaces.

R will need to be installed as will the following packages:
"tidyverse","PeriodicTable","InterpretMSSpectrum"

Takes a .pos and .rrng file.  The .pos file should contain an interface that is perpendicular to one of the "x", "y", or "z" directions.

Will then combine the .pos and .rrng files and, with user input, calculate: 
1) Where the interface starts and ends (based on the statistial distribution of elements in the regions adjacent to the interface)
2) Gibbsian interfacial excess values
3) Composition of the GB region

Will return:
1) A .pos file of the interface region
2) A summary table of GB
3) .txt file with information on location of original .pos file, the width of the interface
