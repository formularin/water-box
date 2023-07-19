#!/bin/bash

# We can use SPC instead of SPC/E for this part because the .gro file only contains coordinates, therefore any three-point water structure is okay.
gmx solvate -cs spc216 -o conf.gro -box 2.3 2.3 2.3 -p topol.top

# First minimization step
# These steps change the positions of waters so that there is a low potential energy (e.g. none are on top of each other)

# This takes the data from the .gro, .top, and .mdp files and converts it into something that can be used as input to a simulation
# It outputs its own commented .mdp file (same values as the one inputted)
# It outptus its own .top file which is also the same as the inputted one, but copies in the included files for the force field and water model
gmx grompp -f mdp/min.mdp -o min.tpr -pp min.top -po min.mdp

# This and subsequent commands start pick up from the previous step using the -c <step>.gro  option. 
# This is minimization, so for these steps the .gro files contain only position, no velocity
gmx mdrun -s min.tpr -o min.trr -x min.xtc -c min.gro -e min.edr -g min.log

# Second minimization step
gmx grompp -f mdp/min2.mdp -o min2.tpr -pp min2.top -po min2.mdp -c min.gro
gmx mdrun -s min2.tpr -o min2.trr -x min2.xtc -c min2.gro -e min2.edr -g min2.log

# These commands parse the energy log and output the potential energy as a function of time in such a way that it can be plotted by gnuplot
gmx energy -f min.edr -o min-energy.xvg
gmx energy -f min2.edr -o min2-energy.xvg
# Before running this command, go into the two XVG files and replace the @ symbols with #
gnuplot -e "set terminal png size 600, 450; set output 'min-energy.png'; plot 'min-energy.xvg' w l"
gnuplot -e "set terminal png size 600, 450; set output 'min2-energy.png'; plot 'min2-energy.xvg' w l"

# NVT Equilibration
# This step adds velocities so that a satisfactory temperature is achieved (298.15K)
gmx grompp -f mdp/eql.mdp -o eql.tpr -pp eql.top -po eql.mdp -c min2.gro
gmx mdrun -s eql.tpr -o eql.trr -x eql.xtc -c eql.gro -e eql.edr -g eql.log

# Save the temperature for gnuplot - replace the @'s again
gmx energy -f eql.edr -o eql-temp.xvg
gnuplot -e "set terminal png size 600, 450; set output 'eql-temp.png'; plot 'eql-temp.xvg' w l"

# NPT Equilibration
# Changes the dimensions of the box so that the theoretical pressure (if the particles were to have collisions with the walls, instead of PBC) is equal to 1 bar.
gmx grompp -f mdp/eql2.mdp -o eql2.tpr -pp eql2.top -po eql2.mdp -c eql.gro
gmx mdrun -s eql2.tpr -o eql2.trr -x eql2.xtc -c eql2.gro -e eql2.edr -g eql2.log
gmx energy -f eql2.edr -o eql-pressure.xvg
gnuplot -e "set terminal png size 600, 450; set output 'eql-pressure.png'; plot 'eql-pressure.xvg' w l"

# Production
gmx grompp -f mdp/prd.mdp -o prd.tpr -pp prd.top -po prd.mdp -c eql2.gro
gmx mdrun -s prd.tpr -o prd.trr -x prd.xtc -c prd.gro -e prd.edr -g prd.log
