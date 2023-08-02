#!/bin/bash

# Create a topol.top file containing the forcefield and water model:
# #include "oplsaa.ff/forcefield.itp"
# #include "oplsaa.ff/spce.itp"
#
# [ System ]
# <name>
#
# [ Molecules ]

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
echo "Potential" | gmx energy -f min.edr -o min-energy.xvg
echo "Potential" | gmx energy -f min2.edr -o min2-energy.xvg
sed -i 's/@/#/g' min-energy.xvg
sed -i 's/@/#/g' min2-energy.xvg
gnuplot -e "set terminal png size 600, 450; set output 'min-energy.png'; plot 'min-energy.xvg' w l"
gnuplot -e "set terminal png size 600, 450; set output 'min2-energy.png'; plot 'min2-energy.xvg' w l"

# NVT Equilibration
# This step adds velocities so that a satisfactory temperature is achieved (298.15K)
gmx grompp -f mdp/nvt.mdp -o nvt.tpr -pp nvt.top -po nvt.mdp -c min2.gro
gmx mdrun -s nvt.tpr -o nvt.trr -x nvt.xtc -c nvt.gro -e nvt.edr -g nvt.log

echo "Temperature" | gmx energy -f nvt.edr -o nvt.xvg
sed -i 's/@/#/g' nvt.xvg
gnuplot -e "set terminal png size 600, 450; set output 'nvt.png'; plot 'nvt.xvg' w l"

# NPT Equilibration
# Changes the dimensions of the box so that the theoretical pressure (if the particles were to have collisions with the walls, instead of PBC) is equal to 1 bar.
gmx grompp -f mdp/npt.mdp -o npt.tpr -pp npt.top -po npt.mdp -c nvt.gro
gmx mdrun -s npt.tpr -o npt.trr -x npt.xtc -c npt.gro -e npt.edr -g npt.log

echo "Pressure" | gmx energy -f npt.edr -o npt.xvg
sed -i 's/@/#/g' npt.xvg
gnuplot -e "set terminal png size 600, 450; set output 'npt.png'; plot 'npt.xvg' w l"

# Production
gmx grompp -f mdp/prd.mdp -o prd.tpr -pp prd.top -po prd.mdp -c npt.gro
gmx mdrun -s prd.tpr -o prd.trr -x prd.xtc -c prd.gro -e prd.edr -g prd.log -v
