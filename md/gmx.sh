gmx grompp -f ../input/min.mdp -c ../system/cellulose_IL.pdb -p ../system/cellulose_IL.top -o min.tpr -maxwarn 8
gmx mdrun -v -s min.tpr -deffnm min
gmx grompp -f ../input/constrain_md.mdp -c min.gro -p ../system/cellulose_IL.top -o constrain_md.tpr -maxwarn 8
gmx mdrun -v -s constrain_md.tpr -deffnm constrain_md
gmx grompp -f ../input/md.mdp -c constrain_md.gro -p ../system/cellulose_IL.top -o md.tpr -maxwarn 8
nohup gmx mdrun -v -s md.tpr -deffnm md -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu  &

