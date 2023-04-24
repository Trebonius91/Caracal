echo "This script performs rate constant calculations for "
echo "for the C-C bond dissociation with different applied "
echo "forces."

forces=(0.0 0.5 1.0 1.5 2.0 2.5 3.0)

for force in "${forces[@]}"
do
   mkdir ${force}nN
   echo "${force}nN"
   cp minimum.qmdff minimum.xyz rate.key ts.xyz ${force}nN/
   cd ${force}nN
   # replace force line with current force
   sed -i "s/ vec1 .*/ vec1 1 ${force}E-9 -6.3 2.5 0/g" rate.key
   sed -i "s/ vec2 .*/ vec2 8 ${force}E-9 6.3 -2.5 0/g" rate.key
   # start the calculation
   mpirun -np 6 ~/bin/calc_rate.x rate.key   > force_rate.log

   # copy the PMF profile back

   cp 300K_1bead/pmf_integration.dat ../pmf_${force}nN.dat
   cd .. 
done
