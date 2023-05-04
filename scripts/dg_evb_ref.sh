#
# This script generates the input for the calculation of 
# reference points for DG-EVB calculations
#


# The parameter for calling is the number x of atoms!

# ./dg_evb_ref.sh [x]

echo " "
echo "usage: dg_evb_ref.sh [number of atoms in structure]"
echo "struc_num.txt, struc.xyz and orca.inp/gauss.com must be present!"
echo "The orca.inp/gauss.com/run.inp files should include only the header"
echo "including charge and multiplicity information!\n"
echo " -------------------------------------------"

echo "Do you use Orca,Gaussian,EVB-QMDFF or CP2K for reference calculations?"
read -p "Insert O, G, E or C!  " method

#echo "Now, say how many command lines the present orca.inp file has!"
#read orca_lines
echo " "

# read in array with structure trajectory
IFS=$'\n' read -d '' -r -a lines < struc_num.txt

# request the user 
echo "Should the calculations be done directly? (y/n)"
read calc
echo " "

rm ref_struc.xyz

# das Array durchgehn und erst die Energien rausschreiben
# danach die Strukturen rausschreiben
j=0
for i in "${lines[@]}"
do
   j=`expr $j + 1`
   echo "Reference-point $j"
   
   von=0
   bis=0


#   zB 3. Struktur mit 4 Atomen:
#	bis = ($i)3 x (Atomanzahl+2)    = 18
#	von = bis - (Atomanzahl+1)  = 13
   zahl1=$(($1+2))
   zahl2=$(($1-1))
   bis=$(($i*$zahl1))
   von=$(($bis-$zahl2))

   mkdir $j

   if [ "$method" = "O" ]
   then 
      cat orca.inp >> $j/run.inp
   fi
   if [ "$method" = "G" ]
   then
      cat gauss.com >> $j/gauss.com
   fi
   if [ "$method" = "E" ]
   then
      cp qmdff.key $j/qmdff.key
   fi
   if [ "$method" = "C" ]
   then
      cp run.inp $j/run.inp
   fi

   if [ "$method" = "O" ]
   then
      sed -n ''$von','$bis' p' struc.xyz >> $j/run.inp
      echo "*" >> $j/ref.inp
   fi
   if [ "$method" = "G" ]
   then 
      sed -n ''$von','$bis' p' struc.xyz >> $j/gauss.com
      echo " " >> $j/gauss.com
   fi
   if [ "$method" = "E" ]
   then 
      sed -n ''1','2' p'  struc.xyz >> $j/ts.xyz
      sed -n ''$von','$bis' p' struc.xyz >> $j/ts.xyz
   fi
   if [ "$method" = "C" ]
   then
      sed -n ''1','2' p'  struc.xyz >> $j/run.xyz
      sed -n ''$von','$bis' p' struc.xyz >> $j/run.xyz
   fi

   if [ "$calc" = "y" ]
   then
      cd $j
      if [ "$method" = "O" ]
      then
         orca run.inp > run.out
      fi
      if [ "$method" = "G" ]
      then
         g09 gauss.com
      fi
      if [ "$method" = "E" ]
      then
         evb_qmdff.x qmdff.key
      fi
      if [ "$method" = "C" ]
      then
         # enter CP2K command...
      fi

      cd ..
   fi
done
echo "dg_evb_ref.sh exits normally..."
