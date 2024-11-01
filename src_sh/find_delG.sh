#!/bin/bash

mol_arr=('TPPD' 'TMED' 'DPPD' 'TDPA' 'TMED')
gas_phase=optim_gasphase.log
charge=('p2' 'p1' '0')
curr_dir=$PWD


#------------------Neutral molecule calculations---------------------------

exec 3>"AllResults.xlsx"
printf "%s\t %s\t %s\t %s\t %s\t %s\n" "molecule SCF_GP GCorr_GP SCF_SolP GCorr_SolP DelG" >&3
for molecule in "${mol_arr[@]}"; do

    echo "Gibbs calculations for $molecule"
    sol_phase=optim_nodiffuse_$molecule.log
    
    for qval in "${charge[@]}"; do

	if [ ! -d "$molecule/DME/q_$qval" ]; then
	   echo "ERROR!: $molecule/DME/q_$qval does exist."
	   continue
	fi

	cd $molecule/DME/q_$qval
    	printf "%s\t" ${molecule}_${qval} >&3

	#----------------------------Gas Phase Calculations----------------------------------
	if [ ! -f $gas_phase ]; then
    	   echo "$gas_phase not found!"
    	   printf "ERROR!: %s\t" "$gas_phase file not found" >&3
   	else
	   #---SCF - Gas-phase---
	   lineGSCF0=$(grep "SCF Done" "$gas_phase")
	   # Extract the SCF value
	   Gas_SCF0=$(echo "$lineGSCF0" | awk -F'=' '{print $2}' | awk '{print $1}' | tr -d ' ' | head -n 1)
	   # Print SCF0
	   echo "Gas-phase SCF0 for $qval: $Gas_SCF0"
	   printf "%s\t" "${Gas_SCF0}" >&3
	
	   #---Gibbs Correlation - Gas Phase---
	   lineGCorr0=$(grep "Gibbs Free Energy" $gas_phase)
	   # Extract the last occurrence Gibbs Correlation 
	   Gas_GCorr0=$(echo "$lineGCorr0" | grep "Thermal correction to Gibbs Free Energy=" | awk -F'=' '{print $2}' | tr -d ' ' | tail -n 1)
	   # Print GCorr0
	   echo "GCorr0 in gas phase for $qval: $Gas_GCorr0"
	   printf "%s\t" "${Gas_GCorr0}" >&3
	fi

	#----------------------------Solvent Phase Calculations----------------------------------
    	if [ ! -f $sol_phase ]; then
    	   echo "$sol_phase not found!"
    	   printf "ERROR!: %s\t" "$sol_phase file not found" >&3
        else
	   #---SCF - Solvent-phase---
	   lineSSCF0=$(grep "SCF Done" "$sol_phase")
	   # Extract the SCF value
	   Sol_SCF0=$(echo "$lineSSCF0" | awk -F'=' '{print $2}' | awk '{print $1}' | tr -d ' ' | head -n 1)
	   # Print SCF0
	   echo "Solvent-phase SCF0 for $qval: $Sol_SCF0"
	   printf "%s\t" "${Sol_SCF0}" >&3
	
	   #---Gibbs Correlation - Solvent Phase---
	   lineGCorr0=$(grep "Gibbs Free Energy" $sol_phase)
	   # Extract the last occurrence Gibbs Correlation 
	   Sol_GCorr0=$(echo "$lineGCorr0" | grep "Thermal correction to Gibbs Free Energy=" | awk -F'=' '{print $2}' | tr -d ' ' | tail -n 1)
	   # Print GCorr0
	   echo "GCorr0 in solvent phase for $qval: $Sol_GCorr0"
	   printf "%s\t" "${Sol_GCorr0}" >&3
        fi

	#---------------------------DelG_Solvation Calculation--------------------------------------
	if [[ -f $gas_phase && -f $sol_phase ]]; then
           #---Compute difference in Gibbs Free Energy---
	   DelG=$(awk "BEGIN {print $Sol_SCF0 - $Gas_SCF0}")
	   echo "DelG for ${molecule}_${qval}: $DelG"
	   printf "%s\n" "${DelG}" >&3
	else
	   echo "ERROR!: Incomplete data for DelG calculation"
	   printf "%s\n" "ERROR!: Incomplete data for DelG calculation"
	fi

	cd $curr_dir

    done
 
    printf "\n" >&3

done

exec 3>&-
