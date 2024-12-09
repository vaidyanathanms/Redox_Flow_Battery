#----To analyze the log files from Gaussian outputs--------
#----Author: Vaidyanathan Sethuraman-----------------------
#----Date: Oct-28-2024-------------------------------------
#----Requirements: my_gaussian_functions.py----------------

import numpy as np
import os
import shutil
import subprocess
import sys
import glob
import cclib
from subprocess import call

#----my_gaussian_functions---------------------------------
from my_gaussian_functions import my_cpy_generic
from my_gaussian_functions import cpy_main_files
from my_gaussian_functions import find_logfile
from my_gaussian_functions import gen_output_files
from my_gaussian_functions import find_all_structs
from my_gaussian_functions import analyze_from_logfile
from my_gaussian_functions import analyze_with_cclib
from my_gaussian_functions import compute_refzero_nmr
from my_gaussian_functions import write_nmr_outputs
from my_gaussian_functions import write_ener_outputs
from my_gaussian_functions import find_redox_type
from my_gaussian_functions import write_solv_outputs
from my_gaussian_functions import close_all_outfiles

# Constants
hartoev = 27.21138505
faraday = 9.64853321233100184e4 # in coulumbs/mol

#----Input flags-------------------------------------------
flag_nmr   = 0  # Flag to compute NMR spectra
flag_shift = 0  # Flag to compute NMR shift
flag_freq  = 0  # Flag to process frequency
flag_nbo   = 0  # Flag to process NBO
flag_solv  = 1  # Flag to process DelG_Solv and Redox potential

nmr_ref_elem    = 'P' # Reference element for NMR spectra
nmr_refzero_dir = 'ref_H3PO4_q0' # Reference solution dir
basis_fun       = 'B3LYP/6-31+G(d,p)' #basis function used
nelectrons      = 1 # number of electrons transferred

#---------Input details------------------------------------
solv_arr     = ['DME']#,'THF']#, 'DiethylEther', 'EthylEthanoate']
# all, structs w/o .cml, single letter or initial to final letters
spec_struct  = ['DPPD']#,TDPA','TMPD','DPPD','TPPD','DME']
phase_arr    = ['gphase','sphase']
charge_arr   = ['Clq_p1','Clq_p2']

#---------Directory info---------------------------------------
maindir    = os.getcwd() #src_py dir
scratchdir = '/scratch/vaidyams/mgrfb_analysis' #output dir
scr_head   = 'w_Mg_diffuse' # head dir for scratch outputs

#---------Generate required output files-----------------------
workdir1 = scratchdir + '/' + scr_head   
if not os.path.isdir(workdir1):
    raise RuntimeError(workdir1 + " not found!")
else:
    if not os.path.isdir(workdir1 + '/all_results'):
        os.mkdir(workdir1 + '/all_results')


if flag_nmr:
    ref_zero_head = workdir1 + '/' + nmr_refzero_dir
    if not os.path.isdir(ref_zero_head):
        raise RuntimeError(ref_zero_head + " not found!")
    
 
fid_nmr,fid_freq,fid_nbo,fid_eall,fid_eeqbm,fid_solv = \
    gen_output_files(workdir1+'/all_results',nmr_ref_elem,\
                     flag_nmr,flag_freq,flag_nbo,flag_solv)

#---------Main analysis---------------------------------------

# Analyze for different solvents
for solvent in solv_arr:

    scf_eqbm = []
        
    if flag_nmr:
        ref_zerosolv_dir = ref_zero_head + '/' + solvent
        if not os.path.isdir(ref_zerosolv_dir):
            print('ERROR: Reference solution directory for ' + \
                  solvent + ' not found!')
            fid_nmr.write(solvent + ' ref directory not found!')
            if flag_shift: continue

        log_file = find_logfile(ref_zerosolv_dir,basis_fun)
        if not log_file:
            print('ERROR! No/bad ref log files found in ' +\
                  ref_zerosolv_dir)
            fid_nmr.write(solvent + ' No/bad ref log file found!')
            if flag_shift: continue
        
        ref_nmrfreq = compute_refzero_nmr(ref_zerosolv_dir,log_file,\
                                          nmr_ref_elem)
        if ref_nmrfreq == -1:
            print('Isotropy keyword found, but no values for: ' + \
                  nmr_ref_elem + ' in ' + ref_zerosolv_dir )
            if flag_shift: continue

        fid_nmr.write(f'\nNMR_Freq for reference in {solvent},'\
                      f'{ref_nmrfreq}\n')
            

    all_structs = find_all_structs(workdir1, spec_struct)
    for structname in all_structs:
        # Structure directory present
        workdir2  = workdir1 + '/' + structname
        if not os.path.isdir(workdir2):
            print('ERROR! ' + workdir2 + ' not found')
            continue
       
        # Check if solvent directory is present
        workdir_sol = workdir2 + '/' + solvent
        if not os.path.isdir(workdir_sol):
            print('ERROR! ' + workdir_sol + ' not found')
            continue

        if flag_solv: Geq_arr = [None]*4
        
        for charge in charge_arr:
            # Check if charge directory is present
            workdir_q = workdir_sol + '/' + charge
            if not os.path.isdir(workdir_q):
                print('ERROR! ' + workdir_q + ' not found')
                continue
                    
            for phase in phase_arr:
                # Check if phase/destination directory is present
                workdir_main = workdir_q + '/' + phase
                if not os.path.isdir(workdir_main):
                    print('ERROR! ' + workdir_main + ' not found')
                    continue

                print('---Starting analysis for ' + structname + ' in '\
                      +  solvent + ' ( ' + phase + ' )------')

                # Check if Gaussian log file is present
                log_file = find_logfile(workdir_main,basis_fun)
                if not log_file:
                    print('ERROR! No/bad log files found in '\
                          + workdir_main)
                    continue


                # Process file 
                nmrvals, scf_all, Tcorr_dum = \
                    analyze_from_logfile(log_file,flag_nmr,nmr_ref_elem)
                                         
                # Process using cclib
                all_data = analyze_with_cclib(log_file)

                # Write energy data
                fid_eall.write(f'{structname}\t{solvent}\n')
                fid_eeqbm.write(f'{structname}\t{solvent}\t')
                write_ener_outputs(fid_eall,fid_eeqbm,scf_all,Tcorr_dum,\
                                   all_data.freeenergy)

                # Write NMR data
                if flag_nmr:
                    fid_nmr.write(f'{structname}, ')
                    write_nmr_outputs(fid_nmr,nmrvals,ref_nmrfreq,\
                                      flag_shift,all_data.freeenergy)

#                if flag_mol: write_mol_outputs(fid_mol,\
#                                               all_data.moenergies)

                if flag_solv:
                    Gid = find_redox_type(phase,charge)
                    if Gid == -1:
                        print("Check folder path naming convention")
                        continue
                    Geq_arr[Gid] = all_data.freeenergy

#                if flag_cspa:
#                    fragcharges = compute_cspa(

            #---End of phase loop---

        #---End of charge loop---
        print(fid_solv)
        if flag_solv: write_solv_outputs(Geq_arr,structname,solvent,\
                                         nelectrons,fid_solv)

    #---End of struct loop---
    
#---End of solvent loop---


#---Close all opened files    
close_all_outfiles(flag_nmr,fid_nmr, flag_freq,fid_freq, \
                   flag_nbo, fid_nbo, fid_eall, fid_eeqbm,\
                   flag_solv, fid_solv)

print('Analysis completed :) ...')

