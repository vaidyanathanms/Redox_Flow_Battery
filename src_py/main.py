#----To copy NaPS structures from Avogadro output----------
#----To equilibrate and analyze NaPS structures------------
#----To analyze NMR spectra--------------------------------
#----Author: Vaidyanathan Sethuraman-----------------------
#----Date: Sept-01-2024------------------------------------
#----Requirements: my_gaussian_functions.py----------------

import numpy as np
import os
import shutil
import subprocess
import sys
import glob
from subprocess import call

#----my_gaussian_functions---------------------------------
from my_gaussian_functions import my_cpy_generic
from my_gaussian_functions import cpy_main_files
from my_gaussian_functions import find_inp_files
from my_gaussian_functions import edit_gen_inp_gauss_files
from my_gaussian_functions import run_gaussian
from my_gaussian_functions import clean_backup_initfiles

#----Input flags-------------------------------------------
num_hrs   = 8  # Total number of hours for run
num_nodes = 1  # Number of nodes
num_cores = 32 # Number of cores per node

#---------Input details------------------------------------
basis_fun    = 'UB3LYP/6-311+G(d,p)'
pop_style    = 'nbo' # nbo/reg
maxcycle     = 200 # Max # of optim steps 
maxstep      = 10  # Max size for an optim step = 0.01*maxstep Bohr
maxEstep     = 300 # maxEstep/1000 Bohr when moving from saddle pt
solv_arr     = ['Water','THF']#,'EthylEthanoate'] #water (default) THF, DiEthylEther
scrf         = 'pcm' # SCRF method pcm/smd/Dipole
multiplicity = 1 # Multiplicity is 1 unless it is a radical

# Add specific structures here or use 'all' keyword
# DO NOT PROVIDE the cml extension
spec_struct  = ['TMPD','TPPD','DPPD','TDPA']
                
#--------File lists--------------------------------------------
com_files = ['optim_var.com']
sh_files  = ['gauss_var.sh'] 

#---------Directory info---------------------------------------
maindir    = os.getcwd() #src_py dir
src_com    = '/home/vaidyams/all_codes/files_naps/src_com' #src_com dir
src_sh     = '/home/vaidyams/all_codes/files_naps/src_sh' #src_sh dir
guess_dir  = '/home/vaidyams/all_codes/files_naps/init_structs' #inp guess
scratchdir = '/projects/iontransport' #output dir
gauss_exe  = 'g16' # lmp executable file
scr_head   = 'naps_analysis_311dpdiff' # head dir for scratch outputs

#--------Find all/specific structure files-----------------------

inp_struct_list = find_inp_files(guess_dir, spec_struct)

#---------Main analysis------------------------------------------
for structs in inp_struct_list:

    workdir1 = scratchdir + '/' + scr_head   
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    inpfyle = os.path.basename(structs)
    inproot = inpfyle.split('.')[0]
    workdir2  = workdir1 + '/' + inproot

    if not os.path.isdir(workdir2):
        os.mkdir(workdir2)
	
    print("----Starting Gaussian calcs for " + inproot +\
          "-----")

    # Analyze for different solvents
    for solvent in solv_arr:

        workdir_main = workdir2 + '/' + solvent
        if not os.path.isdir(workdir_main):
            os.mkdir(workdir_main)

        print("Solvent: " + solvent)
            
        os.chdir(workdir_main)
        destdir = os.getcwd()
        print( "Current dir: ", destdir)
        
        #---Copying files------
        print( "Copying files...")
        
        cpy_main_files(guess_dir,destdir,inpfyle)
                
        for fyllist in range(len(com_files)):
            cpy_main_files(src_com,destdir,com_files[fyllist])
        
        for fyllist in range(len(sh_files)):
            cpy_main_files(src_sh,destdir,sh_files[fyllist])

        #---Write headers to input files------
        print("Generating input Gaussian file...")
        edit_gen_inp_gauss_files(com_files,inpfyle,basis_fun,\
                                 maxcycle,maxstep,maxEstep,scrf,\
                                 solvent,pop_style,multiplicity)


        #---Edit and submit gauss_var.sh----------------
        print("Editing and submitting submission script...")
        run_gaussian(sh_files[0],inproot,com_files[0],\
                     'job_gauss.sh',num_hrs,num_nodes,num_cores)

        #---Cleaning up jobs----------------------------
        clean_backup_initfiles(destdir,inpfyle)    
