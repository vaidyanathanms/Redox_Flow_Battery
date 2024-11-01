# Generic Function Definitions
# Version_1: V_Sept_04_2024

import numpy as np
import os
import shutil
import subprocess
import sys
import glob
import re
import cclib
import xml.etree.ElementTree as ET # For reading Avogadro cml/xml files
from pathlib import Path

#---Generic function to copy files with *different* src/dest names
def my_cpy_generic(srcdir,destdir,inpfylname,destfylname):
    src_fyl  = srcdir  + '/' + inpfylname
    dest_fyl = destdir + '/' + destfylname
    shutil.copy2(src_fyl, dest_fyl)

#---Generic function to copy files with *same* src/dest names
def cpy_main_files(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl, desfyl)

#---Function to find all initital structures with cml extension
def find_inp_files(init_dir, inpstruct = 'all'):

    if not os.path.isdir(init_dir):
        raise RuntimeError(init_dir + " not found!")

    if isinstance(inpstruct,str):
        if inpstruct == 'all':
            inplist = glob.glob(init_dir + '/*.cml') 
            if inpstruct == []:
                raise RuntimeError('No cml input files found in ' + init_dir)
        else:
            if not os.path.exists(init_dir + '/' + inpstruct + '.cml'):
                raise RuntimeError(inpfyle + ' not found in ' + init_dir)
            inplist = [init_dir + '/' + inpstruct + '.cml']
    elif isinstance(inpstruct,list):
        inplist = []
        for structname in inpstruct:                
            fpath = init_dir + '/' + structname + '.cml'
            if not os.path.exists(fpath):
                print('ERROR', fpath, 'not found')
            else:
                inplist.append(fpath)
        if inplist == []:
             raise RuntimeError('No specified input cml found in '\
                                + init_dir)
    else:
        raise RuntimeError('Unknown type for spec_struct variable')
    
    return inplist

#---Edits headers of Gaussian i/p & adds coordinates to the files
def edit_gen_inp_gauss_files(com_files,inpfyle,basis_fun='6-31G**',\
                             maxcycle=100,maxstep=30,maxEstep=600,scrf='pcm',\
                             solvent='water',pop_style='reg',multiplicity=1):

    # Edit headers and options for calculations
    struct_name = inpfyle.split('.')[0] # Define structure name

    # Obtain total charge from Avogadro file or from filename
    tree = ET.parse(inpfyle)  # Read XML data from file
    root = tree.getroot()
    totcharge = root.get('formalCharge')
    if totcharge is None: totcharge=get_charges_from_fname(struct_name)
    print('Net_Charge: ', totcharge)

    # Edit com files
    for fyle in com_files:
        fr  = open(fyle,'r')
        fw  = open(struct_name + '_' + fyle.split('_')[0]+'.com','w')
        fid = fr.read().replace("py_basis",basis_fun)\
            .replace("py_maxcyc",str(maxcycle))\
            .replace("py_maxstep",str(maxstep))\
            .replace("py_maxEstep",str(maxEstep))\
            .replace("py_cavity",scrf)\
            .replace("py_solv",solvent)\
            .replace("py_popstyle",pop_style)\
            .replace("py_struct",struct_name)\
            .replace("py_charge",str(totcharge))\
            .replace("py_mult",str(multiplicity))
        fw.write(fid)
        fr.close()

        for atom in root.find('atomArray').findall('atom'):
            element_type = atom.get('elementType')
            x3 = atom.get('x3')
            y3 = atom.get('y3')
            z3 = atom.get('z3')
        
            # Write in the format: elementType x3 y3 z3
            fw.write(f"{element_type} {x3} {y3} {z3}\n")

        fw.write('\n') # Extra line at the end

#---Functiont to obtain charges from filename
#---Filename should be of the form *_qxm_*.cml where x is the charge
#---and m is optional for minus (negative) charges
def get_charges_from_fname(inpname):
    str_main = inpname.split('_')
    q_strarr = [[i,x] for i,x in enumerate(str_main) if x.startswith('q')]
    if len(q_strarr) == 0:
        raise RuntimeError('Input file format should have *_qxy_* in' \
                           'its filename and should repeat ONLY once!')
    for q_str in q_strarr:
        if q_str[0] == 0: continue # can be the naming
        qval = int(re.findall(r'\d+', q_str[1])[0])
        qval = -qval if 'm' in q_str[1] else qval
    return(qval)

#---Function to edit job submission files and submit jobs
def run_gaussian(inpjob,structname,ginput='optim_var.com',\
                 outjob='submit.sh',tot_hrs=1,tot_nodes=1,tot_cores=1):

    if not os.path.exists(inpjob):
        raise RuntimeError('ERROR: ' + inpjob + ' not found')

    ginp_main = structname + '_' + ginput.split('_')[0]
    if not os.path.exists(ginp_main+'.com'):
        raise RuntimeError('ERROR: ' + ginp_main + ' not found')
        
    jobstr = structname
    fr  = open(inpjob,'r')
    fw  = open(outjob,'w')
    fid = fr.read().replace("py_jname",jobstr).\
          replace("py_tottime",str(tot_hrs)).\
          replace("py_nnodes",str(tot_nodes)).\
          replace("py_ncores",str(tot_cores)).\
          replace("py_ncores",str(tot_cores)).\
          replace("py_nptot",str(tot_cores*tot_nodes)).\
          replace("py_ginput",str(ginp_main))
    fw.write(fid)
    fw.close()
    fr.close()
    subprocess.call(["sbatch", outjob])
    
#---Function to clean and backup files after generating i/p
def clean_backup_initfiles(destdir,structfile):
    
    initdir = destdir + '/init_files'
    if not os.path.isdir(initdir):
        os.mkdir(initdir)

    print(structfile)
    if os.path.exists(structfile):
        cpy_main_files(destdir,initdir,structfile)

    files = glob.glob('*var*')
    for fyl in files:
        if os.path.exists(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)
            
#---Function to find recent file 
def find_recent_file(destdir,keyword): 

    fylnames = destdir + '/' + keyword
    list_of_files = glob.glob(fylnames)
    if list_of_files != []:
        fyl_str = max(list_of_files, key=os.path.getctime)
        fyl_arr = fyl_str.split("/")
        print( "File Name: ", fyl_arr[len(fyl_arr)-1])
        return fyl_arr[len(fyl_arr)-1]
    else:
        return "nil"

#---Function to create files to write outputs and headers
def gen_output_files(outdir,nmr_elem='None', flag_nmr = 0,\
                     flag_freq = 0, flag_nbo = 0, flag_solv = 0):

    fid_nmr = 0; fid_freq = 0; fid_nbo = 0; fid_eall = 0; fid_eeqbm = 0
    fid_solv = 0
    
    fid_eall  = open(outdir + '/all_energy.dat','w')
    fid_eeqbm = open(outdir + '/all_eqbmenergy.dat','w')
    fid_eeqbm.write('%s\t%s\t%s\t%s\t%s\n' %('Structure','Solvent',\
                                             'SCF_Energy',\
                                             'Thermal_Correction',\
                                             'Total'))
    
    if flag_nmr:
        if nmr_elem == 'None':
            raise RuntimeError('ERROR: No reference for NMR analysis')
        fylename = set_filename(outdir, nmr_elem +'_nmr_output','csv')
        fid_nmr  = open(fylename,'w')
        fid_nmr.write('%s, %s, %s, %s, %s, %s, %s' %('Structure','NRef_Cntrs',\
                                                     'NMR_Freqs','NMR_FreqAvg',\
                                                     'NMR_Shifts','NMR_ShiftAvg',\
                                                      'ETot'))
        

    if flag_solv:
        fid_solv = open(outdir + '/all_solvation.dat','w')
        fid_solv.write('%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n'\
                       %('Structure','Solvent','Nelecs',\
                         'dGredox_gas','dGredox_sol',\
                         'dGSolv_neutral','dGSolv_redox','RedoxPot'))        

    if flag_freq:
        fid_freq = open(outdir + '/freq_analysis.dat','w')
    if flag_nbo:
        fid_nbo  = open(outdir + '/freq_analysis.dat','w')
    return fid_nmr, fid_freq, fid_nbo, fid_eall, fid_eeqbm, fid_solv

#---Function to set the filenames depending upon what is already there
# in the directory
def set_filename(dirname, fileroot, fmt='csv'):
    fnames = glob.glob(dirname + '/' + fileroot + '*')
    if fnames == []:
        return dirname + '/' + fileroot + '.' + fmt

    froots  = [Path(item).stem for item in fnames]
    root_id = [item.split('_')[-1] for item in froots if '_' in item]
    num_ids = [int(item) for item in root_id if item.isdigit()]
    
    if num_ids:
        return dirname + '/' + fileroot + '_' + str(max(num_ids)+1) +\
            '.' + fmt
    else:
        return dirname + '/' + fileroot + '_1.' + fmt

#---Main function to read and process the logfile
def analyze_from_logfile(log_file,flag_nmr=0,nmr_ref_elem='P'):
    # Open file and process each line
    efind = -1
    nmr_dum = []; scf_dum = []; Tcorr_dum = []
    with open(log_file,'r') as flog_id:
        for line in flog_id:
            line = line.strip()
            if flag_nmr and 'isotropy' in line and\
               line.split()[1] == nmr_ref_elem:
                nmr_dum.append(float(line.split()[4]))
            if 'SCF Done:' in line:
                E_eqbm = re.search(r"E\((U|R)B3LYP\)\s*=\s*(-?\d+\.\d+)",\
                                   line.strip())
                if E_eqbm:
                    scf_dum.append(float(E_eqbm.group(2)))
            if 'Thermal correction to Gibbs Free Energy=' in line:
                Tcorr_dum.append(line.split(' ')[-1])

    return nmr_dum, scf_dum, Tcorr_dum

#----Main function to analyze using cclib libraries
def analyze_with_cclib(log_file):
    data_parser   = cclib.io.ccopen(log_file)
    all_data      = data_parser.parse()
#    write_all_gendata(all_data)
    return all_data
    
#---Function to compute the nmr frequency of the reference solution
def compute_refzero_nmr(solvdir,log_file,nmr_ref_elem):
    # Open file and process each line
    value = 0; cnt = 0
    with open(log_file,'r') as flog_id:
        for line in flog_id:
            line = line.strip()
            if 'isotropy' in line and line.split()[1] == nmr_ref_elem:
                value += float(line.split()[4]); cnt += 1

    if cnt == 0:
        return -1
    return float(value/cnt)

#---Write NMR outputs to file
def write_nmr_outputs(fid_nmr,nmrvals,ref_nmrfreq,fshift,Geq):
    if len(nmrvals) == 0:
        fid_nmr.write(f'{len(nmrvals)} , Results not converged\n')
    else:
        fid_nmr.write(f'{len(nmrvals)} ,') #Avg
        fid_nmr.write(" , ".join([str(value) for value in nmrvals]))
        fid_nmr.write(f', {sum(nmrvals)/len(nmrvals)} ,') #Avg
        if fshift:
            fid_nmr.write(" , ".join([str(ref_nmrfreq - value) \
                                     for value in nmrvals]))
            fid_nmr.write(' , %g' %(np.sum(ref_nmrfreq - \
                                           np.array(nmrvals))/len(nmrvals)))

        fid_nmr.write(f', {float(Geq)}')#Etot
        fid_nmr.write('\n')
            
#---Write NMR outputs to file
def write_ener_outputs(fid_eall,fid_eeqbm,scf_arr,Tcorr_arr,Geq):
    if len(Tcorr_arr) == 0:
        [fid_eall.write(f'{float(en_val)}\n') for en_val in scf_arr] #All-E
        fid_eall.write(f'Results not converged\n')
        fid_eeqbm.write(f'Results not converged\n')
    else:
        [fid_eall.write(f'{float(en_val)}\n') for en_val in scf_arr] #All-E
        fid_eeqbm.write(f'{float(scf_arr[-1])} \t {float(Tcorr_arr[-1])} \t')
        fid_eeqbm.write(f'{float(scf_arr[-1])+float(Tcorr_arr[-1])}\n')#Etot

#---Find Redox type for calculations
# 0 - unox/unred (gas), 1 - ox/red (gas), 2 - unox/unred(solv)
# 3 - ox/red (solv)
def find_redox_type(phase,charge):
    idval = -1
    if 'sphase' == phase or 'sol' in phase:
        idval = 2 if charge.split('_')[-1].isnumeric() else 3
    if 'gphase' == phase or 'gas' in phase:
        idval = 0 if charge.split('_')[-1].isnumeric() else 1
    return idval
      
#---Write solvation/redox calculation outputs
def write_solv_outputs(Geq_arr,structname,solvent,nelectrons=1,\
                       fid_solv=0):

    fid_solv.write('%s\t%s\t%d\t' %(structname,solvent,nelectrons))

    # Print warnings
    if Geq_arr[0] == None:
        print('WARNING: neutral gas phase GFE not found\n')
    if Geq_arr[1] == None:
        print('WARNING: redox gas phase GFE not found\n')
    if Geq_arr[2] == None:
        print('WARNING: neutral sol phase GFE not found\n')
    if Geq_arr[3] == None:
        print('WARNING: redox sol phase GFE not found\n')

    # Write to outputs
    if Geq_arr[0] != None and Geq_arr[1] != None:
        delG_ox_gasp  = Geq_arr[1]-Geq_arr[0] #dG_redox gas phase
        fid_solv.write(f'{delG_ox_gasp}\t')
    else:
        fid_solv.write('NA\t')

    if Geq_arr[3] != None and Geq_arr[2] != None:
        delG_ox_solp = Geq_arr[3]-Geq_arr[2] #dG_redox sol phase
        fid_solv.write(f'{delG_ox_solp}\t')
    else:
        fid_solv.write('NA\t')

    if Geq_arr[2] != None and Geq_arr[0] != None:
        delG_solv_0  = Geq_arr[2]-Geq_arr[0]
        fid_solv.write(f'{delG_solv_0}\t') #dG_solv neutral species
    else:
        fid_solv.write('NA\t')
        
    if Geq_arr[2] != None and Geq_arr[0] != None:
        delG_solv_redox = Geq_arr[3]-Geq_arr[1]
        redox_pot = delG_solv_redox/nelectrons
        fid_solv.write(f'{delG_solv_redox}\t') #dG_solv redox species
        fid_solv.write(f'{redox_pot}\n') #Redox pot.
    else:
        fid_solv.write('NA\t')
        fid_solv.write('NA\n')   
        
#---Find all possible combinations of structural directories
def find_all_structs(headdir,spec_struct=['e']):
    if len(spec_struct) == 1:
        dirlist = []
        if spec_struct[0] == 'all':
            for dirpath in glob.glob(headdir + '/*/'):
                dirlist.append(os.path.basename(os.path.normpath(dirpath)))
        elif len(spec_struct[0].split('-')) == 0: #single letter
            for dirpath in glob.glob(headdir+'/'+spec_struct[0]+'*/'):
                dirlist.append(os.path.basename(os.path.normpath(dirpath)))
        elif len(spec_struct[0].split('-')) == 1: #range
            alpharange = [chr(i) for i in range(ord(spec_struct[0].split('-')[0])\
                                               ,ord(spec_struct[0].split('-')[1])+1)]
            for char in alpharange:
                allchardir = glob.glob(headdir + '/' + char + '*/')
                for dirpath in allchardir:
                    dirlist.append(os.path.basename(os.path.normpath(dirpath)))
        return dirlist
    else:
        return spec_struct
                
#---Function to check log files are present in the directory 
def find_logfile(destdir,basis_fun):
    log_file = glob.glob(destdir + '/*.log')
    if log_file == []:
        return 0
    elif len(log_file) == 1:
        val = is_basis_present(log_file[0],basis_fun)
        return log_file[0] if val else False
    else:
        val = is_basis_present(max(log_file,key=os.path.getctime),basis_fun)
        return max(log_file,key=os.path.getctime) if val else False

#---Check if the log file corresponds to the basis function of interest    
def is_basis_present(logfile,basis_fun):
    with open(logfile,'r') as fin:
        for line in fin:
            if line.strip().startswith("#"):
                return True
    return False
    
#---Close all output files
def close_all_outfiles(flag_nmr, fid_nmr, flag_freq, fid_freq, \
                       flag_nbo, fid_nbo, fid_eall, fid_eeqbm, \
                       flag_solv, fid_solv):
    
    fid_eall.close(); fid_eeqbm.close()
    if flag_nmr: fid_nmr.close()
    if flag_freq: fid_freq.close()
    if flag_nbo: fid_nbo.close()
    if flag_solv: fid_solv.close()
    

 
