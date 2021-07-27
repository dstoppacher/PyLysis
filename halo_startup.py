"""halo_startup.py

    Is the starup-file for the HaloAnalysis class
"""

# Load packages
from __future__ import print_function
import os
system_info = os.getcwd()
start = system_info.find('anaconda')
if start==-1:
    mycomp='/z/doris/'
else:
    mycomp = system_info[:start]

import compileall
compileall.compile_dir(mycomp+'anaconda/pro/PyLysis/', force=1)

import haloAnalysis as HA

#Enter what you want to analyse!
#choose one of those:
 
task='load_entire_box' 
task='compare_individual_halos'
task='construct_main_prog_tree'
task='get_main_branch'
task='walk_merger_trees'
task='load_ascii_file'
#task='get_halo_info'
#task='get_merger_tree_stats'



##############################################################################################################
# User Configuration Space
   
"""README:
    -self.path_base_name looks for a folder 'halo_files' where halo catalogs according to the standard outputs
    such as ROCKSTAR or GADGET are stored
    -there is the need of file called 'snapidzred.txt', where the redshift and snapshot information are stored
    in such a way, that the following pandas function
    can read it.
    -the total number of snapshot which contain halos might be handy to know!
"""

pathname_dict={       
'path_basename':                '/data/256_hydro_Cholla50Mpc_halo_tests/',#Customize your main path where the data can be found (don't forget to but the halo files in a folder called 'halo_files')
#'path_basename':               '/home/dstoppac/anaconda/pro/data/Cholla256_50Mpc_halo_tests/',   
#'path_basename':               '/data/groups/comp-astro/bruno/cosmo_sims/256_hydro_50Mpc_halo_tests/', #Test files provided by Bruno (01/29/2021) here
'output_filename_merger_trees': 'Cholla256_50Mpc_CT-MB.txt',
'output_filename_data':         'all_data_CTreeID79509.txt',
'CT_merger_trees':              'Cholla256_50Mpc_halo_tests_trees_1.0.txt'
}

pathname_dict.update({'output_filename_data_tree_basename':         pathname_dict['path_basename']+'Cholla256_50Mpc_CT-MB/data_CT_treeID'})




config_dict={
'total_snaps_sim':  194,   #total number of snapshot present in this simulation
'n_files_snapshot': 8,     #each snapshot is split into n files 
'snapid':           65,    #a certain haloid you want to check or find (use: load_entire_box, get_halo_info, compare_individual_halos, construct_main_prog_tree_RS)
'snapid_end':       258,   #Snapshot number of the last snapshot to be read in (use: load_entire_box, compare_individual_halos, construct_main_prog_tree_RS)
'snapid_CT_offset': 104,   #Snapshot number in the original simulation where tree builder sets it first snapid (usually this is set to 0 in tree builder) (use: get_main_branch)

'haloid':           65,    #find a certain haloid or compare this one to another (haloid2) (use: load_entire_box, get_halo_info, compare_individual_halos)
'haloid2':          3,     #find or compare a certain haloid, this is the haloid2 where haloid1 is compared to (use: compare_halo_info)

'treeID ':          79622, #treeID of the first merger tree to be read in (use: walk_merger_trees)
'treeID_end':       79669, #treeID of the last merger tree to be read in (use: walk_merger_trees)
'myprops2read':     ['haloid','x_pos', 'y_pos', 'z_pos', 'mhalo', 'descIndex', 'n_particles','rvir'], #Choose property names which should be processed: The standard set of properties are: 'haloid','x_pos', 'y_pos', 'z_pos', 'mhalo', 'descIndex', 'n_particles','rvir'
}

config_dict.update({'total_snaps2read': config_dict['snapid_end']-config_dict['snapid']+1})

#initalize class and call MAIN!
myHalo_Analysis = HA.HaloAnalysis(pathname_dict,config_dict)
myHalo_Analysis.MAIN(task) 

                
print('EOF - END OF FUN\n-------------------------')                   

 
