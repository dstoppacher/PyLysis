"""halo_startup.py

    Is the starup-file for the HaloAnalysis class
"""

# Load packages
from __future__ import print_function
import os
system_info = os.getcwd()
start = system_info.find('anaconda')
if start==-1:
    mycomp='/home/dstoppac/'
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
#task='walk_merger_trees'
#task='load_ascii_file'
#task='get_halo_info'
task='get_merger_tree_stats'
#task='create_DescScales_file'
#task='create_snapidzred_file'



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

#Cholla simulation configuration
res=3840

pathname_dict={       
            'path_basename':                '/data/'+str(res)+'_VSMDPL_160Mpc/',#Customize your main path where the data can be found (don't forget to but the halo files in a folder called 'halo_files')
            'path_orginal_data_hydro':      '',
            'output_filename_data':         'all_data_VSMDPL.txt'
            }

pathname_dict.update({'output_filename_merger_trees':   pathname_dict['path_basename']+'halo_files/trees/VSMDPL'+str(res)+'_160Mpc_tree_2_5_0_CT-v1.01_MB.txt'})
pathname_dict.update({'CT_merger_trees': pathname_dict['path_basename']+'halo_files/trees/tree_2_5_0.dat'})

if res==3840:

    config_dict={
                'first_snapid_in_sim': 33,  #first/lowest snapid, highest redshift in the simulation
                'total_snaps_sim':  151,   #total number of snapshot present in this simulation
                'n_files_snapshot': 16,     #each snapshot is split into n files 
                'snapid_CT_offset': 0,   #Snapshot number in the original simulation where tree builder sets it first snapid (usually this is set to 0 in tree builder) (use: get_main_branch)   
                }


config_dict.update({
                'res':              res,
                'box_size_h-1Mpc':  160,   #box size in units of h^-1Mpc
                'snapid':           33,    #a certain haloid you want to check or find (use: load_entire_box, get_halo_info, compare_individual_halos, construct_main_prog_tree_RS)
                'snapid_end':       150,   #Snapshot number of the last snapshot to be read in (use: load_entire_box, compare_individual_halos, construct_main_prog_tree_RS)
            
                'haloid':           0,    #find a certain haloid or compare this one to another (haloid2) (use: load_entire_box, get_halo_info, compare_individual_halos)
                'haloid2':          0,     #find or compare a certain haloid, this is the haloid2 where haloid1 is compared to (use: compare_halo_info)
                
                'rootID ':          0, #treeID of the first merger tree to be read in (use: walk_merger_trees)
                'rootID_end':       0, #treeID of the last merger tree to be read in (use: walk_merger_trees)
                'myprops2read':     ['haloid','x_pos', 'y_pos', 'z_pos', 'mhalo', 'descIndex', 'n_particles','rvir'] #Choose property names which should be processed: The standard set of properties are: 'haloid','x_pos', 'y_pos', 'z_pos', 'mhalo', 'descIndex', 'n_particles','rvir'
                })
  

config_dict.update({'total_snaps2read': config_dict['snapid_end']-config_dict['snapid']+1})

#initalize class and call MAIN!
myHalo_Analysis = HA.HaloAnalysis(pathname_dict,config_dict)
data = myHalo_Analysis.MAIN(task)

print('EOF - END OF FUN\n-------------------------')

class ReturnData2Jupyter:
    
    def __init__ (self):
        print('initialising ... ReturnData2Jypyter!')
        
    def return_data2jupyter(self):
        
        return data

           



             

 
