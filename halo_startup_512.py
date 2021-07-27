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

#task='create_DescScales_file'
#task='create_snapidzred_file'
#task='compare_individual_halos'
#task='construct_main_prog_tree'
#task='get_halo_info'
#task='get_tree_statistics'
#task='load_entire_box' 
task='load_ascii_file'
#task='walk_merger_trees'
#task='utilise'
task='identify_merger_trees'



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
res=512

pathname_dict={       
            'path_basename':                '/data/'+str(res)+'_hydro_Cholla50Mpc/',#Customize your main path where the data can be found (don't forget to but the halo files in a folder called 'halo_files')
            #'path_basename':               mycomp+'/anaconda/pro/data/Cholla'+str(res)+'_50Mpc/',   
            #'path_basename':               ', #Test files provided by Bruno (01/29/2021) here
            'path_orginal_data_hydro':      '/data/groups/comp-astro/bruno/cosmo_sims/'+str(res)+'_hydro_50Mpc/output_files_pchw18/',
            'output_filename_data':         'all_data_CTreeID79509.txt'
            }

pathname_dict.update({'output_filename_merger_trees':   pathname_dict['path_basename']+'halo_files/trees/Cholla'+str(res)+'_50Mpc_CT-v1.01.txt'})
pathname_dict.update({'CT_merger_trees': pathname_dict['path_basename']+'halo_files/trees/Cholla'+str(res)+'_50Mpc_trees_1.01.dat'})


if res==256:

    config_dict={
                'first_snapid_in_sim': 0,  #first/lowest snapid, highest redshift in the simulation
                'total_snaps_sim':  194,   #total number of snapshot present in this simulation
                'n_files_snapshot': 8,     #each snapshot is split into n files 
                'snapid_CT_offset': 104,   #Snapshot number in the original simulation where tree builder sets it first snapid (usually this is set to 0 in tree builder) (use: get_main_branch)   
                }

elif res==512:

    config_dict={
                'first_snapid_in_sim': 0,  #first/lowest snapid, highest redshift in the simulation
                'total_snaps_sim':  170,   #total number of snapshot present in this simulation
                'n_files_snapshot': 16,     #each snapshot is split into n files 
                'snapid_CT_offset': 57,   #Snapshot number in the original simulation where tree builder sets it first snapid (usually this is set to 0 in tree builder) (use: get_main_branch)   
                }


config_dict.update({
                'res':              res,
                'box_size_h-1Mpc':  50,   #box size in units of h^-1Mpc
                'snapid':           0,    #a certain haloid you want to check or find (use: load_entire_box, get_halo_info, compare_individual_halos, construct_main_prog_tree_RS)
                'snapid_end':       169,   #Snapshot number of the last snapshot to be read in (use: load_entire_box, compare_individual_halos, construct_main_prog_tree_RS)
            
                'haloid':           65,    #find a certain haloid or compare this one to another (haloid2) (use: load_entire_box, get_halo_info, compare_individual_halos)
                'haloid2':          3,     #find or compare a certain haloid, this is the haloid2 where haloid1 is compared to (use: compare_halo_info)
                
                'rootID ':          79622, #treeID of the first merger tree to be read in (use: walk_merger_trees)
                'rootID_end':       79669, #treeID of the last merger tree to be read in (use: walk_merger_trees)
                'myprops2read':     ['haloid','x_pos', 'y_pos', 'z_pos', 'mhalo', 'descIndex', 'n_particles','rvir'] #Choose property names which should be processed: The standard set of properties are: 'haloid','x_pos', 'y_pos', 'z_pos', 'mhalo', 'descIndex', 'n_particles','rvir'
                })
  

config_dict.update({'total_snaps2read': config_dict['snapid_end']-config_dict['snapid']+1})
config_dict.update({'simulation_name': str(res)+'Cholla', 'halo_finder': 'RS-v0.99', 'tree_builder': 'CT-v1.01'})


#initalize class and call MAIN!
myHalo_Analysis = HA.HaloAnalysis(pathname_dict,config_dict)
data = myHalo_Analysis.MAIN(task)

print('EOF - END OF FUN\n-------------------------')

class ReturnData2Jupyter:
    
    def __init__ (self):
        print('initialising ... ReturnData2Jypyter!')
        
    def return_data2jupyter(self):
        
        return data

           



             

 
