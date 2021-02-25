"""Halo, merger tree handling and analysis Last update 2021-02-24 by DS

"""
from __future__ import print_function
# Load packages
#own
import halo_analysis_lib as ha_lib
import halo_analysis_load_cats as ha_lc

#system
import numpy as np
import pandas as pd
import numpy.lib.recfunctions as rcfuncs
from time import time

#Enter what you want to analyse!
#choose one of those:
 
task='load_entire_box' 
task='compare_individual_halos'
task='construct_merger_trees'
#task='walk_merger_trees'
#task='get_halo_info'
#task='get_merger_tree_stats'

class HaloAnalysis:
    
    def __init__ (self):
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
        #Customize your main path where the data can be found (don't forget to but the halo files in a folder called 'halo_files')
        self.path_basename='/data/256_hydro_Cholla50Mpc_halo_tests/'
        #self.path_basename='/home/dstoppac/anaconda/pro/data/Cholla256_50Mpc_halo_tests/'
        #Test files provided by Bruno (01/29/2021) here
        #self.path_basename='/data/groups/comp-astro/bruno/cosmo_sims/256_hydro_50Mpc_halo_tests/'        

        self.snapid_array = ha_lib.df_to_sarray(pd.read_csv(self.path_basename+'snapidzred.txt', 
                                                   skiprows=2, 
                                                   names=['snapid', 'z', 'a', 't'], 
                                                   sep=' '))
        
        #print(self.snapid_array)

        #Simulation specific input!
        #total number of snapshot present in this simulation
        self.total_snaps_sim    = 194
        self.n_files_snapshot   = 8
     
        


        
        ##############################################################################################################
    
    
    def MAIN(self, task):   
        ##############################################################################################################
        # User Configuration Space       
        snapid                  = 65    #Snapshot number of the first snapshot to be read in
        snapid_end              = 258    #Snapshot number of the last snapshot to be read in
        
        
        treeID                  = 0     #treeID of the first merger tree to be read in
        treeID_end              = 1000   #treeID of the last merger tree to be read in
        
        
        self.total_snaps2read   = snapid_end-snapid
        #Choose property names which should be processed: The standard set of properties are: 'haloid','x_pos', 'y_pos', 'z_pos', 'mhalo', 'descIndex', 'n_particles','rvir'
        self.myprops2read       = ['haloid','x_pos', 'y_pos', 'z_pos', 'mhalo', 'descIndex', 'n_particles','rvir']
               
        def load_entire_box():
            """load all snapshots and particle information between the snapid and the snapid_end"""
            self.myDataHalo, self.particleIDs, self.halo_redshift_info = self.load_data_ROCKSTAR(self.path_basename,
                                       snapid,
                                       snapid_end+1)

        def get_halo_info():
            """get infos of a certain halos at a certain snapshot
                input:
                    haloid  id number of the halo you want the info
                    snapid  snapshot number where the halo was found with haloid
                    
                output:
                    info of halo properties from self.myprops2read printed on the screen
            """
            load_entire_box()
            self.get_halo_info(haloid=3, snapid=snapid)
            

        def compare_individual_halos(): 
            """compare properties and print percentage of particle they share (between snapshots or at one snapshot)
                input:
                    haloid1     id of the first halo you want to compare
                    haloid2     id of the second halo you want to compare
                    snapid1     snapshot number of the first halo you want to compare
                    snapid2     snapshot number of the second halo you want to compare
                    
                keywords:
                    print_on_screen  default: False
                                        info of halo properties from self.myprops2read printed on the screen
                                        difference between the halo properties and particle content are calculate
                                        as absolute vales or as percentage where the properties of haloid2 (lower redshif) are compared
                                        to those of haloid1 (higher redshift)                   
                    
                output:
                    structured array of one column which contains the information                  

                    
            """
            load_entire_box()
            self.compare_halo_info(haloid1=3, haloid2=3, snapid1=snapid, snapid2=snapid+1, print_on_screen=True)        

        def construct_merger_trees():
            start = time() 
    
            data, merger_dict, first_z_dict = self.construct_merger_tree(n_trees=[],
                                                                           snapid=snapid,
                                                                           snapid_end=snapid_end)
            

            print('Time until all halos tracked:', format((time()-start)/60.0/60.0, '0.2f'), 'h /', format((time()-start)/60.0, '0.2f'), 'min /', format((time()-start), '0.1f'), 'sec')

            ha_lib.writeIntoFile(
                       self.path_basename+'all_data_MM.txt',
                       data[['haloid','descIndex', 'rootIndex', 'predIndex', 'mhalo', 'delta_mhalo', 'rvir', 'delta_rvir', 'x_pos', 'y_pos', 'z_pos', 'snapid', 'n_particles', 'npros']],
                       myheader='All halos fro Cholla 50Mpc 256 box from Bruno run 29/01/2021\n'\
                       +'(1) haloid (2) descIndex (3) rootIndex (4) predIndex (5) Mvir [h-1Msun] (6) delta_Mvir [h-1Msun] (7) Rvir [h-1kpc] (8) delt_Rvir [h-1kpc]'\
                       +'(9) X [h-1Mpc] (10) Y [h-1Mpc] (11) Z [h-1Mpc] (12) snapid (13) n_particles (14) n_progs',
                       mydelimiter='\t',
                       data_format='%i\t%i\t%i\t%i\t%0.6e\t%0.6e\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%i\t%i\t%i')


        def walk_merger_trees():
            
            load_entire_box()
            mytree_data = ha_lib.read_unshaped_txt(self.path_basename+'merger_trees.txt', 102, self.total_snaps_sim, delimiter=' ', comments='#', dtype=np.int32) 
            count_error=0
            failed_IDs=[]
            for i in range(0,100,1):
                if self.total_snaps2read>=mytree_data[i,1]-mytree_data[i,0]:
                    end_snap=mytree_data[i,1]-mytree_data[i,0]
                else:
                    end_snap=self.total_snaps2read

                print('i:', i, 'snapid_first:', mytree_data[i,0], 'snapid_last:', mytree_data[i,1], 'n_snaps2read:', end_snap, 'rootIndex:', mytree_data[i,2], 'progIDs:', mytree_data[i, 3:end_snap+3], 'check_sum:', len(mytree_data[i, 3:end_snap+3]))

                #exit()
                tree_structure=self.track_merger_tree(mytree_data[i,0], mytree_data[i,0]+end_snap, mytree_data[i,2], mytree_data[i, 3:end_snap+3])
                if tree_structure==False:
                    count_error+=1
                    failed_IDs.extend([mytree_data[i,2]])
                print('--> completed')

            print('\n-------------------\nSUMMARY:\n', i-count_error+1, 'trees sucessfully walked / ', count_error, 'failed!', '\nfailed treeIDs are:\n', failed_IDs)

        def get_merger_tree_stats(print_into_file=True):

            tree_data=self.load_data_merger_trees(treeID,treeID_end,file_basename='data_comparision_MM')
            #print(np.info(tree_data))
            merger_tree_stats, myheader, mycols, myformat_string, mydt = ha_lib.create_data_array('stats_perc_bucket')
            
            dt_basic = ha_lc.get_dtype_sarray('stats_basic_bucket')
            #print(dt_basic)                  
            myprops=ha_lib.get_props_for_stats_calc(dt_basic, input_is_dtype=True)       
            #print(myprops)
            
            for i, snap in enumerate(range(snapid,snapid_end,1)):
                
                data_stats = self.calc_statistics(tree_data, myprops, sel_name='snapid1', sel_value=snap, dt=mydt)                

                data_stats['snapid'] = snap
                data_stats['z']      = self.snapid_array[np.where(self.snapid_array['snapid']==snap)[:][0]]['z']
                
                if i==0:
                    merger_tree_stats = data_stats
                else:                                     
                    merger_tree_stats = np.append(merger_tree_stats, data_stats, axis=0)
    
                        
            #print(merger_tree_stats)
            #print(np.info(merger_tree_stats))
            if print_into_file==True:
                ha_lib.writeIntoFile(
                       self.path_basename+'merger_tree_stats.txt',
                       merger_tree_stats[mycols],
                       myheader='statistics of properties found in halos from Cholla 50Mpc 256 box from Bruno run 29/01/2021\n'+myheader,
                       mydelimiter='\t',
                       data_format=myformat_string)
            


        def caseSwitcher(task):
            choose = {
                'load_entire_box':          load_entire_box,
                'compare_individual_halos': compare_individual_halos,
                'construct_merger_trees':   construct_merger_trees,
                'walk_merger_trees':        walk_merger_trees,
                'get_halo_info':            get_halo_info,
                'get_merger_tree_stats':    get_merger_tree_stats                
                }
        
            func = choose.get(task)
            
            return func()
    
        caseSwitcher(task)
    
        
        ##############################################################################################################
    
    def assign_progenitor_indices(self,
                                  data_RS_now,
                                  data_RS_before,
                                  snapid_before,
                                  first_assignment=1,
                                  document_progress=False):
        """Function assigns predestant indices to the descdent indices and the haloids from a list of halo properties from ROCKSTAR out_xxx.list halo catalog. xxx stands for the 
        the snapshot number. The function goes through every snapshot and connect the haloid and descdent index of the main progenitor (also most massive progenitor) of a halo with
        the predestant index a snapshot before. It further assigns a unique tree index (treeID -- an ascending number of the first occurence of the halo)
        
        input:
        ==========
            data_RS_now:          (array) structured array of the ROCKSTAR halo properties catalog of the current snapshot
            data_RS_before:       (array) structured array of the ROCKSTAR halo properties catalog of the snapshot before            
            snapid_before_before: (float) redshift of one snapshot before, if the first snapshot is processed it will take this redshift
            
            keywords:
                first_assignment    default 1, then the structured array are expendent to host further properties
                docuemtn_progress   default False, if True prints progess and output data on screen
            
        output:
        ==========        
            data_RS_now:    (array) catalog of main progenitor halos with assigned predestant (predIndex) at the current snapshot, unique first progenitor ID (rootIndex), number of progenitors where the main progenitor
                            is chosen from, difference in halo mass (delta_mhalo) and virial radius (delta_rvir) to the halo on the tree
                            one snapshot before
            data_RS_before: (array) same as data_RS_now
        
        """

        data_RS_now = rcfuncs.append_fields([data_RS_now], 
                                            ['delta_mhalo','delta_rvir','predIndex','rootIndex','npros'], 
                                            [np.zeros(data_RS_now.size,),np.zeros(data_RS_now.size,),np.zeros(data_RS_now.size,),np.zeros(data_RS_now.size,),np.zeros(data_RS_now.size,)], 
                                            dtypes=['f4','f4','i8','i8','i4'], usemask=False)
       
        data_RS_now['npros']    = -99
        data_RS_now['rootIndex']= -99
        data_RS_now['predIndex']= -99

        if first_assignment==0:        
            data_RS_before = rcfuncs.append_fields([data_RS_before], 
                                                   ['npros','rootIndex','predIndex'] , 
                                                   [np.zeros(data_RS_before.size,),np.zeros(data_RS_before.size,),np.zeros(data_RS_before.size,)], 
                                                   dtypes=['i4','i8','i8'], usemask=False)
            data_RS_before['npros']    = -99
            data_RS_before['rootIndex']= -99
            data_RS_before['predIndex']= data_RS_before['haloid']        
    
        #Sort the data you want to assign a predIndex, the halo on top is the main progenitor (also the most massive one)    
        data_RS_before[::-1].sort(order=['descIndex','mhalo','haloid'], axis=0)
       
        #Look for the unique IDs, the functions selects the ID which appreas at first which are also our main progenitors. Store in count how often a certain descIndex is occuring
        #that is the number of progenitors of the current halos. Assign the satellite status as orphan==0 which means it is a central halo (note that this does not matter right now for
        #this analysis)
        non_unique, index, count = np.unique(data_RS_before['descIndex'], return_index=True, return_counts=True)
        data_RS_before['npros'][index]=count
        
        #Intersect now the haloids from this snapshot with the descdent idexes (descIndex) a snapshot before. Store in index_now and index_before where those halos
        #can be found in our structured array
        parentIndices, index_now, index_before = np.intersect1d(data_RS_now['haloid'], data_RS_before['descIndex'], return_indices=True)
    
        #The predestant index (predIndex) is the haloid of the halo on the main progenitor one snapshot before. We set this haloid as predIndex in the current snapshot       
        data_RS_now['predIndex'][index_now]=data_RS_before['haloid'][index_before]
        
        #We calculate the difference in halo mass/virial radius of the same halo found now and one snaphsot before
        data_RS_now['delta_mhalo'][index_now]=data_RS_now['mhalo'][index_now]-data_RS_before['mhalo'][index_before]
        data_RS_now['delta_rvir'][index_now]=data_RS_now['rvir'][index_now]-data_RS_before['rvir'][index_before]
        
        #We transfer further data from the halo before to the current halos e.g. the rootIndex which also serves as unique identification of the merger tree
        data_RS_now['npros'][index_now]=data_RS_before['npros'][index_before]
        data_RS_now['rootIndex'][index_now]=data_RS_before['rootIndex'][index_before]
    
        #and sort them    
        data_RS_now[::-1].sort(order=['predIndex','mhalo','haloid'], axis=0)
        data_RS_now[::-1].sort(order=['descIndex','mhalo','haloid'], axis=0)
    
        #We set rootIndex which is also the root of the merger tree (treeID) and a unique identification of the whole tree!
        #This ID is an ascendent count of all main progenitor trees in the data, if a new halo was find, then the rootIndex is assigned at first -99
        #those IDs need to be assigned now. We find them as data_wo_root (halos in data without a root ID)      
        data_wo_root=data_RS_now[np.where(data_RS_now['rootIndex']==-99)[:][0]]
    
        #We sort them by the haloid which is also ascending      
        data_wo_root.sort(order=['haloid'], axis=0)
        
        #We assign the rootIndex of the newly found halos an ascending number count. That means that newly found halos with lower haloids get lower rootIndex
        data_wo_root['rootIndex']=range(data_RS_now.size-data_wo_root.size,data_RS_now.size,1)
    
        #The last step is to assign now the right halos their new rootIndex (unique identification of the main progrenitor merger tree) by comparing the haloids       
        parentIndices, index_root, index_now = np.intersect1d(data_wo_root['haloid'], data_RS_now['haloid'], return_indices=True)
        data_RS_now['rootIndex'][index_now]=data_wo_root['rootIndex'][index_root]

        if document_progress==True:        
            print('data_RS_before:\n', data_RS_before[['haloid','descIndex', 'rootIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos', 'predIndex']],\
                  '\ndata_RS_now:\n', data_RS_now[['haloid','descIndex', 'rootIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos', 'predIndex']],\
                  'PROGENITOR INFO successfully assigned!\n#################################################')           

        return data_RS_now, data_RS_before
    
    
    def calc_delta_perc(self, prop1, prop2):
        """calculates the decrease or increase of prop2 in comparison with prop1 in percent"""
        
        return 100.0/prop1*prop2-100.0

    def calc_statistics(self,
                        data,
                        props,
                        sel_name='',
                        sel_value=None,
                        dt=None,
                        print_on_screen=False):
        """Function calculates the median and 10th, 32the, 68th and 90th percentiles. The median and 10th and 90th are also printed on the screen.
            A data-bucket as an 1-row structrued array is returned where the statistics of all properties defined in 'props' are calculated.
            Optionally one can select a subsample from a specific property defined as 'sel_name' and the corresponding value 'sel_value' e.g. snapshot
            column name 'snapid' and snsnapshot number '258' as value.            
        """
        #print(np.info(data))
        try:
            if sel_name!=None:
                try:
                    sample=data[np.where(data[sel_name]==sel_value)[:][0]]
                except:
                    sample=data
                    print('Selection of property', sel_name, 'not possible!')
            else:
                sample=data
                
            if print_on_screen==True:
                print('\nsel_name:', sel_name, 'sel_value:', sel_value, 'ngal:', sample.size, '\t', "{0:.0f}".format(100.0/data.size*sample.size),'% of all data was used to calculate statistics!')
                print('////////////////////////////////////////////////')
    
            #print('data selection:\n', sample[['snapid1', 'snapid2','haloid1','haloid2', 'descIndex1', 'descIndex2', 'rootIndex', 'n_particles1', 'n_particles2']])
            #print(np.info(sample))
            bucket_stats = np.zeros((1,), dtype=dt)
            #print(np.info(bucket_stats))
            
            bucket_stats['n_count'] = sample.size                  
            if print_on_screen==True:
                print(('prop').ljust(25), 'median 10th / 90th\t[unit]\n-----------------------------------------------------------------\n')
            for prop in props:
                unit, myformat = ha_lib.find_unit(prop)
                if print_on_screen==True:
                    print((prop).ljust(25), end='')
                if prop.find('perc')==-1 and prop.find('n_')==-1 and prop.find('pos')==-1 and prop.find('rvir')==-1:
                    if print_on_screen==True:
                        print("{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))), '\t',\
                              "{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))-np.log10(np.nanpercentile(sample[prop], 10))),\
                              '/', "{0:.2f}".format(np.log10(np.nanpercentile(sample[prop], 90))-np.log10(np.nanmedian(sample[prop]))), 'log('+unit+')')
                    
                    bucket_stats[prop+'_log50'] = "{0:.2f}".format(np.log10(np.nanmedian(sample[prop])))
                    bucket_stats[prop+'_log10'] = "{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))-np.log10(np.nanpercentile(sample[prop], 10)))
                    bucket_stats[prop+'_log32'] = "{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))-np.log10(np.nanpercentile(sample[prop], 32)))
                    bucket_stats[prop+'_log68'] = "{0:.2f}".format(np.log10(np.nanpercentile(sample[prop], 68))-np.log10(np.nanmedian(sample[prop])))
                    bucket_stats[prop+'_log90'] = "{0:.2f}".format(np.log10(np.nanpercentile(sample[prop], 90))-np.log10(np.nanmedian(sample[prop])))                            
                else:
                    if print_on_screen==True:
                        print("{0:.2f}".format(np.nanmedian(sample[prop])), '\t',\
                              "{0:.2f}".format(np.nanmedian(sample[prop])-np.nanpercentile(sample[prop], 10)), '/',\
                              "{0:.2f}".format(np.nanpercentile(sample[prop], 90)-np.nanmedian(sample[prop])), unit)
    
                    bucket_stats[prop+'_50'] = np.nanmedian(sample[prop])
                    bucket_stats[prop+'_10'] = np.nanmedian(sample[prop])-np.nanpercentile(sample[prop], 10)
                    bucket_stats[prop+'_32'] = np.nanmedian(sample[prop])-np.nanpercentile(sample[prop], 32)
                    bucket_stats[prop+'_68'] = np.nanpercentile(sample[prop], 68)-np.nanmedian(sample[prop])
                    bucket_stats[prop+'_90'] = np.nanpercentile(sample[prop], 90)-np.nanmedian(sample[prop])

        except:
            print('data selection not available! --> SKIPPED!')

        return bucket_stats                        
        


    def compare_halo_info(self,
                           props1=[],
                           props2=[],
                           haloid1=0,
                           haloid2=0,
                           snapid1=0,
                           snapid2=0,
                           print_into_file=False,
                           print_on_screen=False,
                           dt=None):
        """Function compares two halos and writes and returns their properties into a bucket data structed array"""

        print('\nCOMPARE ID1/ID2:', haloid1, '/', haloid2, 'at snapid1/snapid2:', snapid1, '/', snapid2, '\n#######################################################################\n')           
        
        if len(props1)==0:
            props1   = self.myDataHalo[np.where((self.myDataHalo['haloid']==haloid1) & (self.myDataHalo['snapid']==snapid1))[:][0]]
        if len(props2)==0:
            props2   = self.myDataHalo[np.where((self.myDataHalo['haloid']==haloid2) & (self.myDataHalo['snapid']==snapid2))[:][0]]
        
        z1           = self.snapid_array['z'][np.where(self.snapid_array['snapid']==snapid1)][0]
        z2           = self.snapid_array['z'][np.where(self.snapid_array['snapid']==snapid2)][0]
        n_particles1 = self.particleIDs[snapid1][haloid1]['particle_IDs'].size
        n_particles2 = self.particleIDs[snapid2][haloid2]['particle_IDs'].size


        shared_particleIDs = np.intersect1d(self.particleIDs[snapid1][haloid1]['particle_IDs'], self.particleIDs[snapid2][haloid2]['particle_IDs'], return_indices=False)
      
        shared_fraction1=100.0/n_particles1*shared_particleIDs.size
        shared_fraction2=100.0/n_particles2*shared_particleIDs.size

        if print_on_screen==True:
            self.print_on_screen(props1, z1, n_particles1, comparison=True, 
                                 props2=props2, z2=z2, n_particles2=n_particles2, 
                                 shared_fraction=shared_fraction1, shared_fraction2=shared_fraction2)           
            
        
        choose = {
            'haloid1':          haloid1,
            'haloid2':          haloid2,
            'n_particles1':     n_particles1,
            'n_particles2':     n_particles2,
            'snapid1':          snapid1,
            'snapid2':          snapid2,
            'z1':               z1,
            'z2':               z2,
            'n_particles_shared': shared_particleIDs.size,
            'n_particles_shared_perc1': shared_fraction1,
            'n_particles_shared_perc2': shared_fraction2,              
            'mhalo1':            props1['mhalo'],
            'mhalo2':            props2['mhalo'],
            'delta_mhalo':      props2['mhalo']-props1['mhalo'],           
            'delta_mhalo_perc': self.calc_delta_perc(props1['mhalo'],props2['mhalo']),             
            'rvir1':            props1['rvir'],
            'rvir2':            props2['rvir'],
            'delta_rvir':      props2['rvir']-props1['rvir'],
            'delta_rvir_perc':  self.calc_delta_perc(props1['rvir'],props2['rvir']),             
            'x_pos1':           props1['x_pos'],
            'y_pos1':           props1['y_pos'],
            'z_pos1':           props1['z_pos'],
            'x_pos2':           props2['x_pos'],
            'y_pos2':           props2['y_pos'],
            'z_pos2':           props2['z_pos'],              
            'delta_x_pos_perc': self.calc_delta_perc(props1['x_pos'],props2['x_pos']), 
            'delta_y_pos_perc': self.calc_delta_perc(props1['y_pos'],props2['y_pos']), 
            'delta_z_pos_perc': self.calc_delta_perc(props1['z_pos'],props2['z_pos'])             
            }

        bucket=np.zeros((1,), dtype=dt)
       
        for prop in choose.keys():
                bucket[prop] = choose.get(prop)
            
        return bucket

    def crossmatch_Cholla(self,
                          data):
        """Function crossmatches the properties of Cholla hydro with other tools e.g. ROCKSTAR
            UNDER CONSTRUCTION!
        """   
    
        return data
          

    def get_halo_info(self,
                       haloid,
                       snapid=None):
        
        print('\nFIND halo with ID:', haloid, 'at snapid/z:', snapid, '/', format(self.snapid_array['z'][np.where(self.snapid_array['snapid']==snapid)][0],'0.2f'),\
              '\n#######################################################################\n')
            
        props       = self.myDataHalo[np.where((self.myDataHalo['haloid']==haloid) & (self.myDataHalo['snapid']==snapid))[:][0]]
        z           = self.snapid_array['z'][np.where(self.snapid_array['snapid']==snapid)][0]
        n_particles = self.particleIDs[snapid][haloid]['particle_IDs'].size

        for i in range(0,props.size,1):
            self.print_on_screen(props[i], z, n_particles)


    def load_data_GADGET(self,
                          path, 
                          snapid_start, 
                          snapid_end=None, 
                          return_particle_info=True):
        """function reads and returns the GADGET halo catalog and particle information in binary and optional ASCII file at a given snapshot
            UNDER CONSTRUCTION
        """
        exit()

    def load_data_merger_trees(self,
                               treeID,
                               treeID_end,
                               file_basename):
        """Function reads and returns the merger trees catalog from ASCII file for a given filename and between certain tree ID numbers"""        
        
        merger_tree_data, myheader, mycols, myformat_string, mydt = ha_lib.create_data_array('merger_trees_ASCII')
        
        for i, tree_id in enumerate(range(treeID,treeID_end,1)):
            try:
                data = ha_lib.df_to_sarray(pd.read_csv(self.path_basename+file_basename+'_treeID'+str(tree_id)+'.txt', 
                                                       skiprows=2,
                                                       names=mycols,
                                                       sep='\t'))                    
            except:
                print('data for merger tree with treeID:', tree_id, 'not available! --> SKIPPED!')
            
            if i==0:
                merger_tree_data = data
            else:
                merger_tree_data = np.append(merger_tree_data, data, axis=0)
        
        return merger_tree_data
               
    
    def load_data_ROCKSTAR(self,
                          path, 
                          snapid_start, 
                          snapid_end=None, 
                          return_particle_info=True):
        """function reads and returns the ROCKSTAR halo catalog and particle information in binary and optional ASCII file at a given snapshot"""
    
        particleIDs_haloid_snapid_dict  = {}
        haloid_redshift_info            = {}
        
        if snapid_end==None:
            snapid_end=snapid_start+1

        for i, snapid in enumerate(range(snapid_start,snapid_end,1)):

            redshift = self.snapid_array[np.where(self.snapid_array['snapid']==snapid)]['z'][0]
            
            data_subfile, particleIDs, particle_snapid_info = ha_lc.load_binary_ROCKSTAR(path+'halo_files/', 
                                                                                 snapid=snapid, 
                                                                                 filename_prefix='halos_'+str(snapid)+'.', 
                                                                                 filename_suffix='.bin', 
                                                                                 n_files=self.n_files_snapshot,
                                                                                 myprops=self.myprops2read)
    
            haloid_redshift_info.update({snapid: {'redshift': redshift, 'haloid': data_subfile['haloid']}, 'halo_info': particle_snapid_info})
            particleIDs_haloid_snapid_dict.update({snapid: particleIDs})
                    
            if i==0:
                data = np.zeros((data_subfile.size,),dtype=np.dtype([(k, data_subfile.dtype[k]) for k in data_subfile.dtype.names]))
                data = data_subfile
            else:
                data = np.append(data, data_subfile, axis=0)

        if return_particle_info==True:
                        
            return data, particleIDs_haloid_snapid_dict, haloid_redshift_info
        
        else:
            return data
     
    def construct_merger_tree(self,
                              haloid=None,
                              n_trees=[],
                              snapid=None,
                              snapid_end=None,
                              print_merger_tree2file=True):
        """Function constructs merger trees of detected halos using ROCKSTAR halo binary information and descendent IDs 'DescID' from out_x.list
        catalogs (where x stands for the snapshot number)
        
        The main prognitor haloids are printed into an ASCCII file along with the redshift of the first and last detection of that halo and a unique ID
        of the tree
        
        output:
        ===============
            returned by the function
            data    structured array of all halo properties of all halo found in the box
            merger_dict     dictionary of the merger tree of each unique tree ID (rootIndex)
            z_first_app     dictionary of the first appearance of 
            
            printed into file each row give one tree in the format ASCII --> z(first appearance) z(last appearance) rootIndex [all haloids of the main progenitors]
        """
        merger_dict={}
        z_first_app={}
        z_last_app={}        

                        
        if haloid==None and n_trees==[]:
             tree_list2analyse=range(0,1000,1)           
        elif haloid!=None:           
            if type(haloid)==list:
                tree_list2analyse=haloid
            else:
                tree_list2analyse=[haloid]
        else:                
            tree_list2analyse=range(0,n_trees,1)               
            
        if snapid!=None:
            mysnap_array=self.snapid_array[np.where(self.snapid_array['snapid']==snapid_end)[:][0][0]:np.where(self.snapid_array['snapid']==snapid)[:][0][0]+1]
        else:
            mysnap_array=self.snapid_array
        
        data_before = self.load_data_ROCKSTAR(self.path_basename,
                                              np.sort(mysnap_array['snapid'])[0],
                                              return_particle_info=False)
        data_before, data_dummy = self.assign_progenitor_indices(data_before, data_before, snapid, first_assignment=0)
       
        for tree_id in tree_list2analyse:
            progs_tree=[]
            print('tree_id:', tree_id, '/', len(tree_list2analyse)-1, '... processing!')
                
            count_first_app=0                 
            count_last_app=0            
            for i, snapid in enumerate(np.sort(mysnap_array['snapid'])):
                
                if i==0:
                    data=data_before
                else:
                    data = np.append(data, data_before, axis=0)                
                
                new_data, progID, data_before = self.get_progenitor_data(self.path_basename,
                                                                         data_before,
                                                                         snapid,
                                                                         tree_id)
 
                if len(progID)>1:
                    progID=[progID[0]]

                progs_tree.extend(progID)
                
                if i==0 and tree_id==tree_list2analyse[0]:
                    data_tree = new_data
                else:
                    data_tree = np.append(data_tree, new_data, axis=0)
                                       
                if len(new_data)==0:
                    count_last_app+=1
                else:
                    count_first_app+=1
                    snap_last_app=snapid
                  
            if progs_tree!=[]:                          
                merger_dict.update({tree_id:progs_tree})
                z_first_app.update({tree_id:snap_last_app-count_first_app+1})
                z_last_app.update({tree_id:snap_last_app})                
   
        if print_merger_tree2file==True: 
            filename_trees=self.path_basename+'merger_trees_MM.txt'     

            for tree_id in merger_dict.keys():
                if tree_id==tree_list2analyse[0]:
                    ha_lib.writeIntoFile(
                               filename_trees,
                               [str(z_first_app[tree_id])+' '+str(z_last_app[tree_id])+' '+str(tree_id)+' '+" ".join(str(item) for item in merger_dict[tree_id][1:])],
                               myheader=str(self.path_basename)+'Main progenitor merger trees as haloid (ID in the orginal ROCKSTAR files)\n'\
                                            +'(1) snapid of first appearance of this tree with ID (2) unique ID of the merger tree (descendent) (3)-(n) progenitors',
                               append_mytext=False,
                               data_is_string=False,
                               data_format='%s')
                else:
                    ha_lib.writeIntoFile(
                               filename_trees,
                               str(z_first_app[tree_id])+' '+str(z_last_app[tree_id])+' '+str(tree_id)+' '+" ".join(str(item) for item in merger_dict[tree_id][1:])+'\n',
                               append_mytext=True,
                               data_is_string=True,
                               data_format='%s')
        
        return data, merger_dict, z_first_app
    
    
    def get_progenitor_data(self,
                            path,
                            data_before,
                            snapid,
                            tree_id,
                            first_assignment=1):
        """function access the progenitor information finds the correct (most massive) progenitor index to a given haloid"""

        data = self.load_data_ROCKSTAR(self.path_basename,
                                       snapid,
                                       return_particle_info=False)
            
        data, data_before = self.assign_progenitor_indices(data, data_before, snapid, first_assignment=first_assignment)
        
        tree_data=data[np.where(data['rootIndex']==tree_id)[:][0]]
        tree_data[::-1].sort(order='mhalo', axis=0)
    
        import numpy.lib.recfunctions as rcfuncs
        tree_data = rcfuncs.append_fields([tree_data], ['Z'] , [np.zeros(tree_data.size,)], dtypes=['f4'], usemask=False)
        
        tree_data['Z']=self.snapid_array[np.where(self.snapid_array['snapid']==snapid)]['z'][0]         
                         
        return tree_data, tree_data['haloid'], data
 
    
    def print_on_screen(self,
                        props, z, n_particles, comparison=False, 
                        props2=None, z2=None, n_particles2=None, shared_fraction=None, shared_fraction2=None):
        """Function prints halo properties on the screen either of one halo or comparses two haloos"""
        
        if comparison==True:
            print('properties \thaloid1     |\thaloid2     |\tDelta(haloid2-haloid1) / [%]\n--------------------------------------------------------------------------\n',\
                  'z:\t\t', str(format(z, '0.2f')).ljust(10), ' |\t', str(format(z2, '0.2f')).ljust(10), ' |\t')
            
            for item in props.dtype.names:
                if item.find('Index')==1 or item.find('id')==1 or item.find('n_')==1:
                    print(item, ':\t', str(props[item][0]).ljust(10),' |\t', str(props2[item][0]).ljust(10), ' |\t')                        

                print('n_particles:\t', str(n_particles).ljust(10), ' |\t', str(n_particles2).ljust(10), ' |\t',)                
                
            try:
                print(n_particles2-n_particles, ',', format(shared_fraction, '0.1f'), '/', format(shared_fraction2, '0.1f'), '[%] shared IDs')
            except:
                print('unknown [%] shared IDs')  
            
        else:
            print('properties\n--------------------------------------------------------------------------\n',\
                  'z:\t\t', str(format(z, '0.2f')).ljust(10))                
            for item in props.dtype.names:
                if item.find('Index')==1 or item.find('id')==1 or item.find('n_')==1:
                    print(item, ':\t', str(props[item][0]).ljust(10))                        
            
            print('n_particles:\t', str(n_particles).ljust(10))                

        for item in props.dtype.names:

            if item.find('Index')==-1 and item.find('id')==-1 and item.find('n_')==-1:
                if comparison==True:
                    print(item, ':\t\t', format(props[item], '0.4e'),' |\t', format(props2[item], '0.4e'), ' |\t', format(props2[item]-props[item], '0.4e'), '/', format(100.0-100.0/props2[item]*props[item], '0.2f'), '[%]')
                else:
                    print(item, ':\t\t', format(props[item], '0.4e'))
        
        print('--------------------------------------------------------------------------\n')            

    def track_merger_tree(self, snapid_first, snapid_last, rootID, progIDs):

        data, myheader, mycols, myformat_string, mydt = ha_lib.create_data_array('merger_trees_ASCII')        

        for i, snapid in enumerate(range(snapid_first,snapid_last,1)):
            
            halo1 = self.myDataHalo[np.where((self.myDataHalo['descIndex']==progIDs[i]) & (self.myDataHalo['snapid']==snapid))[:][0]]
            halo2 = self.myDataHalo[np.where((self.myDataHalo['haloid']==progIDs[i]) & (self.myDataHalo['snapid']==snapid+1))[:][0]]
            
            halo1[::-1].sort(order=['mhalo'], axis=0)
            halo2[::-1].sort(order=['mhalo'], axis=0) 
           
            print('halo1:', halo1[['haloid','descIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos', 'snapid']],\
                      '\nhalo2:', halo2[['haloid','descIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos', 'snapid']])
            try:
                mydata=self.compare_halo_info(props1=halo1[0], haloid1=halo1['haloid'][0], props2=halo2[0], haloid2=halo2['haloid'][0], snapid1=snapid,snapid2=snapid+1, dt=mydt)
                mydata['rootIndex']=rootID
               
                mydata['descIndex1']=halo1['descIndex'][0]
                if i==0:
                    data = mydata
    
                else:
                    data = np.append(data, mydata, axis=0)
            except:
                return False

        if len(progIDs)>0:         
            ha_lib.writeIntoFile(
                       self.path_basename+'data_comparision_treeID'+str(rootID)+'.txt',
                       data[mycols],
                       myheader='All halos fro Cholla 50Mpc 256 box from Bruno run 29/01/2021\n'+myheader,
                       mydelimiter='\t',
                       data_format=myformat_string)      
            
        return data
       
       
    def verify_merger_trees(self):  
        """Function compares the particle content in the Cholla particle.h5 files with the halo mass found by ROCKSTARS
            UNDER CONSTRUCTION!
       """       
        cube_size=25000.0
        dx=195.3125
        n_rvir=1.5
        perc_delta_halo=10
        
        filename_check_trees='anaconda/pro/data/Cholla256_50Mpc/Cholla256_50Mpc_check_merger_trees_'+str(n_rvir)+'Rvir.txt'
        
        header_prefix='ROCKSTAR run on Cholla256 50h-1Mpc, 0.81 < z < 11.51, main progenitors crossmatch with particle field in Cholla particle file. '\
                        'Particles around '+str(n_rvir)+' x rvir [h-1comvkpc] where used.\n'
                        
        header_endfix=' (15) n_particle [count] (Cholla particle data) (16) halo mass [h-1Msun] (n_particle x particle mass (Cholla particle data)) '\
                      '(17) delta(mhalo_particle-mhalo_Rockstar) [h-1Msun] (18) percentage delta(mhalo_particle-mhalo_Rockstar) [%] '\
                      '(19) check delta(mhalo_particle-mhalo_Rockstar) [flag: (18)<'+str(perc_delta_halo)+'%==1,else==0]'                    
        
        for snapid,z in zip(self.snapid_array['snapid'][10:], self.snapid_array['z'][10:]):        
            print('snapid:', snapid, 'redshift:', z,)
                    
            data_prog = self.HaloData[np.where(self.HaloData['Z']==z)[:][0]]
            
            print('n_progs:', data_prog.size, '\n////////////////////////////////////////////////////////\n')
            
            for k, prog in enumerate(data_prog[['x_pos', 'y_pos', 'z_pos']]):
                
                print('prog count:', k+1,'/', data_prog.size, 'haloid:', data_prog['haloid'][k], 'treeID:', data_prog['firstProgenitorID'][k],  'position(X,Y,Z) [h-1comvkpc]:', prog, '\n===============================================\n')
        
                #which cell is it?
                cell_loc = ha_lib.find_cell_location(prog, dx)
             
                #In which sub-file is the cell stored!'
                subcube_loc = ha_lib.find_cell_location(prog, cube_size)
        
                print('in sub-cube:', subcube_loc, '--> decimal:',)
                binary_rep=int(str(subcube_loc[0])+str(subcube_loc[1])+str(subcube_loc[2]), base=2)
                print(binary_rep)
            
                if binary_rep==4:
                    binary_rep=1
                elif binary_rep==1:
                    binary_rep=4
                elif binary_rep==3:
                    binary_rep=6
                elif binary_rep==6:
                    binary_rep=3     
         
                #load particle data -- positions

                path_p=str(snapid)+'_particles.h5.'+str(ha_lib.binaryToDecimal(binary_rep))

                data_particle_pos= ha_lib.simple_hdf52struct_array(path_p, ['pos_x', 'pos_y', 'pos_z'])
                
                #load particle data -- density grid
                #data_particle_grid= simple_hdf52struct_array(path_p, ['density', 'grav_potential'])
                
                #load hydro data
                #path_f = '/data/256_hydro_50Mpc/'+snapid+'.h5.'+str(binaryToDecimal(binary_rep))
                #data_hydro= simple_hdf52struct_array(path_f)
                
                #print np.info(data_particle_pos)    
                #print np.info(data_particle_grid)    
                #print np.info(data_hydro)
                
                #scan_file_format_Cholla(f_p, path_p)
                   
                print('r_vir:', format(data_prog['rvir'][k], '0.2f'), '[h-1kpc]',)        
                
                n_particle = len(ha_lib.find_members(data_particle_pos[['pos_x','pos_y','pos_z']], prog, n_rvir*data_prog['rvir'][k])[:][0])
                
                print('n_particles:', n_particle, 'mhalo_particle:', format(n_particle*m_particle , '0.5e'), '[h-1Msun]')
                
                output_string=''
                myheader=''
                
                myprops={'0': {'name': 'haloid', 'unit':'-'}, '1': {'name':'descIndex', 'unit':'-'},'2': {'name':'predIndex','unit':'-'}, '3': {'name':'firstProgenitorID','unit':'-'},\
                         '4': {'name':'npros','unit':'count'}, '5': {'name':'mhalo','unit':'h-1Msun'}, '6': {'name':'delta_mhalo','unit':'h-1Msun'}, '7': {'name':'rvir','unit':'h-1comkpc'},\
                         '8': {'name':'delta_rvir','unit':'h-1comvkpc'}, '9': {'name':'x_pos','unit':'h-1comvkpc'}, '10': {'name':'y_pos','unit':'h-1comvkpc'}, '11': {'name':'z_pos','unit':'h-1comvkpc'},\
                         '12': {'name':'Z','unit':'-'}}
        
                for i in range(len(myprops)):
                    #print 'i:', i, myprops[str(i)]['name']
                    myheader+='('+str(i+1)+') '+myprops[str(i)]['name']+' ['+myprops[str(i)]['unit']+'] '
                    if myprops[str(i)]['name'].find('mhal')!=-1:
                        output_string+=str(format(data_prog[k][myprops[str(i)]['name']], '0.5e'))+' '
                    else:
                        output_string+=str(data_prog[k][myprops[str(i)]['name']])+' '
                
                #print myheader
                #print output_string
                
                delta_mhalo = n_particle*self.m_particle-data_prog[k]['mhalo']
                
                delta_mhalo_per = 100.0/data_prog[k]['mhalo']*delta_mhalo
                
                if delta_mhalo_per>-perc_delta_halo and delta_mhalo_per<perc_delta_halo:
                    check_delta_mhalo=1
                else:
                    check_delta_mhalo=0                
        
                if str(snapid)=='21' and k==0:
                   ha_lib.writeIntoFile(
                               filename_check_trees,
                               [output_string+' '+str(n_particle)+' '+str(format(n_particle*m_particle, '0.5e'))+' '+str(format(delta_mhalo, '0.5e'))+' '+str(format(delta_mhalo_per, '0.1f'))+' '+str(check_delta_mhalo)],
                               myheader=header_prefix+myheader+header_endfix,
                               append_mytext=False,
                               data_is_string=False,
                               data_format='%s')
                else:
                    ha_lib.writeIntoFile(
                               filename_check_trees,
                               output_string+' '+str(n_particle)+' '+str(format(n_particle*m_particle, '0.5e'))+' '+str(format(delta_mhalo, '0.5e'))+' '+str(format(delta_mhalo_per, '0.1f'))+' '+str(check_delta_mhalo)+'\n',
                               append_mytext=True,
                               data_is_string=True,
                               data_format='%s')
                    
                print('\n===============================================\n')
                                
#initalize class and call MAIN!
HaloAnalysis().MAIN(task)      