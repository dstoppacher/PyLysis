"""Halo, merger tree handling and analysis Last update 2021-02-26 by DS

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

class HaloAnalysis:
    
    def __init__ (self,
                  pathname_dict,
                  config_dict):
        self.pathname_dict = pathname_dict
        self.config_dict = config_dict
        try:
            self.snapid_array = ha_lib.df_to_sarray(pd.read_csv(self.pathname_dict['path_basename']+'snapidzred.txt', 
                                               skiprows=2, 
                                               names=['snapid', 'a', 'z', 't'], 
                                               sep='\t'))

        except:
            print('Could not read snapidzred.txt!')
        
        print('\nInitiating HaloAnalysis pipeline ....\n++++++++++++++++++++++++++++++\n')
  
    def MAIN(self, task):   
               
        def create_DescScales_file():
            """helper function: access simulation data and on basis of that creates and stores a list desc scales as needed by Consitent-Trees to generate merger trees"""            
            ha_lib.create_DescScales_file(self.pathname_dict['path_basename']+'halo_files/',
                                          self.config_dict['snapid'],
                                          self.config_dict['snapid_end'])
        def create_snapidzred_file():
            """helper function: access simulation data and on basis of that creates and stores a list of snapshots, redshift, scale factors present in the simulation"""
            ha_lib.create_snapidzred_file(self.pathname_dict['path_basename'],
                                          self.pathname_dict['path_orginal_data_hydro'],
                                          self.config_dict['snapid'],
                                          self.config_dict['snapid_end'])            
        
        def load_entire_box():
            """load all snapshots and particle information between the snapid and the snapid_end"""
            self.myDataHalo, self.particleIDs, self.halo_redshift_info = self.load_data_ROCKSTAR(self.pathname_dict['path_basename'],
                                                                                                 self.config_dict['snapid'],
                                                                                                 self.config_dict['snapid_end']+1)

        def get_halo_info():
            """get infos of a certain halos at a certain snapshot
                input:
                    haloid  id number of the halo you want the info
                    snapid  snapshot number where the halo was found with haloid
                    
                output:
                    info of halo properties from self.config_dict[myprops2read] printed on the screen
            """
            load_entire_box()
            self.get_halo_info(haloid=self.config_dict['haloid'], snapid=self.config_dict['snapid'])

        def get_tree_statistics():
            
            #fname=self.pathname_dict['path_basename']+'halo_files/trees/Cholla512_50Mpc_CT-v1.01_bin'
            fname=self.pathname_dict['path_basename']+'halo_files/trees/VSMDPL3840_160Mpc_tree_2_5_0_CT-v1.01_MB_'            
            
            count_trees=0
            import numpy.lib.recfunctions as rcfuncs
            
            for i in range(0,4,1):
                print('loading default:', fname+str(i)+'.txt')
                mycolnames=['mhalo']
                data = load_file(fname+str(i)+'.txt',get_col_names('ConsistentTrees_basic_ASCII'))

                for prop in mycolnames:                
                    data = rcfuncs.append_fields([data], ['delta_'+prop,'delta_'+prop+'_perc'] , [np.zeros(data.size,),np.zeros(data.size,)], dtypes=[ha_lib.get_property_dict()[prop]['dtype'],np.float32], usemask=False)
    
                data_dict=self.create_dict_unique_IDs(data,'rootIndex')
                
                for count,ids in enumerate(data_dict.keys()):                            
                    index_rootID = np.where(data['rootIndex']==ids)
                    #print('count:', count, 'ids:', ids)

                    data_index=data[index_rootID]
                    for j in np.unique(data_index['subTreeID']):
                        index_subTreeID=np.where(data_index['subTreeID']==j)                                
                        #print('\tsubTreeID:', j)
                        
                        data_index[index_subTreeID]=self.calc_MAH_from_tree(data_index[index_subTreeID],[prop])
                        count_trees+=j
                        
                    data[index_rootID]=data_index                              

                    count_trees+=count
                print(count,'rootIndices processed! Total amount of trees are ', count_trees)
                #print(data[['snapid','mhalo','delta_mhalo','delta_mhalo_perc','haloid_CT','rootIndex','subTreeID']][0:200])      
                    
                self.get_merger_tree_stats(data,
                                           filename_endfix='_bin'+str(i),
                                           colnames=['mhalo','delta_mhalo','delta_mhalo_perc'])
            

        def compare_individual_halos(): 
            """compare properties and print percentage of particle they share (between snapshots or at one snapshot)
                input:
                ------------
                    loads entire box!
                    haloid1     id of the first halo you want to compare
                    haloid2     id of the second halo you want to compare
                    snapid1     snapshot number of the first halo you want to compare
                    snapid2     snapshot number of the second halo you want to compare
                    
                keywords:
                 ------------                   
                    print_on_screen  default: False
                                        info of halo properties from self.config_dict[myprops2read] printed on the screen
                                        difference between the halo properties and particle content are calculate
                                        as absolute vales or as percentage where the properties of haloid2 (lower redshif) are compared
                                        to those of haloid1 (higher redshift)                   
                    
                output:
                ------------                    
                    structured array of one column which contains the information                  

                    
            """
            load_entire_box()
            self.compare_halo_info(haloid1=self.config_dict['haloid'], 
                                   haloid2=self.config_dict['haloid2'], 
                                   snapid1=self.config_dict['snapid'], 
                                   snapid2=self.config_dict['snapid']+1, 
                                   print_on_screen=True)        

        def construct_main_prog_tree_RS():
            """Function constructs main progenitor tree form ROCKSTAR data using the 'DescID'
                This approach does not use a tree builder to verify the connection between halos and fixes them!
                BEWARE USING THIS FUNCTION!
            """
            start = time() 
    
            data, merger_dict, first_z_dict = self.construct_main_prog_tree(haloid=45,
                                                                            snapid=self.config_dict['snapid'],
                                                                            snapid_end=self.config_dict['snapid_end'])
            

            print('Time until all halos tracked:', format((time()-start)/60.0/60.0, '0.2f'), 'h /', format((time()-start)/60.0, '0.2f'), 'min /', format((time()-start), '0.1f'), 'sec')

            ha_lib.writeIntoFile(
                       self.pathname_dict['path_basename']+self.pathname_dict['output_filename_data'],
                       data[['haloid','descIndex', 'rootIndex', 'predIndex', 'mhalo', 'delta_mhalo', 'rvir', 'delta_rvir', 'x_pos', 'y_pos', 'z_pos', 'snapid', 'n_particles', 'npros']],
                       myheader='Cholla'+str(self.config_dict['res'])+' '+str(self.config_dict['box_size_h-1Mpc'])+'h-1Mpc ROCKSTAR Main Branch all halos\n'\
                       +'(1) haloid (2) descIndex (3) rootIndex (4) predIndex (5) Mvir [h-1Msun] (6) delta_Mvir [h-1Msun] (7) Rvir [h-1kpc] (8) delt_Rvir [h-1kpc]'\
                       +'(9) X [h-1Mpc] (10) Y [h-1Mpc] (11) Z [h-1Mpc] (12) snapid (13) n_particles (14) n_progs',
                       mydelimiter='\t',
                       data_format='%i\t%i\t%i\t%i\t%0.6e\t%0.6e\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%i\t%i\t%i')     



        def identify_merger_trees():
            """Function hubs to self.get_progenitor_tree_branches_CT"""
            
            tree_data=load_CT()
            #print(tree_data)
            #print(np.info(tree_data))
            #test_ids=tree_data['rootIndex'][0:10]
            #print(test_ids)
            self.get_progenitor_tree_branches_CT(tree_data,
                                                   rootIndex=[],
                                                   print_all_trees=True,
                                                   print_rootIDs_sep=False,
                                                   print_main_branch=True,
                                                   sample_endfix='')
            


        def walk_merger_trees(generate_merger_tree_file=True):
            """Function identifies and walks the merger trees of all progenitor branches in the merger tree calalog of choice.
            
            Thereby all trees linked to a certain unique topnode id are found and subsequently processed. The merger tree halo catalog can be usually
            found here: self.pathname_dict['path_basename']+self.pathname_dict['output_filename_merger_trees'] and be loaded to the structured array 'mytree_data'.
            Alternatively, the merger tree catalog can be generated on the fly from the halo catalog if 'generate_merger_tree_file' is set to 'True'. Then the
            function 'self.generate_merger_tree_dict' is called.
            
            Each merger tree is stored in a dictionary 'mytree_dict' where the keys are the ascending number count of the total number of trees found
            in the simulation box
            
            input:
            ------------                
                loads entire box
                
            keywords:
            ------------                
                generate_merger_tree_file  default False, meaning that if keyword is set to TRUE, the merger tree which should be tracked is greated on the fly otherwise
                                           the a file is loaded which contrain that information --> self.pathname_dict['path_basename']+self.pathname_dict['output_filename_merger_trees']
            output:
            ------------ 
                prints trackes merger trees into a file if  'generate_merger_tree_file' is set TRUE, checks how many merger trees have been correctly tracked.
            """
            
            load_entire_box()
            tree_indices=self.get_progenitor_tree_branches_CT(print_into_file=False)
           
            unique_rootIndices = np.unique(tree_indices['rootIndex'])
            
            if generate_merger_tree_file==True:
                #generate tree data from scratch
                mytree_dict = self.generate_merger_tree_dict(unique_rootIndices, tree_indices)               
            
            else:
                #Load tree data
                mytree_data = ha_lib.read_unshaped_txt(self.pathname_dict['path_basename']+self.pathname_dict['output_filename_merger_trees'], self.config_dict['treeID_end']+10, self.config_dict['total_snaps_sim'], delimiter=' ', comments='#', dtype=np.int32)
                mytree_dict={}                                      
                for i in range(0,self.config_dict['treeID_end']+1,1):
                    if self.total_snaps2read>=mytree_data[i,1]-mytree_data[i,0]:
                        end_snap=mytree_data[i,1]-mytree_data[i,0]
                    else:
                        end_snap=self.config_dict['total_snaps2read']
                        
                    print('i:', i, 'snapid_first:', mytree_data[i,0], 'snapid_last:', mytree_data[i,1], 'n_snaps2read:', end_snap, 'rootIndex:', mytree_data[i,2], 'progIDs:', mytree_data[i, 3:end_snap+3], 'check_sum:', len(mytree_data[i, 3:end_snap+3]))
                       
                    mytree_dict.update({i: (mytree_data[i,0], mytree_data[i,0]+end_snap, mytree_data[i, 3:end_snap+3])})                   
                
            #print(mytree_dict)
                                              
            count_error=0
            failed_IDs=[]

            
            for index in unique_rootIndices:            

                print('rootIndex:', index, 'snapid_first:', mytree_dict[index][0], 'snapid_last:', mytree_dict[index][1], 'n_snaps2read:', len(mytree_dict[index][2]), 'progIDs:',mytree_dict[index][2])
                
                first_sn_found  = mytree_dict[index][0]
                last_sn_found   = mytree_dict[index][1] 
                progIndices     = mytree_dict[index][2]

                tree_structure=self.track_merger_tree(first_sn_found, last_sn_found, index, progIndices)

                if tree_structure==False:
                    count_error+=1
                    failed_IDs.extend([index])
                print('--> completed')

            print('\n-------------------\nSUMMARY:\n', len(unique_rootIndices)-count_error+1, 'trees sucessfully walked / ', count_error, 'failed!', '\nfailed treeIDs are:\n', failed_IDs)      

        def get_col_names(name,cols=[]):
            
            dt = ha_lc.get_dtype_sarray(name,colnames=cols)            
            #print(dt)
    
            return [k[0] for k in dt]
        
            
        def load_file(filename,col_names):
            
            return ha_lib.df_to_sarray(pd.read_csv(filename,
                                       skiprows=2, 
                                       names=col_names,
                                       sep='\t'))            
        
        def load_CT():            
            
            return ha_lc.load_ASCII_ConsitentTrees(self.pathname_dict['CT_merger_trees'])

        def load_CT_MB():            
            #load all main branches of all rootIndices
            return ha_lc.load_ASCII_ConsitentTrees(self.pathname_dict['CT_merger_trees'])            
        
            
        def load_merger_trees(filename):
                            
            return load_file(filename,get_col_names('merger_trees_ASCII'))
           
           
        def load_stats(filename):
            
            return load_file(filename,get_col_names('stats_perc_bucket'))            



        def load_ascii_file():
            """load a random text file in ascii code into a pandas data frame and then transform it to a structured array. Get the data types
                from a dictionary (dt) if ascii file was also reated by PyLysis
                    
            output:
            ------------                   
                structured array with data types from dt which can be easily transfered to jupyter-lab!
            """
          
            key='rootIDr'
            filenames = ['data_treeID79509','data_comparision_MM_treeID3', 'data_comparision_MM_treeID6', 'data_comparision_treeID54', 'data_comparision_treeID231']   #key: MT --> main branch data standard merger trees       
            filenames = ['308328'] #rootIDs   use key:rootID
            filenames = [self.pathname_dict['output_filename_merger_trees']] #load all main branches   use key:MB
            #filenames = ['509', '561', '604', '619','622']               
            #filenames = ['Cholla256_50Mpc_CT-MB_stats'] #key: stats --> CT main branch stats
            #filenames = ['halo_files/trees/tree_0_0_0.dat'] #key: CT --> CT original files
            #filenames = [self.pathname_dict['path_basename']+'halo_files/trees/Cholla512_50Mpc_CT-v1.01_bin'] #key: CT --> CT original files or and other file in this format
            
            filenames = ['0','1','2','3']
            data_dict={}
            for i, fname in enumerate(filenames):
                #print('i:', i, 'fname:', fname)
 
                if key=='CT':
                    print('loading Consitent tree file:', fname)
                    data = load_CT()
                    
                elif key=='MT':
                    print('loading standard merger tree file:', fname)
                    data=load_merger_trees(fname)
                    
                elif key=='stats':
                    print('loading stats file:', fname)
                    data=load_stats(fname)

                elif key=='MB':
                    print('loading all main progenitor trees from file:', fname)
                    
                    data = load_file(fname,get_col_names('ConsistentTrees_basic_ASCII'))
                    #Sort by a certain properties e.g. if the first 10 most massive rooIndex-halos should be found
                    data[::-1].sort(order=['mhalo'], axis=0)
                    #print(data[['rootIndex', 'mhalo']][0:100])
                    data_dict=self.create_dict_unique_IDs(data,'rootIndex')                  
                   
                elif key=='rootID':
                    fname=self.pathname_dict['output_filename_merger_trees'][:-4]+'_rootID'+fname+'.txt'
                    print('loading ...', fname)
                    
                    data_dict=self.create_dict_unique_IDs(load_file(fname,get_col_names('ConsistentTrees_basic_ASCII')),'subTreeID')
                else:

                    filename=self.pathname_dict['path_basename']+'halo_files/trees/Cholla512_50Mpc_CT-v1.01_bin'+fname+'_stats.txt'
                    
                    
                    
                    mydata=load_file(filename,get_col_names('stats_custom_bucket',cols=['mhalo','delta_mhalo','delta_mhalo_perc']))
                    

                    data_dict.update({fname: mydata})
                    
                        
            #print('data_dict:', data_dict)

            return data_dict
                  
            
        def utilise():
            """playground for anykind of manipulating and arranging data"""

            #data=ha_lc.load_ASCII_ConsitentTrees(self.pathname_dict['CT_merger_trees'])
            #load_entire_box()
            #calculate histogramm
            #snap_list = [154, 150, 140, 130, 100, 50, 10] #256_hydro
            #snap_list = [169,130,106,90,74,62] #512_hydro
    
            #self.calc_histo(self.myDataHalo, 'mhalo_200c', 'h-1Msun', snap_list, 'HMF')
            #exit()
            ######################################################################################################


            #Get tree statistics
                 
            tree_data=self.filter_trees(snapid=169,
                                        create_bins=True,
                                        prop='mhalo',
                                        bins=[0,1e11,1e12,1e13])
            #bins 512-Cholla
#            bin0: >=0 & <1e+11 	n_halos: 1534
#            bin1: >=1e+11 & <1e+12 	n_halos: 2697
#            bin2: >=1e+12 & <1e+13 	n_halos: 199
#            bin3: >1e+13 		n_halos: 8               

            #Load treedata f
#            tree_data=self.load_data_merger_trees_from_text(self.config_dict['rootID'],self.config_dict['rootID_end'],file_basename='/Cholla'+str(self.config_dict['res'])+'_'+str(self.config_dict['treeID'])+'Mpc_CT-MB/data_CT')
#            
            for item in tree_data.keys():
                print('item:', item, tree_data[item])           

        def caseSwitcher(task):
            choose = {
                    'compare_individual_halos':     compare_individual_halos,
                    'construct_main_prog_tree_RS':  construct_main_prog_tree_RS,
                    'create_DescScales_file':       create_DescScales_file,
                    'create_snapidzred_file':       create_snapidzred_file,
                    'identify_merger_trees':        identify_merger_trees,                      
                    'get_halo_info':                get_halo_info,
                    'get_tree_statistics':          get_tree_statistics,                    
                    'load_ascii_file':              load_ascii_file,
                    'load_entire_box':              load_entire_box,       
                    'walk_merger_trees':            walk_merger_trees,
                    'utilise':                      utilise
                    }
        
            func = choose.get(task)
            
            return func()
    
        return caseSwitcher(task)
    
        
        ##############################################################################################################
      
    def assign_progenitor_indices(self,
                                  data_RS_now,
                                  data_RS_before,
                                  snapid_before,
                                  first_assignment=1,
                                  show_info=False):
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
                                            dtypes=['f4','f4','i8','i8','i4'], 
                                            usemask=False)
       
        data_RS_now['npros']    = -99
        data_RS_now['rootIndex']= -99
        data_RS_now['predIndex']= -99

        if first_assignment==0:        
            data_RS_before = rcfuncs.append_fields([data_RS_before], 
                                                   ['npros','rootIndex','predIndex'] , 
                                                   [np.zeros(data_RS_before.size,),np.zeros(data_RS_before.size,),np.zeros(data_RS_before.size,)], 
                                                   dtypes=['i4','i8','i8'], 
                                                   usemask=False)
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

        if show_info==True:        
            print('data_RS_before:\n', data_RS_before[['haloid','descIndex', 'predIndex', 'rootIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos']],\
                  '\ndata_RS_now:\n', data_RS_now[['haloid','descIndex', 'predIndex', 'rootIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos', ]],\
                  'PROGENITOR INFO successfully assigned!\n#################################################')           

        return data_RS_now, data_RS_before
    
    
    def calc_delta_perc(self, prop1, prop2):
        """calculates the decrease or increase of prop2 in comparison with prop1 in percent"""
        
        return 100.0/prop1*prop2-100.0

    def calc_histo(self,
                   data,
                   prop,
                   unit,
                   snap_list,
                   histo_name):
        """caculates histogram with 'calcFastHisto' function
        
        input:
            prop    galaxy/halo property the histogram should be calculated (e.g. mhalo_200c, mstar, sfr)
            unit    the current unit of the property the histogram should be calculated (e.g. 'h-1Msun', 'yr-1Msun')
            snap_list   list of snapshot number the histogram should be calculated, they will be converted to redshifts
            histo_name  preffix of the histogram filename e.g. 'HMF' for 'halo mass function', 'SFRF' for star formation rate function etc.
            
        output:
            histogram is stored into 'self.pathname_dict['path_basename']+'histos/[histo_name]_CHOLLA_[self.config_dict['res']_self.config_dict['box_size_h-1Mpc']Mpc_z_[redshift]_[prop].txt'
        """

        for snap in snap_list:
            data_histo=data[np.where(data['snapid']==snap-self.config_dict['snapid_CT_offset'])[:][0]]
            redshift=format(ha_lib.expfactor_to_redshift(data_histo['a'][0]), '0.2f')
            print('snap:', snap, 'CT-snap:', snap-self.config_dict['snapid_CT_offset'], 'a:', data_histo['a'][0], 'z:', redshift, 'nhalos:', data_histo.size)
            ha_lib.calcFastHisto(data_histo[prop],\
                                 self.pathname_dict['path_basename']+'histos/'+histo_name+'_Cholla'+str(self.config_dict['res'])+'_'+str(self.config_dict['box_size_h-1Mpc'])+'Mpc_z_'+str(redshift)+'_'+prop+'.txt',\
                                 prop, unit, \
                                 self.config_dict['box_size_h-1Mpc'], \
                                 comment='Cholla-'+str(self.config_dict['res'])+' '+str(self.config_dict['box_size_h-1Mpc'])+'h-1Mpc all halos z='+str(redshift)+' nhalos: '+str(data_histo.size))

    def calc_MAH_from_tree(self,
                         data,
                         colnames):
        
        """function calculation the evolution of properties with redshift and compares the property at redshift snapid+1 with the value at sanpid. Where the 
        snapids decrease with increasing redshift
        
        input:
            data        a structured array
            colnames    the property where the redshift evolution or mass accrection history shuld be calculated (e.g. mhalo, rvir)
            
        output:
            data        same structrued array but with an additional column 'delta_'+[name of property]
            
        """
        
        #print('HERE 547\n///////////////////////////////////////////////////////////////////////////////////////////////////////////')
        #print(np.info(data))
        data[::-1].sort(order=['snapid'], axis=0)

        for prop in colnames:
            #print('prop:',data[prop])
            data['delta_'+prop][0:data.size-1] = data[prop][0:data.size-1]-data[prop][1:data.size]
            data['delta_'+prop+'_perc'][0:data.size-1]=self.calc_delta_perc(data[prop][1:data.size],data[prop][0:data.size-1])
            
        #print(data[['snapid',prop,'delta_'+prop,'delta_'+prop+'_perc','haloid_CT','subTreeID']])

        return data


    def give_stats_info(self,
                         data,
                         props):
        """The method return info such as max/min-values, median/mean and error estimation of data from a structred array
        
        input:
        ------------            
            data    a structured array
            props   a list of properties which should be displayes
            
        output:
        ------------            
            the statistics is printed on screen
        """
        
        print('##############################################\n')              
        for prop in props:
            unit, myformat = ha_lib.find_unit(prop)           
            print(('prop').ljust(25), 'median 10th / 90th\tmax/min\t[unit]\n-----------------------------------------------------------------\n')
            print((prop).ljust(25),end='')
            if prop.find('perc')==-1 and prop.find('n_')==-1 and prop.find('pos')==-1 and prop.find('rvir')==-1:
                print("{0:.2f}".format(np.log10(np.nanmedian(data[prop]))),'\t',\
                      "{0:.2f}".format(np.log10(np.nanmedian(data[prop]))-np.log10(np.nanpercentile(data[prop], 10))),'/',\
                      "{0:.2f}".format(np.log10(np.nanpercentile(data[prop], 90))-np.log10(np.nanmedian(data[prop]))),'\t',\
                      "{0:.2f}".format(min(np.log10(data[prop]))),'/',\
                      "{0:.2f}".format(max(np.log10(data[prop]))),\
                      'log('+unit+')')
                                      
            else:
                print("{0:.2f}".format(np.nanmedian(data[prop])),'\t',\
                      "{0:.2f}".format(np.nanmedian(data[prop])-np.nanpercentile(data[prop], 10)),'/',\
                      "{0:.2f}".format(np.nanpercentile(data[prop], 90)-np.nanmedian(data[prop])),\
                      "{0:.2f}".format(min(data[prop])),\
                      "{0:.2f}".format(max(data[prop])),'\t',\
                      unit)   
        
        print('##############################################\n')
              
    def calc_statistics(self,
                        data,
                        props,
                        sel_name='',
                        sel_value=None,
                        dt=None,
                        print_on_screen=True):
        """Function calculates the median and 10th, 32the, 68th and 90th percentiles. The median and 10th and 90th are also printed on the screen.
            A data-bucket as an 1-row structrued array is returned where the statistics of all properties defined in 'props' are calculated.
            Optionally one can select a subsample from a specific property defined as 'sel_name' and the corresponding value 'sel_value' e.g. snapshot
            column name 'snapid' and snsnapshot number '258' as value.            
        """
        from statsmodels import robust
        #print(np.info(data))
        #print(data[props])
        #print(dt)
        #exit()
        #try:
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
        
        prop_info=ha_lib.get_property_dict()
        
        bucket_stats['n_count'] = sample.size                  
        if print_on_screen==True:
            print(('prop').ljust(25), 'median 10th / 90th\t[unit]\n-----------------------------------------------------------------\n')
        for prop in props:
            unit, myformat = ha_lib.find_unit(prop)

            if print_on_screen==True:
                print((prop).ljust(25), end='')
            if prop_info[prop]['output_prop_as']=='log':
                if print_on_screen==True:
                    print("{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))), '\t',\
                          "{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))-np.log10(np.nanpercentile(sample[prop], 10))),\
                          '/', "{0:.2f}".format(np.log10(np.nanpercentile(sample[prop], 90))-np.log10(np.nanmedian(sample[prop]))), 'log('+unit+')')
                
                bucket_stats[prop+'_log50'] = "{0:.2f}".format(np.log10(np.nanmedian(sample[prop])))
                bucket_stats[prop+'_log10'] = "{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))-np.log10(np.nanpercentile(sample[prop], 10)))
                bucket_stats[prop+'_log32'] = "{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))-np.log10(np.nanpercentile(sample[prop], 32)))
                bucket_stats[prop+'_log68'] = "{0:.2f}".format(np.log10(np.nanpercentile(sample[prop], 68))-np.log10(np.nanmedian(sample[prop])))
                bucket_stats[prop+'_log90'] = "{0:.2f}".format(np.log10(np.nanpercentile(sample[prop], 90))-np.log10(np.nanmedian(sample[prop])))
                bucket_stats[prop+'_logMAD'] = "{0:.2f}".format(np.log10(robust.mad(sample[prop])))                    
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
                bucket_stats[prop+'_MAD'] = robust.mad(sample[prop])

#        except:
#            print('data selection not available! --> SKIPPED!')

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

    def filter_trees(self,
                     snapid=0,
                     create_bins=False,
                     prop='mhalo',
                     bins=[],
                     sample_endfix=''):
        """functions helps in finding, arranging, and storing trees and main branches of those trees into file and returns a dictionary of all trees processes
        where the keys are the top-node halo IDs (rootIndex) at a certain snapshot chosen
        
        There is the possibility to select halos in certain bins and trac the trees for that sample of halos seperately. Set keyword 'bin_data' to 'True' for that.
        
        keywords:
        ------------
            snapid         if subsamples should selected, then snapid sets the snapshot number where the top-node halos should be selected
            create_bins    default is 'False', set 'True' to select subsamples of a certain property with funtion: self.find_subsample()
            prop           property which should be binned
            bins           list of bin edges
            sample_endfix  default is '' string which specify the filename to be printed
            
        output:
        ------------    
            tree_data_bin   dictionary of structured arrays        
        
        """
        
        #load a certain Consitent-tree file
        data = ha_lc.load_ASCII_ConsitentTrees(self.pathname_dict['CT_merger_trees'])
       
        tree_bin_dict={}
        if create_bins==True:
            #choose a subsample if needed
            data_sample=data[np.where(data['snapid']==snapid-self.config_dict['snapid_CT_offset'])[:][0]]
            
            #bin the data in the subsample
            data_bin = self.find_subsample(data_sample,
                                           prop,
                                           bins)              
            
            #get merger trees for halos in certain bins
            for k, mybin in enumerate(data_bin.keys()):
                print('mybin:', mybin, 'unique IDs: ',end='')
                unique_mysample_rootIDs, index, count = np.unique(data_bin[mybin], return_index=True, return_counts=True)
                print(len(unique_mysample_rootIDs),'\n ...processing!')

                
                #Get all tree branches of a certain bin
                tree_data_bin = self.get_progenitor_tree_branches_CT(data,
                                                                     rootIndex=unique_mysample_rootIDs,
                                                                     print_rootIDs_sep=True,
                                                                     print_all_trees=True,
                                                                     print_main_branch=True,
                                                                     sample_endfix='_bin'+str(k))
                
                
                
                tree_bin_dict.update({k: tree_data_bin})
                
            tree_data=tree_bin_dict

        else:
            
            unique_mysample_rootIDs, index, count = np.unique(data['rootIndex'], return_index=True, return_counts=True)

            #Get all tree branches
            tree_data = self.get_progenitor_tree_branches_CT(data,
                                                             rootIndex=[],
                                                             print_rootIDs_sep=True,
                                                             print_all_trees=True,
                                                             print_main_branch=True,
                                                             sample_endfix=sample_endfix)
                

        return tree_data

    def find_subsample(self,
                       data,
                       prop,
                       bins,
                       all_last_gt=True):
        """function selects top-node halo ID ('rootIndex') from a structured array where some condition given by the list bins=[] by certain property 'prop' (e.g. mhalo) 
        is fullfilled. The IDs are stored into a dictionary with the keys corresponding to the accending bin number (starting with zero).

        input:
        ------------    
            data    structured array
            prop    the property of a halo which should be used to for binning (e.g. mhalo)
            bins    a list of values which should be serve a bin edges e.g. [1e11,1e12,1e13,1e14]
            
        keywords:
        ------------    
            all_last_gt     default: True, means that an additional bin is created selecting all IDs where the halos properties exceeds the value of the last bin '>'
        
        output:
        ------------    
            dict_bin    a dictionary with the bin-numbers as keys e.g. dict_bin={0: [ID1, ID1 ... IDn]}
        
        """
        dict_bin={}
        print('create', len(bins), end='')
        if all_last_gt==True:
            print('+1',end='')

        print(' bins with halo property:', prop,'\n+++++++++++++++++++++++++++++++++++++')
            
        for i, mybin in enumerate(bins[:-1]):
            print('bin'+str(i)+': >='+str(mybin)+' & <'+str(bins[i+1]), '\tn_halos: ', end='')
            dict_bin.update({i: data[np.where((data[prop]>=float(mybin)) & (data[prop]<float(bins[i+1])))[:][0]]['rootIndex']})
            print(len(dict_bin[i]))
    
        if all_last_gt==True:
            print('bin'+str(i+1)+': >'+str(bins[-1]),'\t\tn_halos: ', end='')
            dict_bin.update({i+1: data[np.where(data[prop]>=float(bins[-1]))[:][0]]['rootIndex']})
            print(len(dict_bin[i+1]))
        
        return dict_bin
 

    def generate_merger_tree_dict(self,
                                  unique_indices,
                                  data):
        """Function stores all halos on the main branch with there snapids of first and last occurence of a certain tree root to a dictionary. The
        keys of the dictionary are the 'rootIndex'
        
        input:
        ------------                
            unique_indices  list of the unique indices of the rootIds which should be found
            data            a certain catalog
            
        output:
        ------------                
            merger_tree_dict    a dictionary with the 'rootIndex' as key which stores {roodIndex: (first snapid,last snapid, all haloid on the main brainch)}
        """
        merger_tree_dict={}
        
        for rootID in unique_indices:

            data_tree = data[np.where(data['rootIndex']==rootID)[:][0]]                
            data_tree.sort(order=['snapid'], axis=0)
            
            merger_tree_dict.update({rootID: (data_tree['snapid'][0], data_tree['snapid'][-1], data_tree['haloid'])})
        
        return merger_tree_dict         

    def create_dict_unique_IDs(self,
                               data,
                               prop):
        """create a dictionary where the keys are unique IDs of a certain property (e.g. rootIndex, subTreeID, haloid etc.)
        
        input:
        ------------                   
            data:   a structured data array
            prop:   name of the index the data file should be filtered by and put into a dictionary
            
        output:
        ------------                   
            data_dict: dictionary with the unique IDs from prop as keys                     
        """
        data_dict={}
        unique_IDs, index = np.unique(data[prop], return_index=True)
        
        #That line perserves the order in data when the unique indices are found, otherwise, the unique_IDs are sorted (handy e.g. if
        #first 10 most massive main branches should be found etc.)
        unique_IDs=unique_IDs[np.argsort(index)] 

        #mask = np.asarray([random.sample(range(0,unique_IDs.size,1),3)])
   
        #unique_IDs=unique_IDs[mask[:][0]]
        print('Arrange data into a dictionary with unique IDs by -->', prop, 'total amount of unique IDs is:', len(unique_IDs))
        
        for index in unique_IDs: 
            #print('index:', index)                                       
            data_tree = data[np.where(data[prop]==index)[:][0]]
            data_tree.sort(order='a', axis=0)
            
            #print(data_tree[['rootIndex','a','haloid_CT','mhalo']])
            data_dict.update({index: data_tree})
            
        return data_dict

    def get_merger_tree_stats(self,
                              tree_data,
                              print_into_file=True,
                              filename_endfix='',
                              colnames=[]):
        """Function provides statistics such as the median value of the population of halos of a certain halo property
        within uncertainties for each snapshot and prints it to a file
        
        input:
        ------------
            tree_data            structured array of a certain tree
            
        keywords:
        ------------
             print_into_file     default is TRUE, prints statistics of all halo in a merger tree into a file
             filename_endfix     specify the filename to be printed
             
        output:
            
            merger_tree_stats   a certain predefined array of statistic properties of a merger tree or subsample of trees at each snapshot number 
        """              

        merger_tree_stats, myheader, mycols, myformat_string, mydt = ha_lib.create_data_array('stats_custom_bucket',colnames=colnames)            
        #print(mydt)

        snapid = min(tree_data['snapid'])
        snapid_end = max(tree_data['snapid'])
        
        print('snapid:', snapid, 'snapid_end:', snapid_end, 'colsnames:', colnames)
        
        for i, snap in enumerate(range(snapid,snapid_end+1,1)):
            
            data_stats = self.calc_statistics(tree_data, 
                                              colnames, 
                                              sel_name='snapid', 
                                              sel_value=snap, 
                                              dt=mydt,
                                              print_on_screen=False)                

            data_stats['snapid'] = snap
            data_stats['z']      = self.snapid_array[np.where(self.snapid_array['snapid']==snap)[:][0]]['z']
            data_stats['a']      = self.snapid_array[np.where(self.snapid_array['snapid']==snap)[:][0]]['a']
            #print(data_stats)
            
            if i==0:
                merger_tree_stats = data_stats
            else:                                     
                merger_tree_stats = np.append(merger_tree_stats, data_stats, axis=0)
              
        if print_into_file==True:
            ha_lib.writeIntoFile(
                   self.pathname_dict['output_filename_merger_trees'][:-4]+filename_endfix+'_stats.txt',
                   merger_tree_stats[mycols],
                   myheader=self.config_dict['simulation_name']+' '+\
                           str(self.config_dict['box_size_h-1Mpc'])+'h-1Mpc, '\
                           'halo finder: '+self.config_dict['halo_finder']+', '\
                           'tree builder: '+self.config_dict['tree_builder']+', '\
                           ' statistics\n'+myheader,
                   mydelimiter='\t',
                   data_format=myformat_string)
        
        return merger_tree_stats


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

    def get_progenitor_tree_branches_CT(self,
                                       data,
                                       rootIndex=[],
                                       print_all_trees=True,
                                       print_rootIDs_sep=False,
                                       print_main_branch=False,
                                       sample_endfix=''):
        """Function uses Consistent-Tree data run on ROCKSTAR halo catalogs and gets the progenitor tree branches for each top node halo by using the
        a sorting alorithm because the main branch is loacted on the lowest 'DFirstIDs'
        
        input:
        ------------                
            data        structured array (e.g. load CONSITENT TREES standard files 'self.pathname_dict['CT_merger_trees']' and plug into function)
            rootIndex   list of unique indices, default is [], if default then all unique top-node root indices are going to be processed
            
        keywords:
        ------------                
            print_all_trees     default=True --> the halos are stored to a file with filname default: 'self.pathname_dict['path_basename']+CHOLLA_[self.config_dict['res']]_[self.config_dict['box_size_h-1Mpc']]Mpc_MB-CT+[endfix].txt'
            print_rootIDs_sep   default=False --> each tree with the top-node halo ID 'rootIndex' is going to be printed in a seperate text-file with ending '.txt',
                                 the file name can be customized with 'sample_endfix' or/and 'filename_endfix, same name convention as 'print_all_trees'
            print_main_branch   default=False --> each main proginitor tree with the top-node halo ID 'rootIndex' is going to be printed in a seperate text-file with ending '.txt',
                                 the file name can be customized with 'sample_endfix' or/and 'filename_endfix,  same name convention as 'print_all_trees'                                  
            sample_endfix       default='' --> endfix on filename before '.txt'

        output:
        ------------            
            data_CT_MB      structured array

        """

        if rootIndex!=[]:
            #if only certain 'rootIndex' are interesting
            unique_rootIndices, index1, index2 = np.intersect1d(data['rootIndex'], rootIndex, return_indices=True)
            #print statistics of the subsample on screen
            self.give_stats_info(data[index1],props=['mhalo'])
            #print(len(unique_rootIndices), unique_rootIndices)
        else:
            #Find unique rootIndicies (top node of the merger tree)
            unique_rootIndices, index = np.unique(data['rootIndex'], return_index=True)              
        
        data.sort(order=['rootIndex','DFirstID'], axis=0)
        #print('sort by rootIndex & Last main leaf DFirst ID desc -->', data[['snapid','DFirstID', 'LastMLDFirstID', 'rootIndex']][0:10])              

        data = rcfuncs.append_fields([data], 
                                    ['subTreeID'], 
                                    [np.zeros(data.size,)], 
                                    dtypes=['i8'], 
                                    usemask=False)

       
        for a,rootID in enumerate(unique_rootIndices):
            
            data_rootID=data[np.where(data['rootIndex']==rootID)[0][:]]
            print('rootID:', rootID)
            unique_mainLeafIDIndices, index = np.unique(data_rootID['LastMLDFirstID'], return_index=True)
            for k, mainLeafID in enumerate(unique_mainLeafIDIndices):                            
                #print('subtree:', k, 'mailLeafID:', mainLeafID)
                data_subtree=data_rootID[np.where(data_rootID['LastMLDFirstID']==mainLeafID)[0][:]]

                data_subtree['subTreeID']=k
                #print(data_subtree[['snapid','DFirstID', 'LastMLDFirstID', 'rootIndex', 'subTreeID']])
                
                if k==0:
                    data_tree = data_subtree
                else:                                     
                    data_tree = np.append(data_tree, data_subtree, axis=0)                 
              
            filename_endfix='rootID'+str(rootID)
            data_tree['snapid']+=self.config_dict['snapid_CT_offset']
            if print_rootIDs_sep==True:
                #write main branch of each rootIndex into file
                ha_lib.writeIntoFile(
                           self.pathname_dict['output_filename_merger_trees'][:-4]+sample_endfix+'_'+filename_endfix+'.txt',
                           data_tree[['snapid', 'a', 'a_desc', 'a_lastMM', 'rootIndex','descID','DFirstID', 'LastDFirstID', 'LastMLDFirstID','haloid_CT', 'haloid',\
                                    'x_pos', 'y_pos', 'z_pos', 'mhalo', 'vmax', 'rhalf_mass', 'rvir', 'T_U', 'subTreeID']],
                           myheader='snapid(1) a(2) a_desc(3) a_lastMM(4) rootIndex(5) descID(6) DFirstID(7) LastDFirstID(8) LastMLDFirstID(9) haloid_CT(10) orignal_haloid_RS(11) '\
                                     'x_pos(12) y_pos(13) z_pos(14) mhalo(15) vmax(16) rhalf_mass(17) rvir(18) T_U(19) subTreeID(20)',
                           mydelimiter='\t',
                           data_format='%i\t%0.5f\t%0.5f\t%0.5f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%0.5f\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%i') 

            if a==0:
                data_output = data_tree
            else:                                     
                data_output = np.append(data_output, data_tree, axis=0) 

        if print_all_trees==True:

            ha_lib.writeIntoFile(
                               self.pathname_dict['output_filename_merger_trees'][:-4]+sample_endfix+'.txt',
                               data_output[['snapid', 'a', 'a_desc', 'a_lastMM', 'rootIndex','descID','DFirstID', 'LastDFirstID', 'LastMLDFirstID','haloid_CT', 'haloid',\
                                        'x_pos', 'y_pos', 'z_pos', 'mhalo', 'vmax', 'rhalf_mass', 'rvir', 'T_U', 'subTreeID']],
                               myheader='snapid(1) a(2) a_desc(3) a_lastMM(4) rootIndex(5) descID(6) DFirstID(7) LastDFirstID(8) LastMLDFirstID(9) haloid_CT(10) orignal_haloid_RS(11) '\
                                         'x_pos(12) y_pos(13) z_pos(14) mhalo(15) vmax(16) rhalf_mass(17) rvir(18) T_U(19) subTreeID(20)',
                               mydelimiter='\t',
                               data_format='%i\t%0.5f\t%0.5f\t%0.5f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%0.5f\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%i')



        
        if print_main_branch==True:
            data_output_MB=data_output[np.where(data_output['subTreeID']==0)[0][:]]
            #wirte alle main branches into file
            ha_lib.writeIntoFile(
                               self.pathname_dict['output_filename_merger_trees'][:-4]+sample_endfix+'_MB.txt',
                               data_output_MB[['snapid', 'a', 'a_desc', 'a_lastMM', 'rootIndex','descID','DFirstID', 'LastDFirstID', 'LastMLDFirstID','haloid_CT', 'haloid',\
                                        'x_pos', 'y_pos', 'z_pos', 'mhalo', 'vmax', 'rhalf_mass', 'rvir', 'T_U', 'subTreeID']],
                               myheader='snapid(1) a(2) a_desc(3) a_lastMM(4) rootIndex(5) descID(6) DFirstID(7) LastDFirstID(8) LastMLDFirstID(9) haloid_CT(10) orignal_haloid_RS(11) '\
                                         'x_pos(12) y_pos(13) z_pos(14) mhalo(15) vmax(16) rhalf_mass(17) rvir(18) T_U(19) subTreeID(20)',
                               mydelimiter='\t',
                               data_format='%i\t%0.5f\t%0.5f\t%0.5f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%0.5f\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%i')

            
                       
        return data_output



    def load_data_GADGET(self,
                          path, 
                          snapid_start, 
                          snapid_end=None, 
                          return_particle_info=True):
        """function reads and returns the GADGET halo catalog and particle information in binary and optional ASCII file at a given snapshot
            UNDER CONSTRUCTION
        """
        exit()

    def load_data_merger_trees_from_text(self,
                                        rootID,
                                        rootID_end,
                                        file_basename):
        """Function reads and returns the merger trees catalog from ASCII file for a given filename and between certain tree ID numbers"""        
        
        merger_tree_data, myheader, mycols, myformat_string, mydt = ha_lib.create_data_array(dt_name='merger_trees_ASCII')
        
        for i, tree_id in enumerate(range(rootID,rootID_end,1)):
            try:
                data = ha_lib.df_to_sarray(pd.read_csv(self.pathname_dict['path_basename']+file_basename+'_rootID'+str(tree_id)+'.txt', 
                                                       skiprows=2,
                                                       names=mycols,
                                                       sep='\t'))                    
            except:
                print('data for merger tree with rootID:', rootID, 'not available! --> SKIPPED!')
            
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
                                                                                 n_files=self.config_dict['n_files_snapshot'],
                                                                                 myprops=self.config_dict['myprops2read'])
    
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
     
    def construct_main_prog_tree(self,
                              haloid=None,
                              n_trees=1,
                              snapid=None,
                              snapid_end=None,
                              print_merger_tree2file=True,
                              show_infos=True,
                              identify_MB_by='dFirstID'):
        
        """Function constructs merger trees of detected halos using ROCKSTAR halo binary information and descendent IDs 'DescID' from out_x.list
        catalogs (where x stands for the snapshot number)
        
        The main prognitor haloids are printed into an ASCCII file along with the redshift of the first and last detection of that halo and a unique ID
        of the tree
        
        output:
        ===============
            data            structured array of all halo properties of all halo found in the box
            merger_dict     dictionary of the merger tree of each unique tree ID (rootIndex)
            z_first_app     dictionary of the first appearance of 
            
            printed into file each row give one tree in the format ASCII --> z(first appearance) z(last appearance) rootIndex [all haloids of the main progenitors]
        """
        merger_dict={}
        z_first_app={}
        z_last_app={}        

                        
        if haloid==None and n_trees==[]:
             tree_list2analyse=range(0,1,1)           
#        elif haloid!=None:           
#            if type(haloid)==list:
#                tree_list2analyse=haloid
#            else:
#                tree_list2analyse=[haloid]
        else:                
            tree_list2analyse=range(0,n_trees,1)
            
        if snapid!=None:
            mysnap_array=self.snapid_array[np.where(self.snapid_array['snapid']==snapid_end)[:][0][0]:np.where(self.snapid_array['snapid']==snapid)[:][0][0]+1]
        else:
            mysnap_array=self.snapid_array
        
        data_before = self.load_data_ROCKSTAR(self.pathname_dict['path_basename'],
                                              np.sort(mysnap_array['snapid'])[0],
                                              return_particle_info=False)
        if haloid!=None:
            data_before=data_before[np.where(data_before['haloid']==haloid)[:][0]]
            print(data_before[['haloid','descIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos']])

        data_before, data_dummy = self.assign_progenitor_indices(data_before, data_before, np.sort(mysnap_array['snapid'])[0], first_assignment=0)
        print(data_before[['haloid','descIndex', 'rootIndex', 'predIndex', 'mhalo', 'delta_mhalo', 'rvir', 'delta_rvir', 'x_pos', 'y_pos', 'z_pos', 'snapid', 'n_particles', 'npros']], '\n', data_dummy)

        for tree_id in tree_list2analyse:
            progs_tree=[]
            print('tree_id:', tree_id, '/', len(tree_list2analyse)-1, '... processing!')
                
            count_first_app=0                 
            count_last_app=0            
            for i, snapid in enumerate(np.sort(mysnap_array['snapid'])[1::]):
                
                if i==0:
                    data=data_before
 
                else:
                    data = np.append(data, data_before, axis=0)                
                print('BEFORE desc indentification:\n', data_before[['haloid','descIndex', 'rootIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos']])
                new_data, progID, data_before = self.get_MMP_data(self.pathname_dict['path_basename'],
                                                                 data_before,
                                                                 snapid,
                                                                 tree_id,
                                                                 sort_prop_by=identify_MB_by)
 
                if show_infos==True:
                    print('processing treeID:', tree_id, '--> progIDs found at snapshot', snapid,':', progID,\
                          #'data_before:', data_before[['haloid','descIndex', 'rootIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos', 'predIndex']],\
                          '\nnew_data:', new_data[['haloid','descIndex', 'rootIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos', 'predIndex']])

    
    
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
            filename_trees=self.pathname_dict['path_basename']+'merger_tree_tracks.txt'    

            for tree_id in merger_dict.keys():
                if tree_id==tree_list2analyse[0]:
                    ha_lib.writeIntoFile(
                               filename_trees,
                               [str(z_first_app[tree_id])+' '+str(z_last_app[tree_id])+' '+str(tree_id)+' '+" ".join(str(item) for item in merger_dict[tree_id][1:])],
                               myheader=str(self.pathname_dict['path_basename'])+'Main progenitor merger trees as haloid (ID in the orginal ROCKSTAR files)\n'\
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
        
        return data_tree, merger_dict, z_first_app
    
    
    def get_MMP_data(self,
                            path,
                            data_before,
                            snapid,
                            tree_id,
                            first_assignment=1,
                            show_infos=True,
                            sort_by_prop='dFirstID'):
        """function access the progenitor information finds the correct (most massive) progenitor index to a given haloid"""
        print('HERE 721! data_before:\n', data_before[['haloid','descIndex', 'predIndex', 'rootIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos']])

        data = self.load_data_ROCKSTAR(self.pathname_dict['path_basename'],
                                       snapid,
                                       return_particle_info=False)
            
        data, data_before = self.assign_progenitor_indices(data, data_before, snapid, first_assignment=first_assignment)
        print('HERE 728! data_before:\n', data_before[['haloid','descIndex', 'predIndex', 'rootIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos']])       
        tree_data=data[np.where(data['rootIndex']==tree_id)[:][0]]
        if show_infos==True:
            print('HERE 731! tree_data:\n',tree_data[['haloid','descIndex', 'predIndex', 'rootIndex', 'n_particles', 'mhalo', 'x_pos', 'y_pos', 'z_pos']])
            
        tree_data[::-1].sort(order=sort_by_prop, axis=0)
    
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

    def track_merger_tree(self, snapid_first, snapid_last, rootID, progIDs, if_CT_data=True):

        data, myheader, mycols, myformat_string, mydt = ha_lib.create_data_array(dt_name='merger_trees_ASCII')        

        for i, snapid in enumerate(range(snapid_first,snapid_last,1)):

            if if_CT_data==True:
                print('halo1:', progIDs[i], 'SN1:', snapid, 'halo2:', progIDs[i+1], 'SN2:', snapid+1)
                halo1 = self.myDataHalo[np.where((self.myDataHalo['haloid']==progIDs[i]) & (self.myDataHalo['snapid']==snapid))[:][0]]
                halo2 = self.myDataHalo[np.where((self.myDataHalo['haloid']==progIDs[i+1]) & (self.myDataHalo['snapid']==snapid+1))[:][0]]
            else:
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
                print('SKIPPED!')
                #return False
        
        if len(progIDs)>0:         
            ha_lib.writeIntoFile(
                       self.pathname_dict['path_basename']+self.pathname_dict['output_filename_data_tree_basename']+str(rootID)+'.txt',
                       data[mycols],
                       myheader='Cholla'+str(self.config_dict['res'])+' '+str(self.config_dict['box_size_h-1Mpc'])+'h-1Mpc all halos in maind branch with rootID: '+str(rootID)+'\n'+myheader,
                       mydelimiter='\t',
                       data_format=myformat_string)      
            
        return data
       