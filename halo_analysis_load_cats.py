"""halo_analysis_load_cats.py

This file contains fuctions which load various halo catalogs later to be analysed by halo_analysis.py

"""
from __future__ import print_function
# Load packages
#own
import halo_analysis_lib as ha_lib

#system
import pandas as pd
import numpy as np

def get_dtype_sarray(name):
    """Function returns the a dictionary of the labels and dtypes of the properties stored in certain output files such as ROCKSTAR or Gadget.
        cfc stands for 'common filename convention'
    
        input
        =========
            
        name:   name of the output files such as ROCKSTAR_ASCII (cfc: halos_x.x.ascii), 
                                                 ROCKSTAR_binary (cfc: halos_x.x.bin), 
                                                 ROCKSTAR_ASCII_list (cfc: out_x.list files whole halo catalog at snapshot)
                                                 Gadget_binary (cfc: ?)               
        output
        =========
        
        dt:     dictionary of {(label: dtype)}
                usage e.g. mydata=np.fromfile(myfile, dt)
    
    """
    
    def ConsistentTrees_ASCII_099():
        """
        Property names and column in Consitent-Tree ascii (CT) files (after runing it on ROCKSTAR), original labels of properties in the file are commanded. 
        Label names are converted to internal universal name convention. This is the bucket for version 0.99 of CT

       
        List of property labels in Rockstar ascii files:    
            scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14)
            scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29)
            Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Rs_Klypin Mvir_all M200b M200c M500c M2500c Xoff Voff
            Spin_Bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius
        """
        dt=    [('a'               , np.float32), #scale(0)
                ('haloid_CT'       , np.int64),   #id(1)
                ('a_desc'          , np.float32), #desc_scale(2)
                ('descID'          , np.int64),   #desc_id(3)
                ('n_prog'          , np.int32),   #num_prog(4)
                ('hostid_LM'       , np.int64),   #pid(5)              
                ('hostid_MM'       , np.int64),   #upid(6)
                ('desc_hostid_MM'  , np.int64),   #desc_pid(7)
                ('is_phantom'      , np.int32),    #phantom(8)             
                ('sam_mvir'        , np.float32), #sam_mvir(9)
                ('mhalo'           , np.float32), #mvir(10)                
                ('rvir'            , np.float32), #rvir(11)
                ('rscale'          , np.float32), #rs(12)
                ('vrms'            , np.float32), #vrms(13)                               
                ('is_mmp'          , np.int32),    #mmp?(14)
                ('a_lastMM'        , np.float32), #scale_of_last_MM(15)
                ('vmax'            , np.float32), #vmax(16)
                ('x_pos'           , np.float32), #x(17)
                ('y_pos'           , np.float32), #y(18)
                ('z_pos'           , np.float32), #z(19)
                ('x_vel'           , np.float32), #vx(20)
                ('y_vel'           , np.float32), #vy(21)
                ('z_vel'           , np.float32), #vz(22)
                ('x_ang'           , np.float32), #Jx(23)
                ('y_ang'           , np.float32), #Jy(24)
                ('z_ang'           , np.float32), #Jz(25)
                ('spinParameter'   , np.float32), #Spin(26)
                ('BFirstID'        , np.int64),   #Breadth_first_ID(27)
                ('DFirstID'        , np.int64),   #Depth_first_ID(28)               
                ('rootIndex'       , np.int64),   #Tree_root_ID(29)                               
                ('haloid'          , np.int64),   #Orig_halo_ID(30)
                ('snapid'          , np.int32),   #Snap_num(31)                
                ('NextCoDFirstID'  , np.int64),   #Next_coprogenitor_depthfirst_ID(32)
                ('LastDFirstID'    , np.int64),   #Last_progenitor_depthfirst_ID(33)                 
                ('rscale_Klypin'   , np.float32), #rs_Klypin
                ('mhalo+unbound'   , np.float32), #mvir_all                            
                ('mhalo_200b'      , np.float32), #m200b   
                ('mhalo_200c'      , np.float32), #m200c
                ('mhalo_500c'      , np.float32), #m500c
                ('mhalo_2500c'     , np.float32), #m2500c
                ('x_off'           , np.float32), #Xoff
                ('v_off'           , np.float32), #Yoff
                ('spin_Bullock'    , np.float32), #spin_bullock
                ('b_to_a'          , np.float32), #b_to_a  
                ('c_to_a'          , np.float32), #c_to_a
                ('x_a'             , np.float32), #A[x]
                ('y_a'             , np.float32), #A[y]
                ('z_a'             , np.float32), #A[z] 
                ('b_to_a_500c'     , np.float32), #b_to_a(500c)
                ('c_to_a_500c'     , np.float32), #c_to_a(500c)    
                ('x_a_500c'        , np.float32), #A[x](500c)    
                ('y_a_500c'        , np.float32), #A[y](500c) 
                ('z_a_500c'        , np.float32), #A[z](500c)
                ('T_U'             , np.float32), #T/|U|
                ('Mpseudo_Behroozi', np.float32), #M_pe_Behroozi
                ('Mpseudo_Diemer'  , np.float32), #M_pe_Diemer
                ('rhalf_mass'      , np.float32), #Halfmass_Radius           
                 ] 
        return dt

    def ConsistentTrees_ASCII_101():
        """
        Property names and column in Consitent-Tree ascii (CT) files (after runing it on ROCKSTAR), original labels of properties in the file are commanded. 
        Label names are converted to internal universal name convention. This is the bucket for version 1.01 of CT. This version has the very import
        mainLeaf_depthFirstIDincluded

       
        List of property labels in Rockstar ascii files:    
            scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14)
            scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29)
            Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Rs_Klypin Mvir_all M200b M200c M500c M2500c Xoff Voff
            Spin_Bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius
            
            
            #scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_Mvir(9) Mvir(10) Rvir(11) rs(12) vrms(13) mmp?(14)
            scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29)
            Orig_halo_ID(30) Snap_idx(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Tidal_Force(35) Tidal_ID(36) 
            Rs_Klypin(37) Mvir_all(38) M200b(39) M200c(40) M500c(41) M2500c(42) Xoff(43) Voff(44) Spin_Bullock(45) b_to_a(46) c_to_a A[x](47) A[y](48) A[z](49) b_to_a(500c)(50) c_to_a(500c)(51) A[x](500c)(52) A[y](500c)(53) A[z](500c)(54) T/|U|(55) M_pe_Behroozi(56) M_pe_Diemer(57) Halfmass_Radius(58)
        """
        dt=    [('a'               , np.float32), #scale(0)
                ('haloid_CT'       , np.int64),   #id(1)
                ('a_desc'          , np.float32), #desc_scale(2)
                ('descID'          , np.int64),   #desc_id(3)
                ('n_prog'          , np.int32),   #num_prog(4)
                ('hostid_LM'       , np.int64),   #pid(5)              
                ('hostid_MM'       , np.int64),   #upid(6)
                ('desc_hostid_MM'  , np.int64),   #desc_pid(7)
                ('is_phantom'      , np.int32),    #phantom(8)             
                ('sam_mvir'        , np.float32), #sam_mvir(9)
                ('mhalo'           , np.float32), #mvir(10)                
                ('rvir'            , np.float32), #rvir(11)
                ('rscale'          , np.float32), #rs(12)
                ('vrms'            , np.float32), #vrms(13)                               
                ('is_mmp'          , np.int32),    #mmp?(14)
                ('a_lastMM'        , np.float32), #scale_of_last_MM(15)
                ('vmax'            , np.float32), #vmax(16)
                ('x_pos'           , np.float32), #x(17)
                ('y_pos'           , np.float32), #y(18)
                ('z_pos'           , np.float32), #z(19)
                ('x_vel'           , np.float32), #vx(20)
                ('y_vel'           , np.float32), #vy(21)
                ('z_vel'           , np.float32), #vz(22)
                ('x_ang'           , np.float32), #Jx(23)
                ('y_ang'           , np.float32), #Jy(24)
                ('z_ang'           , np.float32), #Jz(25)
                ('spinParameter'   , np.float32), #Spin(26)
                ('BFirstID'        , np.int64),   #Breadth_first_ID(27)
                ('DFirstID'        , np.int64),   #Depth_first_ID(28)               
                ('rootIndex'       , np.int64),   #Tree_root_ID(29)                               
                ('haloid'          , np.int64),   #Orig_halo_ID(30)
                ('snapid'          , np.int32),   #Snap_num(31)                
                ('NextCoDFirstID'  , np.int64),   #Next_coprogenitor_depthfirst_ID(32)
                ('LastDFirstID'    , np.int64),   #Last_progenitor_depthfirst_ID(33)
                ('LastMLDFirstID', np.int64),#Last_mainleaf_depthfirst_ID(34)
                ('tidalForce'      , np.float32), #Tidal_Force
                ('tidalID'         , np.int64), #Tidal_ID                 
                ('rscale_Klypin'   , np.float32), #rs_Klypin
                ('mhalo+unbound'   , np.float32), #mvir_all                            
                ('mhalo_200b'      , np.float32), #m200b   
                ('mhalo_200c'      , np.float32), #m200c
                ('mhalo_500c'      , np.float32), #m500c
                ('mhalo_2500c'     , np.float32), #m2500c
                ('x_off'           , np.float32), #Xoff
                ('v_off'           , np.float32), #Yoff
                ('spin_Bullock'    , np.float32), #spin_bullock
                ('b_to_a'          , np.float32), #b_to_a  
                ('c_to_a'          , np.float32), #c_to_a
                ('x_a'             , np.float32), #A[x]
                ('y_a'             , np.float32), #A[y]
                ('z_a'             , np.float32), #A[z] 
                ('b_to_a_500c'     , np.float32), #b_to_a(500c)
                ('c_to_a_500c'     , np.float32), #c_to_a(500c)    
                ('x_a_500c'        , np.float32), #A[x](500c)    
                ('y_a_500c'        , np.float32), #A[y](500c) 
                ('z_a_500c'        , np.float32), #A[z](500c)
                ('T_U'             , np.float32), #T/|U|
                ('Mpseudo_Behroozi', np.float32), #M_pe_Behroozi
                ('Mpseudo_Diemer'  , np.float32), #M_pe_Diemer
                ('rhalf_mass'      , np.float32), #Halfmass_Radius           
                 ] 
        return dt

    

    def ConsistentTrees_basic_ASCII():
        """
        Property names and column in Consitent-Tree ascii files (after runing it on ROCKSTAR), original labels of properties in the file are commanded. 
        Label names are converted to internal universal name convention.

       
        Basic list of property labels in ConsistentTree ascii files:    
            # snapid(1) a(2) a_desc(3) a_lastMM(4) rootIndex(5) descID(6) DFirstID(7) LastDFirstID(8) LastMLDFirstID(9) haloid_CT(10) orignal_haloid_RS(11) '
            x_pos(12) y_pos(13) z_pos(14) mhalo(15) vmax(16) rhalf_mass(17) rvir(18) T_U(19) subTreeID(20)
        """
        dt=    [('snapid'          , np.int32), #(1) haloid 
                ('a'          , np.float32),   #(2) descIndex
                ('a_desc'          , np.float32), #(3) rootIndex 
                ('a_lastMM'          , np.float32),   #(4) predIndex
                ('rootIndex'           , np.int64), #(5) Mvir [h-1Msun]                
                ('descID'     , np.int64), #(6) delta_Mvir [h-1Msun]
                ('DFirstID'            , np.int64), #(7) Rvir [h-1kpc]
                ('LastDFirstID'      , np.int64), #(8) delt_Rvir [h-1kpc]              
                ('LastMLDFirstID'          , np.int64), #(12) snapid
                ('haloid_CT'     , np.int64), #(13) n_particles
                ('orignal_haloid_RS'         , np.int64), #(14) n_progs
                ('x_pos'           , np.float32), #(9) X [h-1Mpc]
                ('y_pos'           , np.float32), #(10) Y [h-1Mpc]
                ('z_pos'           , np.float32), #(11) Z [h-1Mpc]
                ('mhalo'           , np.float32), #(9) X [h-1Mpc]
                ('vmax'           , np.float32), #(10) Y [h-1Mpc]
                ('rhalf_mass'           , np.float32), #(11) Z [h-1Mpc]
                ('rvir'           , np.float32), #(9) X [h-1Mpc]
                ('T_U'           , np.float32), #(10) Y [h-1Mpc]
                ('subTreeID'           , np.int32), #(11) Z [h-1Mpc]                
                 ] 
        return dt

    
    

    def stats_basic_bucket():

        dt =    [
                ('mhalo1'                      , np.float32),
                ('delta_mhalo'                 , np.float32),
                ('delta_mhalo_perc'            , np.float32),                
                ('rvir1'                       , np.float32),
                ('delta_rvir'                  , np.float32),
                ('delta_rvir_perc'             , np.float32),                
                ('n_particles1'                , np.float32),
                ('n_particles2'                , np.float32),
                ('n_particles_shared'          , np.float32),
                ('n_particles_shared_perc1'    , np.float32), 
                ('n_particles_shared_perc2'    , np.float32),
                ('delta_x_pos_perc'            , np.float32), 
                ('delta_y_pos_perc'            , np.float32), 
                ('delta_z_pos_perc'            , np.float32)
                ]
        
        return dt
        
    def stats_perc_bucket():
      
        dt =    [
                ('snapid'                      , np.int32),  
                ('z'                           , np.float32),                
                ('mhalo1_log50'                , np.float32),
                ('mhalo1_log10'                , np.float32),
                ('mhalo1_log32'                , np.float32),                 
                ('mhalo1_log68'                   , np.float32),
                ('mhalo1_log90'                   , np.float32),                                 
                ('delta_mhalo_log50'              , np.float32),
                ('delta_mhalo_log10'              , np.float32),
                ('delta_mhalo_log32'              , np.float32),                 
                ('delta_mhalo_log68'              , np.float32),
                ('delta_mhalo_log90'              , np.float32),               
                ('delta_mhalo_perc_50'         , np.float32),
                ('delta_mhalo_perc_10'         , np.float32),
                ('delta_mhalo_perc_32'         , np.float32),                 
                ('delta_mhalo_perc_68'         , np.float32),
                ('delta_mhalo_perc_90'         , np.float32),                
                ('rvir1_50'                    , np.float32),
                ('rvir1_10'                    , np.float32),
                ('rvir1_32'                    , np.float32),                 
                ('rvir1_68'                    , np.float32),
                ('rvir1_90'                    , np.float32),               
                ('delta_rvir_50'               , np.float32),
                ('delta_rvir_10'               , np.float32),
                ('delta_rvir_32'               , np.float32),                 
                ('delta_rvir_68'               , np.float32),
                ('delta_rvir_90'               , np.float32),              
                ('delta_rvir_perc_50'          , np.float32),
                ('delta_rvir_perc_10'          , np.float32),
                ('delta_rvir_perc_32'          , np.float32),                 
                ('delta_rvir_perc_68'          , np.float32),
                ('delta_rvir_perc_90'          , np.float32),                               
                ('n_particles1_50'             , np.float32),
                ('n_particles1_10'             , np.float32),
                ('n_particles1_32'             , np.float32),                 
                ('n_particles1_68'             , np.float32),
                ('n_particles1_90'             , np.float32),                   
                ('n_particles2_50'             , np.float32),
                ('n_particles2_10'             , np.float32),
                ('n_particles2_32'             , np.float32),                 
                ('n_particles2_68'             , np.float32),
                ('n_particles2_90'             , np.float32),                
                ('n_particles_shared_50'       , np.float32),
                ('n_particles_shared_10'       , np.float32),
                ('n_particles_shared_32'       , np.float32),                 
                ('n_particles_shared_68'       , np.float32),
                ('n_particles_shared_90'       , np.float32),                
                ('n_particles_shared_perc1_50' , np.float32),
                ('n_particles_shared_perc1_10' , np.float32),
                ('n_particles_shared_perc1_32' , np.float32),                 
                ('n_particles_shared_perc1_68' , np.float32),
                ('n_particles_shared_perc1_90' , np.float32),                  
                ('n_particles_shared_perc2_50' , np.float32),
                ('n_particles_shared_perc2_10' , np.float32),
                ('n_particles_shared_perc2_32' , np.float32),                 
                ('n_particles_shared_perc2_68' , np.float32),
                ('n_particles_shared_perc2_90' , np.float32),                      
                ('delta_x_pos_perc_50'         , np.float32),              
                ('delta_x_pos_perc_10'         , np.float32),
                ('delta_x_pos_perc_32'         , np.float32),                 
                ('delta_x_pos_perc_68'         , np.float32),
                ('delta_x_pos_perc_90'         , np.float32),                  
                ('delta_y_pos_perc_50'         , np.float32),                 
                ('delta_y_pos_perc_10'         , np.float32),
                ('delta_y_pos_perc_32'         , np.float32),                 
                ('delta_y_pos_perc_68'         , np.float32),
                ('delta_y_pos_perc_90'         , np.float32),                                                    
                ('delta_z_pos_perc_50'         , np.float32),
                ('delta_z_pos_perc_10'         , np.float32),
                ('delta_z_pos_perc_32'         , np.float32),                 
                ('delta_z_pos_perc_68'         , np.float32),
                ('delta_z_pos_perc_90'         , np.float32),
                ('n_count'                     , np.int64)
                ]
        
        return dt


    def merger_trees_ASCII():

        dt =    [
                ('haloid1'          , np.int64),  
                ('haloid2'          , np.int64),                     
                ('descIndex1'       , np.int64), 
                ('descIndex2'       , np.int64),                  
                ('rootIndex'           , np.int64),
                ('snapid1'            , np.int32), 
                ('snapid2'            , np.int32),  
                ('z1'            , np.float32), 
                ('z2'            , np.float32),                   
                ('mhalo1'            , np.float32),
                ('mhalo2'            , np.float32),
                ('delta_mhalo'      , np.float32),
                ('delta_mhalo_perc'   , np.float32),                    
                ('rvir1'            , np.float32),
                ('rvir2'            , np.float32), 
                ('delta_rvir',        np.float32),
                ('delta_rvir_perc'   , np.float32),                     
                ('n_particles1'      , np.int64),
                ('n_particles2'      , np.int64),
                ('n_particles_shared', np.float32),
                ('n_particles_shared_perc1', np.float32),
                ('n_particles_shared_perc2', np.float32),                      
                ('x_pos1'             , np.float32),
                ('y_pos1'             , np.float32), 
                ('z_pos1'             , np.float32),
                ('x_pos2'             , np.float32),
                ('y_pos2'             , np.float32), 
                ('z_pos2'             , np.float32),                      
                ('delta_x_pos_perc'   , np.float32),                              
                ('delta_y_pos_perc'   , np.float32),                     
                ('delta_z_pos_perc'   , np.float32),
                ('subTreeID'          , np.int64)
                ]
        
        return dt

    
    def ROCKSTAR_ASCII():
        """
        Property names and column in Rockstar ascii files, original labels of properties in the file are commanded. 
        Label names are converted to internal universal name convention.
        
        List of property labels in Rockstar ascii files:    
             id num_p mvir mbound_vir rvir vmax rvmax vrms x y z vx vy vz Jx Jy Jz E Spin PosUncertainty VelUncertainty bulk_vx bulk_vy bulk_vz BulkVelUnc n_core
             m200b m200c m500c m2500c Xoff Voff spin_bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) 
             Rs Rs_Klypin T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius idx i_so i_ph num_cp mmetric
        """
        dt=    [
                ('haloid'          , np.int64),  #id
                ('n_particles'      , np.int64),  #num_p
                ('mhalo'           , np.float32), #mvir                   
                ('mbasic'          , np.float32), #mbound_vir             
                ('rvir'            , np.float32), #rvir
                ('vmax'            , np.float32), #vmax
                ('rvmax'           , np.float32), #rvmax
                ('vrms'            , np.float32), #vrms                  
                ('x_pos'           , np.float32), #x
                ('y_pos'           , np.float32), #y
                ('z_pos'           , np.float32), #z
                ('x_vel'           , np.float32), #vx
                ('y_vel'           , np.float32), #vy
                ('z_vel'           , np.float32), #vz
                ('x_ang'           , np.float32), #Jx
                ('y_ang'           , np.float32), #Jy
                ('z_ang'           , np.float32), #Jz
                ('engery'          , np.float32), #E
                ('spinParameter'   , np.float32), #Spin
                ('unc_pos'         , np.float32), #PosUncertainty
                ('unc_vel'         , np.float32), #VelUncertainty
                ('x_vel_bulk'      , np.float32), #bulx_vx
                ('y_vel_bulk'      , np.float32), #bulx_vy
                ('z_vel_bulk'      , np.float32), #bulx_vz
                ('unc_vel_bulk'    , np.float32), #BulkVelUnc
                ('n_core'          , np.int64),   #n_core
                ('mhalo_200b'      , np.float32), #m200b   
                ('mhalo_200c'      , np.float32), #m200c
                ('mhalo_500c'      , np.float32), #m500c
                ('mhalo_2500c'     , np.float32), #m2500c
                ('x_off'           , np.float32), #Xoff
                ('v_off'           , np.float32), #Yoff
                ('spin_Bullock'    , np.float32), #spin_bullock
                ('b_to_a'          , np.float32), #b_to_a  
                ('c_to_a'          , np.float32), #c_to_a
                ('x_a'             , np.float32), #A[x]
                ('y_a'             , np.float32), #A[y]
                ('z_a'             , np.float32), #A[z] 
                ('b_to_a_500c'     , np.float32), #b_to_a(500c)
                ('c_to_a_500c'     , np.float32), #c_to_a(500c)    
                ('x_a_500c'        , np.float32), #A[x](500c)    
                ('y_a_500c'        , np.float32), #A[y](500c) 
                ('z_a_500c'        , np.float32), #A[z](500c)
                ('rscale'          , np.float32), #rs 
                ('rscale_Klypin'   , np.float32), #rs_Klypin
                ('T_U'             , np.float32), #T/|U|
                ('Mpseudo_Behroozi', np.float32), #M_pe_Behroozi
                ('Mpseudo_Diemer'  , np.float32), #M_pe_Diemer
                ('rhalf_mass'      , np.float32), #Halfmass_Radius
                ('idx'             , np.int64), #internal debugging quantity
                ('i_so'            , np.int64), #internal debugging quantity
                ('i_ph'            , np.int64), #internal debugging quantity
                ('n_particle_child', np.int64), #internal debugging quantity
                ('mmetric'         , np.float32) #internal debugging quantity              
                 ] 
        return dt

    def ROCKSTAR_binary():
        """Function returns two dtype structure one for the halo properties (dt) and on for header information (dt_halo_info --> number of halo, number of particles, particle mass,
            particle type (dark matter, gas, etc.), header size (header_size), bytes of information per halo (halo_struct_size), location of bytes where to start the reading from
            (bytes_to_header_info)      
            
            halo structure of the binary dump from halo.h in ROCKSTAR (https://bitbucket.org/gfcstanford/rockstar/src/main/)
            
                    struct halo {
                      int64_t id;
                      float pos[6], corevel[3], bulkvel[3];
                      float m, r, child_r, vmax_r, mgrav, vmax, rvmax, rs, klypin_rs, vrms,
                        J[3], energy, spin, alt_m[4], Xoff, Voff, b_to_a, c_to_a, A[3],
                        b_to_a2, c_to_a2, A2[3],
                        bullock_spin, kin_to_pot, m_pe_b, m_pe_d, halfmass_radius;
                      int64_t num_p, num_child_particles, p_start, desc, flags, n_core;
                      float min_pos_err, min_vel_err, min_bulkvel_err;
                    };
                    
                    struct extra_halo_info {
                       int64_t child, next_cochild, prev_cochild;
                       int64_t sub_of, ph;
                       float max_metric;
            };
        """
        header_size          = 256 #Bytes, size of the header
        halo_struct_size     = 264 #Bytes, properties stored for one halo using dtype structure dt (260 from struct 'halo' in halo.h from ROCKSTAR and 
                                   #4 bytes probably from max_metric from struct 'extra_halo_info' in halo.h)
        bytes_to_header_info = 64  #bytes until the header info starts
        
        dt_header_info = [ 
                            ('n_halos'          , np.int64), #total number of halos in this file
                            ('tot_n_particles'  , np.int64), #total number of particles in this file                            
                            ('box_size'         , np.float32), #side lenght in Mpc/h of simulation box
                            ('m_particles'      , np.float32), #mass of one particle in h-1Msun
                            ('type_particles'   , np.int64) #type of particle (either 1=halo, star, gas etc.)                            
                            ]
        
        dt =    [
                ('haloid'          , np.int64), #int64_t id
                ('x_pos'           , np.float32), #float pos[6], 1
                ('y_pos'           , np.float32), #float pos[6], 2
                ('z_pos'           , np.float32), #float pos[6], 3
                ('pos4'            , np.float32), #float pos[6], 4
                ('pos5'            , np.float32), #float pos[6], 5
                ('pos6'            , np.float32), #float pos[6], 6               
                ('x_corevel'       , np.float32), #float corevel[3], 1
                ('y_corevel'       , np.float32), #float corevel[3], 2
                ('z_corevel'       , np.float32), #float corevel[3], 3                 
                ('x_vel_bulk'      , np.float32), #float bulkvel[3], 1
                ('y_vel_bulk'      , np.float32), #float bulkvel[3], 2
                ('z_vel_bulk'      , np.float32), #float bulkvel[3], 3
                ('mhalo'           , np.float32), #float m            
                ('rvir'            , np.float32), #float r          
                ('rvir_child'      , np.float32), #float child_r
                ('vmax_r'          , np.float32), #float vmax_r
                ('mhalo_bound'     , np.float32), #float mgrav
                ('vmax'            , np.float32), #float vmax
                ('vpeak'           , np.float32), #float rvmax
                ('rscale'          , np.float32), #float rs
                ('rscale_Klypin'   , np.float32), #float klypin_rs
                ('vrms'            , np.float32), #float vrms
                ('x_ang'           , np.float32), #float J[3], 1
                ('y_ang'           , np.float32), #float J[3], 2
                ('z_ang'           , np.float32), #float J[3], 3
                ('energy'          , np.float32), #float energy                       
                ('spinParameter'   , np.float32), #float spin
                ('mhalo_200b'      , np.float32), #float alt_m[4], 1 
                ('mhalo_200c'      , np.float32), #float alt_m[4], 2 
                ('mhalo_500c'      , np.float32), #float alt_m[4], 3 
                ('mhalo_2500c'     , np.float32), #float alt_m[4], 4 
                ('x_off'           , np.float32), #float Xoff
                ('v_off'           , np.float32), #float Voff
                ('b_to_a'          , np.float32), #float b_to_a  
                ('c_to_a'          , np.float32), #float c_to_a
                ('x_a'             , np.float32), #float A[3], 1
                ('y_a'             , np.float32), #float A[3], 2
                ('z_a'             , np.float32), #float A[3], 3            
                ('b_to_a_500c'     , np.float32), #float b_to_a2
                ('c_to_a_500c'     , np.float32), #float c_to_a2
                ('x_a_500c'        , np.float32), #float A2[3], 1   
                ('y_a_500c'        , np.float32), #float A2[3], 2
                ('z_a_500c'        , np.float32), #float A2[3], 3           
                ('spin_Bullock'    , np.float32), #float bullock_spin
                ('T_U'             , np.float32), #float kin_to_pot
                ('Mpseudo_Behroozi', np.float32), #float m_pe_b 
                ('Mpseudo_Diemer'  , np.float32), #float m_pe_d
                ('rhalf_mass'      , np.float32), #float halfmass_radius
                ('n_particles'      , np.int64), #int64_t num_p
                ('n_particles_child', np.int64), #int64_t num_child_particles         
                ('p_start'         , np.int64), #int64_t p_start
                ('descIndex'       , np.int64), #int64_t desc
                ('flags'           , np.int64), #int64_t flags
                ('n_core'          , np.int64), #int64_t n_core
                ('PosUncertainty'  , np.float32), #float min_pos_err
                ('VelUncertainty'  , np.float32), #float min_vel_err
                ('BulkVelUnc'      , np.float32),  #float min_bulkvel_err
                ('mmetric'          , np.float32)  #unclear where it comes from, it might be mmetric             
                 ]
        
        return header_size, halo_struct_size, dt, dt_header_info, bytes_to_header_info

    def ROCKSTAR_ASCII_list():
        """Property names and column in ROCKSTAR out_x.list files (where x stands for the snapshot number id), original labels of properties in the file are commanded. 
        Label names are converted to internal universal name convention.
        
        In this file descendent information can befound (descIndex)!
        
        List of property labels in ROCKSTAR out_x.list files:        
            ID DescID Mvir Vmax Vrms Rvir Rs Np X Y Z VX VY VZ JX JY JZ Spin rs_klypin Mvir_all M200b M200c M500c M2500c Xoff Voff spin_bullock
            b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius
        """
        
        dt=    [
                ('haloid'          , np.int64),  #ID
                ('descIndex'       , np.int64),  #DescID
                ('mhalo'           , np.float32), #Mvir
                ('vmax'            , np.float32), #Vmax
                ('vrms'            , np.float32), #Vrms                         
                ('rvir'            , np.float32), #Rvir
                ('rscale'          , np.float32), #Rs
                ('n_particles'     , np.int64), #Np                                               
                ('x_pos'           , np.float32), #x
                ('y_pos'           , np.float32), #y
                ('z_pos'           , np.float32), #z
                ('x_vel'           , np.float32), #vx
                ('y_vel'           , np.float32), #vy
                ('z_vel'           , np.float32), #vz
                ('x_ang'           , np.float32), #Jx
                ('y_ang'           , np.float32), #Jy
                ('z_ang'           , np.float32), #Jz
                ('spinParameter'   , np.float32), #Spin
                ('rscale_Klypin'   , np.float32), #rs_Klypin               
                ('mbasic'          , np.float32), #Mvir_all
                ('mhalo_200b'      , np.float32), #m200b   
                ('mhalo_200c'      , np.float32), #m200c
                ('mhalo_500c'      , np.float32), #m500c
                ('mhalo_2500c'     , np.float32), #m2500c
                ('x_off'           , np.float32), #Xoff
                ('v_off'           , np.float32), #Yoff
                ('spin_Bullock'    , np.float32), #spin_bullock
                ('b_to_a'          , np.float32), #b_to_a  
                ('c_to_a'          , np.float32), #c_to_a
                ('x_a'             , np.float32), #A[x]
                ('y_a'             , np.float32), #A[y]
                ('z_a'             , np.float32), #A[z] 
                ('b_to_a_500c'     , np.float32), #b_to_a(500c)
                ('c_to_a_500c'     , np.float32), #c_to_a(500c)    
                ('x_a_500c'        , np.float32), #A[x](500c)    
                ('y_a_500c'        , np.float32), #A[y](500c) 
                ('z_a_500c'        , np.float32), #A[z](500c)
                ('T_U'             , np.float32), #T/|U|
                ('Mpseudo_Behroozi', np.float32), #M_pe_Behroozi
                ('Mpseudo_Diemer'  , np.float32), #M_pe_Diemer
                ('rhalf_mass'      , np.float32)  #Halfmass_Radius         
                 ] 
        return dt

    def Gadget_binary():
        dt={}
        return dt

    choose = {
        'merger_trees_ASCII':   merger_trees_ASCII,
        'ConsistentTrees_basic_ASCII':   ConsistentTrees_basic_ASCII,        
        'ConsistentTrees_ASCII_099': ConsistentTrees_ASCII_099,
        'ConsistentTrees_ASCII_101': ConsistentTrees_ASCII_101,
        'stats_basic_bucket':   stats_basic_bucket,
        'stats_perc_bucket':    stats_perc_bucket,
        'ROCKSTAR_ASCII':       ROCKSTAR_ASCII,
        'ROCKSTAR_binary':      ROCKSTAR_binary,
        'ROCKSTAR_ASCII_list':  ROCKSTAR_ASCII_list,
        'Gadget_binary':        Gadget_binary 
        }

    func = choose.get(name)
    return func()


def load_ASCII_AHF():   
    path='/data/AHF_115Mpc_Wfirst/snapshot_110.z9.745.AHF_halos'

    #ID(1)	hostHalo(2)	numSubStruct(3)	Mvir(4)	npart(5)	Xc(6)	Yc(7)	Zc(8)	VXc(9)	VYc(10)	VZc(11)	Rvir(12)	Rmax(13)	r2(14)	
    #mbp_offset(15)	com_offset(16)	Vmax(17)	v_esc(18)	sigV(19)	lambda(20)	lambdaE(21)	Lx(22)	Ly(23)	Lz(24)	b(25)	c(26)	
    #Eax(27)	Eay(28)	Eaz(29)	Ebx(30)	Eby(31)	Ebz(32)	Ecx(33)	Ecy(34)	Ecz(35)	ovdens(36)	nbins(37)	fMhires(38)	Ekin(39)	Epot(40)	
    #SurfP(41)	Phi0(42)	cNFW(43)	

    data_pandas= pd.read_csv(path, skiprows=1, \
                      names=['haloid','hostid','numSubStruct','mhalo','n_particle','x_pos','y_pos','z_pos','x_vel','y_vel', 'z_vel', 'rvir', 'rmax', 'r2', 'mbp_offset',\
                             'com_offset','vmax', 'vesc', 'vdisp', 'spin_Bullock', 'spin_parameter', 'x_ang', 'y_ang', 'z_ang',\
                             'b_to_a', 'c_to_a','x_Ea', 'y_Ea', 'z_Ea', 'x_Eb', 'y_Eb', 'z_Eb', 'x_Ec', 'y_Ec', 'z_Ec',\
                             'ovdens', 'nbins', 'fMhries', 'Ekin', 'Epot', 'SurfP', 'Phi0', 'NFW_con'], skip_blank_lines=True, delim_whitespace=True)       

    data = ha_lib.df_to_sarray(data_pandas)
    
    return data

def load_ASCII_ConsitentTrees(path):   

    dt   = get_dtype_sarray('ConsistentTrees_ASCII_101')
   
    data = ha_lib.df_to_sarray(pd.read_csv(path, comment='#', names=[k[0] for k in dt], delim_whitespace=True, dtype=dt))                              
        
    return data

 
def load_binary_GADGET():
    import ReadGadget as rG
    path='/data/512_GadgetCats/snapshot_050'
    
    header, numFiles=rG.loadgadget_header(path)

    start_byte=256

    dt = [
        ('Pos'             , (np.float32, 3)),           
        ('Vel'             , (np.float32, 3)),
        ('haloid'          , (np.uint32,1))
        ]
 
    f = open(path,'rb') 

    f.seek(start_byte)
    type_name=np.fromfile(f, np.int32, 1)

    start_byte+=12
    n_files=1
    buffer_byte=4
    
    data=np.zeros((header['N'].values[0]*n_files,), dtype=dt)

    for item in dt:
        f.seek(start_byte)
        #print 'start_byte:', start_byte, 'item:', item[0], item[1][0], 'dtype dims:', item[1][1],

        data[item[0]] = np.fromfile(f, item[1], header['N'].values[0])
        #print 'min/max:', min(data[item[0]]), '/', max(data[item[0]]), 'dtype uses',

        if str(np.dtype(item[1][0])).find('32')!=-1:
            #print '--> use 4bytes'
            size=4
        elif str(np.dtype(item[1][0])).find('64')!=-1:
            #print '--> use 8bytes'
            size=8           
        else:
            #print '--> DEFAULT use 4bytes'            
            size=4
                    
        start_byte=start_byte+size*header['N'].values[0]*item[1][1]+2*buffer_byte
    
    return data

def load_binary_ROCKSTAR(path, 
                         snapid=0, 
                         filename_prefix='', 
                         filename_suffix='', 
                         n_files=1, 
                         read_ascii=False,
                         myprops='all',
                         include_desc_info=True):
    
    """Function loads binary files from ROCKSTAR into a structured array, optionally print halo catalog from ascii data
    
    input:
    ================
    
        path:                   string, filename of the binary files from ROCKSTAR (common: halo_xx.bin where xx are the snapshot id)
        
        **keywords
            snapid:             int,             snapshot identification number (default is 0)
            filename_prefix:    string,          prefix of the filename which is identical in all files to read
            filename_suffix:    string,          suffix of the filename which is identical in all files to read
            n_files:            integer,         number of files to read per snapshot
            read_ascii:         bool,            if 'True' the halo information of the ascii file is also read and stored as structed array 'data_ascii' (default is False)
            myprops:            list of strings, either 'all' then all properties available or only those stated in this list (default is 'all')
            include_desc_info:  bool,            if True the descendent index of the halos will be included in the output data structure (beware of the limitations of this property!
                                                 see ROCKSTAR documentation)
        
        
    output:
    ================

        data:               structured array of the halo information (number of properties to be returned can be customized)
        particle_IDs_dict:  dictonary which stores the particle IDs per haloid {haloid1: array([ID1, ID2, ... IDn]), haloid2: array([IDn+1, ... IDn+m]), ... }
        particle_info_dict: dictonary which stores number of halos and particle per snapshot {snapid: {'n_halos': n_halos, 'n_particles': n_particles}
    """
    #print '\nREADING ROCKSTAR BINARY HALO CATALOG\n####################################\n'

    #get header_size (where to start when reading binary), size of halo information for one halo, data type structure, data type of particle information
    #such asnumber of halos, numer of particles, box size, particle mass and particle type
    header_size, halo_struct_size, dt, dt_header_info, bytes_to_header_info = get_dtype_sarray('ROCKSTAR_binary')
    
    particle_ID_dict = {}
    count=0
    #Load binary files
    for nfile in range(0,n_files,1):
        #read the ascii file format as pandas object (optional)        
        if read_ascii==True:
            print('++++++++++++++++++\nReading ... ', path+filename_prefix+str(nfile)+'.ascii',)               
            dt_ascii=ha_lib.get_dtype_sarray('ROCKSTAR_ASCII')
            data_ascii = ha_lib.df_to_sarray(pd.read_csv(path+filename_prefix+str(nfile)+'.ascii', comment='#', names=[k[0] for k in dt_ascii], sep=' ', dtype=dt_ascii))
                                                         
            print('--> succesfully!\nHalo catalog:\n', data_ascii, '\n')                
            
        #print 'Reading ... ', path+filename_prefix+str(nfile)+filename_suffix, 
        f = open(path+filename_prefix+str(nfile)+filename_suffix,'rb')
       
        f.seek(bytes_to_header_info)
        particle_info = np.fromfile(f, dt_header_info, 1)[0]

        #print '--> succesfully!'
        if particle_info['n_halos']>0:       
            #print '\ntot num halos:\t\t', particle_info['n_halos'], '\ntot num particles:\t', particle_info['tot_n_particles'], '\nparticle mass:\t\t', format(particle_info['m_particles'], '0.8e'), '\n'

            f.seek(header_size)
            data_this_file = np.fromfile(f, dt, particle_info['n_halos'])
                   
            #print data_this_file, '\n'
            byte_loc=header_size+halo_struct_size*particle_info['n_halos']
            
            f.seek(byte_loc)
        
            #Read particle IDs as block of n_particles x 8 bytes
            #print 'Reading particle IDs at byte loc:', byte_loc, 'bytes',
            particle_IDs = np.fromfile(f, np.int64, particle_info['tot_n_particles'])
            #print '--> successfully!\n'
        
            #connect haloid with particles bound to that halo
            particle_ID_dict.update({haloid: {'particle_IDs': particle_IDs[sum(data_this_file['n_particles'][:i]) : sum(data_this_file['n_particles'][:i])+data_this_file['n_particles'][i]]} for (i, haloid, descID) in zip(range(0,data_this_file['haloid'].size,1),data_this_file['haloid'],data_this_file['descIndex'])})

            if count==0:
                #print 'here create data!'
                data = np.zeros((data_this_file.size,),dtype=np.dtype([(k, data_this_file.dtype[k]) for k in data_this_file.dtype.names]))
                data = data_this_file
                count+=1
            else:
                data = np.append(data, data_this_file, axis=0)   

        f.close()
        
    #print('SNAPID:', snapid, 'SUCESSFULLY READ!\n')        

    if myprops!='all':
        data=data[myprops]
        
    import numpy.lib.recfunctions as rcfuncs
    data          = rcfuncs.append_fields([data], ['snapid'] , [np.zeros(data.size,)], dtypes=['i4'], usemask=False)
    data['snapid']= snapid
    
    if include_desc_info==True:
        #get descIndex
        data.sort(order=['haloid'], axis=0)
        data['descIndex']= get_descIndex(path, snapid)      

    particle_info_dict={snapid: {'n_halos': data.size, 'n_particles': sum(data['n_particles'])}}

    return data, particle_ID_dict, particle_info_dict

def get_descIndex(path, snapid):
    
    dt   = get_dtype_sarray('ROCKSTAR_ASCII_list')
    data = ha_lib.df_to_sarray(pd.read_csv(path+'out_'+str(snapid)+'.list', comment='#', names=[k[0] for k in dt], sep=' ', dtype=dt))
                                             
    data.sort(order=['haloid'], axis=0)
    
    return data['descIndex']

def load_merger_tree(path):
    
    dt=get_dtype_sarray('merger_trees_ASCII')
    
    return ha_lib.df_to_sarray(pd.read_csv(path, comment='#', names=[k[0] for k in dt], sep='\t', dtype=dt))
    