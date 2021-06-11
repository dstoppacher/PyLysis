"""
Helper function library and useful function to halo_analysis.py

"""
from __future__ import print_function
# Load packages

# own
import halo_analysis_load_cats as ha_lc

# system
import numpy as np

def calcFastHisto(data,
                  filename,
                  col_name,
                  col_unit,
                  box_size,
                  binning='log',
                  nbins=15,
                  custom_min_max=False,
                  print_to_file=True,
                  comment=''):
          
    print('fast histogramm -->', col_name, '[', col_unit, ']',end=' ')

    if binning=='log':           
        data = np.log10(data)
        data = data[np.where(np.isfinite(data))[:][0]] 

    if custom_min_max==True:
        if col_name=='sfr':
            data_min = 1e-4
            data_max = 1000


        elif col_name=='ssfr':
            data_min = 1e-13
            data_max = 1e-8

        elif col_name=='mstar':
            data_min = 1e9
            data_max = 1e12                                          

        elif col_name=='mbh':
            data_min = 1e5
            data_max = 1e9
  
        elif col_name.find('mhalo')!=-1:
            data_min = 1e5
            data_max = 1e15
        elif col_name.find('age')!=-1:
            data_min = 0.1
            data_max = 10
        elif col_name.find('rhalfmass')!=-1:
            data_min = 1e-4
            data_max = 1                
        elif col_name=='vmax' or col_name=='vdisp':
            data_min = 10
            data_max = 10000
        elif col_name=='cgf':
            data_min = 1e-5
            data_max = 10                   
    else:
        data_min=min(data)
        data_max=max(data)
        
    data = data[np.where(data>=data_min)[:][0]]
    data = data[np.where(data<data_max)[:][0]]
          
        
#        data = data[np.where(data>=np.percentile(data,0.01))[:][0]]
#        data = data[np.where(data<np.percentile(data,99.9))[:][0]]
              
    print('data min/max:', min(data), '/', max(data),end='')
    
    binsize = (data_max-data_min)/nbins
    print('nbins:', nbins, 'binsize:', binsize)
    
    bins = np.linspace(data_min, data_max, nbins+1)
    print('bins:', bins)
   
    counts, edges = np.histogram(data, bins)

    histo=np.zeros((nbins, 7), dtype=np.float32)
   
    if binning =='log':
        histo[:,0]= (10**edges[1:]*10**edges[:-1])**0.5
    else:
        histo[:,0]= (edges[1:]+edges[:-1])/2
        
    print(edges)

    volume=box_size**3
    print(' binsize:', binsize, 'box_size:', box_size, 'volume:', volume)
    
    histo[:,1]=counts/binsize/volume
    histo[:,2]=binsize
    histo[:,3]=histo[:,2]
    histo[:,4]=histo[:,1]/counts**0.5
    histo[:,5]=histo[:,4]
    histo[:,6]=counts

    if col_unit.find('h-1')!=-1:
        y_unit='h^3 Mpc-3 dex-1'
    else:
        y_unit='Mpc-3 dex-1'
        
    header='\n(1) '+col_name+' ['+col_unit+'] (2) Phi ['+y_unit+'] (3) dx \t(4) -dx (5) -dy\t(6) +dy\t(7) N count'    
    
    if print_to_file==True:
        writeIntoFile(filename,
                       histo,
                       myheader='FAST HISTOGRAMM! cumulative: NO'+comment+header,
                       data_format="%0.8e",
                       mydelimiter='\t')

    return histo          


def create_data_array(dt_name, delimiter='\t'):
    """create a common structured array header information and string of formats (with delimiter specified) which can be used to print the data structure into file"""

    dt = ha_lc.get_dtype_sarray(dt_name)
                    
    data_array = np.zeros((1,),dtype=dt)

    header=''
    cols=[]
    format_string=''

    for i, item in enumerate(data_array.dtype.names):
        #print('item:', item, end='')
        unit, prop_format = find_unit(item)
        if item=='z':
            prop_format='%0.2f'
        #print(unit, prop_format)
        string='('+str(i+1)+') '+item+' '+unit+' '
        cols+=[item]
        format_string+=prop_format+delimiter
        header+=string
        
      
    return data_array, header, cols, format_string[0:-len(delimiter)], dt
 
def df_to_sarray(df):
    """
    Function converts a pandas data frame to a numpy stuctured array
    """    

    v = df.to_numpy()
    cols = df.columns
    #print(cols)
    
    types = [(k, df[k].dtype.type) for (i, k) in enumerate(cols)]
    #print(types)
    
    z = np.zeros((v[:,0].size,), dtype=types)
    
    for (i, k) in enumerate(z.dtype.names):
        z[k] = v[:, i]
    return z    

def find_cell_location(pos_vector, grid_res):
        
    cell_vector=[0,0,0]
    for i, pos in enumerate(pos_vector):
        pos=int(np.floor(pos_vector[i]/grid_res))
        if pos>=128:
            pos-=128
        cell_vector[i]=pos
    
    return cell_vector 

def find_members(data,
                 vector,
                 radius):
       
    return np.where(((data['pos_x']-vector[0])**2+(data['pos_y']-vector[1])**2+(data['pos_z']-vector[2])**2)**0.5 < radius)

def find_unit(prop):
    """return the common unit and format for each data type by searchning an comparing string"""

    props_unit = {'rvir':   ('[h-1Mpc]', '%0.5f'),
                  'mhalo':  ('[h-1Msun]','%0.5e'),
                  'pos1':   ('[h-1Mpc]','%0.5f'),
                  'pos2':   ('[h-1Mpc]', '%0.5f'),                    
                  'perc1':  ('[%]','%0.1f'),
                  'perc2':  ('[%]','%0.1f'),
                  'perc':   ('[%]','%0.1f'),                      
                  'ns':     ('[-]','%i'),
                  'shared': ('[-]','%i'),
                  'z':      ('[-]', '%0.2f'),                  
                  'log50':  ('log10(50th)', '%0.2f'),                 
                  'log10':  ('log10(10th)', '%0.2f'),
                  'log32':  ('log10(32th)', '%0.2f'),
                  'log68':  ('log10(68th)', '%0.2f'),
                  'log90':  ('log10(90th)', '%0.2f'),
                  '50':     ('50th', '%0.2f'),                 
                  '10':     ('10th', '%0.2f'),
                  '32':     ('32th', '%0.2f'),
                  '68':     ('68th', '%0.2f'),
                  '90':     ('90th', '%0.2f'),
                  'count':  ('[-]','%i')
                }    

    test = prop[::-1].find('_')

    if test==-1:
        try:
            unit=props_unit[prop]
        except:
            try:
                unit=props_unit[prop[:-1]]
            except: 
                unit=('[-]','%i')
    else:
        try:
            unit=props_unit[prop[:-1]]
        except:
            try:
                unit=props_unit[prop[::test]]
            except:
                unit=props_unit[prop[::-1][0:test][::-1]]                
                
    return unit[0], unit[1]

def get_props_for_stats_calc(cols,
                             input_is_dtype=False):
    """Function filters properties which allow startistical calculation from those which not. ID numbers or snapid are e.g. properties where
    no statistical information such as the median or percentiles can be caluclated!"""    
    
    try:
        cols=[k[0] for i,k in enumerate(cols)]
    except:
        pass

    props=[]
    for i in cols:
        if i.find('id')==-1 and i.find('Index')==-1 and i!='z1' and i!='z' and i!='z2' and i.find('n_count')==-1 and\
           i.find('50')==-1 and i.find('10')==-1 and i.find('32')==-1 and i.find('68')==-1 and i.find('90')==-1:
            props.extend([i])        
            
    return props

def read_simple_hdf52struct_array(path,
                             cols=[]):
    """Reads a set of columns dedicated by col from and hdf5-file to a structured array
    thereby it uses the dtype of each column to generate the array:
        
        input
        =========
            
        path:   the path to the hdf5-file
        col:    keyword, exact name of the column to be accessed in the hdf5-file,
                if not set, the column name from the hdf5-file is used and all columns
                in the hdf5-file are read to the structured array
                
        output
        =========
        
        returns a structured array with column names and dtypes and size as in the hdf5-file.
        The array dimension is taken from the hdf5-file. The shape of the structured array
        needs to be identical for all column which should be read in.
    
    """
    import h5py as hdf5
    
    f=hdf5.File(path, "r")
    print('\nSimle read-in hdf5 data to structured array! --> \npath:', path)
    
    if cols==[]:
        cols=f.keys()

    dt = np.dtype([(k, f[k].dtype) for k in cols])

    data = np.zeros((f[cols[0]].shape,), dtype=dt)
   
    #print(np.info(data))

    for name in cols:
        #print('property:', name,)
        for dim in range(len(f[cols[0]].shape)):
            #print('dim:', dim)
            try:
                data[name][dim,:]=f[name][dim]
            except:
               data[name]=f[name][:] 
    
    f.close()
        
    return data



def read_unshaped_csv(path, delimiter=' ', skip_n_rows=0):
        """function reads an unshaped text file meaning that each rows as a different amount of columns into a pandas data frame.
           code by P-S on Starkoverflow (question from 19/10/2018)
           https://stackoverflow.com/questions/27020216/import-csv-with-different-number-of-columns-per-row-using-pandas/27020945
           
           modified by me!
        """
        import pandas as pd
        
        # The max column count a line in the file could have
        largest_column_count = 0
        
        # Loop the data lines
        with open(path, 'r') as temp_f:
            # Read the lines
            lines = temp_f.readlines()
        
            for l in lines:
                # Count the column count for the current line
                column_count = len(l.split(delimiter)) + 1
        
                # Set the new most column count
                largest_column_count = column_count if largest_column_count < column_count else largest_column_count
        
        # Close file
        temp_f.close()
        
        # Generate column names (will be 0, 1, 2, ..., largest_column_count - 1)
        column_names = [i for i in range(0, largest_column_count)]
        
        # return csv
        return pd.read_csv(path, skiprows=skip_n_rows, delimiter=delimiter, names=column_names)

def read_unshaped_txt(path, nr_rows, nr_cols, start_from_row=0, delimiter=' ', comments='#', dtype=np.int32):
    """Reads unshaped (meaning that the amount of column can vary for each row) text files
    
    start_from_row=0 means that the reading starts from the 1st row (comments are excluded!)
    
    TODO! Needs to be fixed, last two column are not read in!
    """

    i = 0
    count = 0
    line_count = 0
    a=0
    data_array = np.zeros((nr_rows, nr_cols), dtype=dtype)
                      
    with open(path, 'r') as f: 
        
        for line in f:

            if line_count<nr_rows:
                data = np.fromstring(line, dtype=np.double, sep='\n')
                j = 0
                #print 'i start:', i
                if data.size==start_from_row and i!=0: 
                    data_array[i-1,:]=data_array[i,:]
                    i-=1
                    
                while (j<data.size):
                    
                    if data.size==start_from_row:
                        data_array[i,j] = data[j]
                        if j==data.size: count+=1
                        
                    else:
                        data_array[i,0:start_from_row] = data_array[i-1,0:start_from_row]
                        data_array[i-1,j+start_from_row] = data[j]
                       
                    j+=1
           
                i+=1
                line_count+=1 
                a+=1
 
    return data_array

def redshift_to_expfactor(redshift):
    
    """Converts the redshift to the expansion factor"""
    
    return 1.0/(redshift + 1.0) 


def expfactor_to_redshift(exp_factor):
    
    """Converts the expansion factor to the redshift""" 
   
    return 1.0 / exp_factor - 1

def scan_file_format_Cholla(f,path):   
    
    print('\n1) Attributes\n-------------------------')
        
    for k in range(0,len(f.attrs.keys()),1):
        print('\t', f.attrs.keys()[k], ':\t\t', f.attrs.values()[k][0])            
            
    print('\n+++++++++++++++++++++\nSCANNING FILE FORMAT of ...', path,'\n',\
              '\n2) Data\n-------------------------')                
    for name in f:
        print('name:\t', name,)            
        try:
            print('size:\t', f[name].size, f[name],)
            try:
                print('min/max:\t', format(min(f[name][:]), '0.2f'), '/', format(max(f[name][:]), '0.2f'))
            except:
                for k in [0,1,2]:
                    for g in [0,1,2]:
                        print('\tsize:', f[name][:,k][:,g].size, '\tdim: [:,'+str(k)+'][:,'+str(g)+']', format(min(f[name][:,k][:,g]), '0.2f'), '/', format(max(f[name][:,k][:,g]), '0.2f'))

        except:
            print('--> failed!\n')
   
def writeIntoFile(filename,
                    data2write,
                    mydelimiter='',
                    myheader='',
                    data_format='%.6e',
                    append_mytext=False,
                    data_is_string=False,
                    verbose=True):
    
    if verbose==False:
        print('wirteINotFile:', filename) 
            
    if append_mytext!=True:
        f_handle = open(filename, 'w')  
        if data_is_string==True:
            f_handle.write(data2write)
        else:
            np.savetxt(filename, data2write, fmt=data_format, header=myheader, delimiter=mydelimiter)
    else:
        f_handle = open(filename, 'a+')
        if data_is_string==True:
            f_handle.write(data2write)
        else:
            np.savetxt(f_handle, data2write, fmt=data_format, delimiter=mydelimiter)
    f_handle.close()
    
def create_DescScales_file(path_to_data,
                           snapid_start,
                           snapid_end):
    """Input file to run Consitent-Trees.
        Consistent-Trees needs a file which states the snapid and the scalefactor in an increasing manner
        e.g. >0 0.100
             >1 0.120
             >2 0.140             
        This function scan a certain folder called 'path_to_data' for 'out_*.list' files and stores the
        snapid and scale factors in '/path_to_data/DescScales.txt'.
        It assumes that the lowest snapid is 0 and by default scans 1000 
    """
    import pandas as pd
    
    for snap in range(snapid_start,snapid_end+1,1):

        path = path_to_data+'out_'+str(snap)+'.list'
        #print('snapid:', snap, end='')

        #read only the row where the scale factor is shown in the header of Rockstar output file!
        data = pd.read_csv(path, skiprows=1, nrows=1, dtype=str, index_col=False, header=None).at[0,0]
                                     
        scale_factor=data[data.find('= ')+1::]
        #print(', a:', scale_factor)
        if snap==0:
            writeIntoFile(
                        path_to_data+'DescScales.txt',
                        [str(snap)+' '+str(scale_factor)],
                        myheader='',
                        append_mytext=False,
                        data_is_string=False,
                        data_format='%s',
                        verbose=True)
        else:
            writeIntoFile(
                       path_to_data+'DescScales.txt',
                       str(snap)+' '+str(scale_factor)+'\n',
                       append_mytext=True,
                       data_is_string=True,
                       data_format='%s',
                       verbose=True)
            
    print('File successfully created -->', path_to_data+'DescScales.txt')
    
def create_snapidzred_file(path_to_data,
                           path_orginal_data_hydro,
                           snapid_start,
                           snapid_end):
    
    import h5py as hdf5
    
    for snap in range(snapid_start,snapid_end+1,1):
        print('snapid:', snap, end='')
        path = path_orginal_data_hydro+str(snap)+'_particles.h5.'+str(0)
        print('-->', path)
        print(path_to_data)
        f=hdf5.File(path, "r")
        
#        for k in range(0,len(f.attrs.keys()),1):
#            print('\t', f.attrs.keys()[k], ':\t\t', f.attrs.values()[k][0]) 
        
        if snap==snapid_start:
            writeIntoFile(
                        path_to_data+'snapidzred.txt',
                        [str(snap)+' '+str(f.attrs['Current_a'][0])+' '+str(f.attrs['Current_z'][0])+' '+str(f.attrs['t'][0])],
                        myheader=path_to_data+'\n(1) snapshot number (2) scale factor a (2) redshift z (4) t',
                        append_mytext=False,
                        data_is_string=False,
                        data_format='%s',
                        verbose=True)
        else:
            writeIntoFile(
                       path_to_data+'snapidzred.txt',
                       str(snap)+' '+str(f.attrs['Current_a'][0])+' '+str(f.attrs['Current_z'][0])+' '+str(f.attrs['t'][0])+'\n',
                       append_mytext=True,
                       data_is_string=True,
                       data_format='%s',
                       verbose=True)    
    

    print('File successfully created -->', path_to_data+'snapidzred.txt')
    
