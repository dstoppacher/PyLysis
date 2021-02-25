"""
Helper function library and useful function to halo_analysis.py

"""
from __future__ import print_function
# Load packages

# own
import halo_analysis_load_cats as ha_lc

# system
import numpy as np

def find_cell_location(pos_vector, grid_res):
        
    cell_vector=[0,0,0]
    for i, pos in enumerate(pos_vector):
        pos=int(np.floor(pos_vector[i]/grid_res))
        if pos>=128:
            pos-=128
        cell_vector[i]=pos
    
    return cell_vector 

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
                    data_is_string=False):
    
    print('wirteINotFile:', filename) 
            
    if append_mytext!=True:
        f_handle = file(filename, 'w')  
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
