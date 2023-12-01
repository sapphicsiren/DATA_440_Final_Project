#Imports 
import pandas as pd
import re
import sys
import os 
import ast 
import sqlite3
import shutil

def Connect():
    #Create the database
    conn = sqlite3.connect('dipeptide_database')
    c = conn.cursor()

    #turn foreign keys on 
    c.execute("PRAGMA foreign_keys=ON;")
    
    return(conn, c)

#read the dipeptide workups
def read_workup(dipeptide):
    PATHWAY = 'data_private/' + dipeptide + '_workup.xlsx'
    workup = pd.read_excel(PATHWAY)
    
    #all the columns have different names, for simplicity manually subset before
    #uploading the file, and then change column names after 
    
    workup.columns = ['conformer', 'description']
    return(workup)


#description subsetting function 
#some of the descriptions for neutrals contain different amounts of info
#for ease of use, this function will go ahead and try to get only the first 
#two descriptors, cis/trans and ring in/out
def split_description(file):
    df_dat = []

    for dat in range(len(file)):
        dat_point = file.iloc[dat]
        dat_isomer = dat_point['description'].split(',')[0]
        dat_ring = dat_point['description'].split(',')[1]

        #we should also find a way to isolate the conformer id from the file name
        #as the log files are f.log but the conformers are c.log and a.log
        dat_conf = dat_point['conformer'][:-6]
    
        df_dat.append((dat_point['conformer'], dat_conf, dat_isomer, dat_ring))
    
    #populate the two dataframes 
    df = pd.DataFrame(df_dat)
    
    #What each column is:
    #conformer_file: the name of the log file that conformer is from
    #conformer_id: the actual identity of the conformer (same among all log file types)
    #isomer: cis or trans isomer
    #ring: proline ring has endo carbon or N (which position in/out (out of plane in ring)
    
    df.columns = ['conformer_file', 'conformer_id', 'isomer', 'ring']
    
    return df


#turn conformer_df into an sql table 
def conformer_to_sql(dipeptide, dipeptide_df):
 
    conn, c = Connect()
    
    tDipeptide = 't' + dipeptide
    tDipeptide_copy = tDipeptide + '_copy'
    
    dipeptide_df.to_sql(tDipeptide, conn, if_exists = 'replace', index = False)
    
    c.execute(f'''
    CREATE TABLE {tDipeptide_copy}( 
        conformer_id TEXT,
        conformer_file TEXT, 
        isomer TEXT,
        ring TEXT,
        PRIMARY KEY (conformer_id)
    )
    ''')

    c.execute(f'''
    INSERT INTO {tDipeptide_copy} (conformer_id, conformer_file, isomer, ring)
       SELECT conformer_id, conformer_file, isomer, ring FROM {tDipeptide}
    ''')

    c.execute(f'''
    DROP TABLE {tDipeptide}
    ''')

    c.execute(f'''
    ALTER TABLE {tDipeptide_copy} RENAME TO {tDipeptide}
    ''')
    
    conn.commit()
    conn.close()
    return(True)


def coord_to_sql(dipeptide_conformer, conformer_df):
    conn, c = Connect()
    
    tConformer = 't' + dipeptide_conformer
    tConformer_copy = tConformer + '_copy'
    
    conformer_df.to_sql(tConformer, conn, if_exists = 'replace', index = False)
    
    c.execute(f'''
    CREATE TABLE {tConformer_copy}( 
        atom_idx TEXT,
        element TEXT, 
        x_coord FLOAT,
        y_coord FLOAT,
        z_coord FLOAT,
        PRIMARY KEY (atom_idx)
    )
    ''')

    c.execute(f'''
    INSERT INTO {tConformer_copy} (atom_idx, element, x_coord, y_coord, z_coord)
       SELECT atom_idx, element, x_coord, y_coord, z_coord FROM {tConformer}
    ''')

    c.execute(f'''
    DROP TABLE {tConformer}
    ''')

    c.execute(f'''
    ALTER TABLE {tConformer_copy} RENAME TO {tConformer}
    ''')
    
    conn.commit()
    conn.close()
    return(True)


#subset the .log file to just return the coordinates 
#dipeptide: input must be ProXxx
def extract_coords(dipeptide, conformer_id):
    start = 0
    end = 0

    DATA_PATH = './data_private/'+ dipeptide + '_neutrals/'

    filename = conformer_id + "f.log"
    newfile = dipeptide + '_clean/' + filename[:-5] + '.final.opt.xyz'

    openold = open(DATA_PATH + filename, "r")

    os.makedirs(os.path.dirname(DATA_PATH + newfile), exist_ok=True)
    opennew = open(DATA_PATH + newfile, "w")

    rline = openold.readlines()
    for i in range (len(rline)):
        if "Standard orientation:" in rline[i]:
            start = i

    for m in range (start + 5, len(rline)):
        if "---" in rline[m]:
            end = m
            break

    ## Convert to Cartesian coordinates format
    ## convert atomic number to atomic symbol

    for line in rline[start+5 : end] :
        words = line.split()
        
        word0 = str(words[0])     #Index
        word1 = int(words[1])     #Identity of atom
        word3 = str(words[3])     #X-coords
        word4 = str(words[4])     #Y-coords
        word5 = str(words[5])     #Z-coords

        #only include relevant atoms to amino acids 
        if   word1 ==   1 : 
            word1 = "H"
        elif word1 ==   6 : 
            word1 = "C"
        elif word1 ==   7 : 
            word1 = "N"
        elif word1 ==   8 : 
            word1 = "O"
        elif word1 ==  16 : 
            word1 = "S"

        #print(word0, word1, word3)
        print([word0, word1, word3, word4, word5], file=opennew)

    openold.close()
    opennew.close()

    #let's also save a clean copy of the data to github
    # Specify the path of the destination directory you want to copy to
    destination_directory = 'DATA_440_Final_Project/data/' + dipeptide + '_neutrals/'

    if not os.path.exists(destination_directory):
        os.makedirs(destination_directory)
    
    # Use the shutil.copy() method to copy the file to the destination directory
    shutil.copy(DATA_PATH + newfile, destination_directory)
    
    print("#"*10 + " extract_coords() Done " + "#"*10)
    return(DATA_PATH + newfile)


#now let's read from the newly outputted file and turn it into a dataframe
def read_coords(dipeptide, conformer_id):
    dipep = dipeptide
    dipep_conformer = conformer_id
    
    filename = conformer_id + "f.log"
    dipep_file = 'data/' + dipeptide + '_neutrals/' + filename[:-5] + '.final.opt.xyz'

    f = open(dipep_file, "r")
    lines = f.readlines()

    coords_list = []
    for line in lines:
        coords_list.append(ast.literal_eval(line[:-1]))
    
    #generate a dataframe using the coordinates data 
    coords_df = pd.DataFrame(coords_list, columns = ['atom_idx', 'element', 'x_coord', 'y_coord', 'z_coord'])
    
    print("#"*10 + " read_coords() Done " + "#"*10)
    return(coords_df)