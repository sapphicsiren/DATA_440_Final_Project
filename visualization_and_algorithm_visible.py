import data_preprocessing_visible as dpp 
from importlib import reload

reload(dpp)

import plotly.express as px
import pandas as pd

import numpy as np 
from numpy import arccos, array
from numpy.linalg import norm
import math

conn, c = dpp.Connect()


def make_conf_tables(dipeptide, dipeptide_conformer):
    '''
    Reads coordinate dataframe and converts it into
    a sql table
    '''
    
    dipep = dipeptide
    dipep_conformer = dipeptide_conformer
    
    conformer_df = dpp.read_coords(dipep, dipep_conformer)
    dpp.coord_to_sql(dipep_conformer, conformer_df)
    
    return(True)


def conf_visualization(dipeptide_conformer):
    '''
    Takes the coordinate information of the specified conformer
    and uses it to produce a 3-dimensional scatterplot
    of the atoms of interest. 
    
    Note: All hydrogens have been omitted to make
    visualization cleaner and easier to interpret 
    '''
    
    conn, c = dpp.Connect()
    
    sql = f'''
    SELECT *
    FROM t{dipeptide_conformer}
    WHERE element != 'H'
    '''

    conformer_noH = pd.read_sql(sql, conn)
    
    conn.commit()
    conn.close()
    
    fig = px.scatter_3d(conformer_noH, 
                    x='x_coord', y='y_coord', z='z_coord', 
                    color = 'element',
                    hover_data = ['atom_idx', 'element', 'x_coord', 'y_coord', 'z_coord'],
                    labels = {
                        'atom_idx': 'Atom Index',
                        'element': 'Element',
                        'x_coord': 'X-coordinate',
                        'y_coord': 'Y-coordinate',
                        'z_coord': 'Z-coordinate'})
    
    return(conformer_noH, fig)


def store_all_confs(dipeptide):
    '''
    Store the conformation coordinates (without H) and
    the corresponding figure as dynamically named variables
    that live within a larger dictionary corresponding
    to the dipeptide of interest 
    
    Naming:
    - ProXxx_df: the dipeptide dataframe
    - ProXxx_dict: the dipeptide dictionary
    '''
    
    dummy_dict = {}
    
    sql = f'''
    SELECT conformer_id, isomer, ring
    FROM t{dipeptide}
    '''
    
    #naming convention for dipeptide_df -> ProXxx_df
    globals()[f'{dipeptide}_df'] = pd.read_sql(sql, conn)
    
    for i in range(len(globals()[f'{dipeptide}_df'])):
        dipep_conformer = globals()[f'{dipeptide}_df'].iloc[i]['conformer_id'] 
        
        make_conf_tables(dipeptide, dipep_conformer)
        
        globals()[f'{dipep_conformer}_df'], globals()[f'{dipep_conformer}_fig']= conf_visualization(dipep_conformer)
        
        dummy_dict[dipep_conformer] = {dataframe: globals()[f'{dipep_conformer}_df'],
                                       figure: globals()[f'{dipep_conformer}_fig']}
        
    globals()[f'{dipeptide}_dict'] = dummy_dict
        
    return(globals()[f'{dipeptide}_dict'])

##########################################################################################

def get_relavent_atoms(dipep_conformer):
    '''
    Isolates the specific atoms of interest. 
    For ProXxx dipeptides, based on the way we process
    the data, we always know we are interested in the 
    atoms at atom_idx 3, 4, 9, and 10 
    
    Naming: 
    - {dipep_conformer}_isomer_df: dataframe of relevant atoms
    '''
    
    sql = f'''
    SELECT *
    FROM t{dipep_conformer}
    WHERE atom_idx IN (3, 4, 9, 10)
    '''

    globals()[f'{dipep_conformer}_isomer_df'] = pd.read_sql(sql, conn)
    
    return


def get_atom_coords(dipep_conformer, origin_point, end_point):
    '''
    Retreives the cartesian coordinates for the atoms of interest.
    Inputs are origin point and end point because we 
    are ultimately interested in comparing vectors that
    simulate the bond between two atoms (the one at origin and
    the one at the end)
    '''
    
    conf_df = globals()[f'{dipep_conformer}_isomer_df']
    
    x1 = float(conf_df.loc[conf_df['atom_idx'] == origin_point]['x_coord'])
    y1 = float(conf_df.loc[conf_df['atom_idx'] == origin_point]['y_coord'])
    z1 = float(conf_df.loc[conf_df['atom_idx'] == origin_point]['z_coord'])
    
    x2 = float(conf_df.loc[conf_df['atom_idx'] == end_point]['x_coord'])
    y2 = float(conf_df.loc[conf_df['atom_idx'] == end_point]['y_coord'])
    z2 = float(conf_df.loc[conf_df['atom_idx'] == end_point]['z_coord'])
    
    return(x1, y1, z1, x2, y2, z2)


def get_plane_vector(dipep_conformer, origin_point = '4', end_point = '9'):
    '''
    Default inputs set to work with ProXxx coordinates.
    Is used to get the vector that is orthogonal 
    to the plane we want to project our bonds onto. 
    
    Said bond is the "peptide bond" between C and N 
    '''
    
    x1, y1, z1, x2, y2, z2 = get_atom_coords(dipep_conformer, origin_point, end_point)
    vector = np.array([x2-x1, y2-y1, z2-z1])
    
    return(vector)


def projection(p,a):
    '''
    Vector projection function 
    '''
    
    lambda_val = np.dot(p,a)/np.dot(a,a)
    return p - lambda_val * a


def get_atom_vector(dipep_conformer, origin_point, end_point):
    '''
    Retrieves the vector projection of the bond of interest.
    The plane that is being projected onto is orthogonal to 
    the "peptide bond" vector. 
    '''
    
    x1, y1, z1, x2, y2, z2 = get_atom_coords(dipep_conformer, origin_point, end_point)
    
    plane_vector = get_plane_vector(dipep_conformer, '4', '9')
    
    atom1 = [x1, y1, z1]
    atom2 = [x2, y2, z2]
    
    atom1_p = projection(atom1, plane_vector)
    atom2_p = projection(atom2, plane_vector)
    
    vector = np.array(np.subtract(atom2_p, atom1_p))
    
    return(vector)


def angles(u, v): 
    '''
    Using the arccos function from numpy, get
    the angle between two vectors
    '''
    
    return arccos(u.dot(v)/(norm(u)*norm(v)))


def c_t_isomer(dipep_conformer):
    '''
    Classify conformer as cis or trans based on
    the angle between the two bonds of interest. 
    
    If the angle is equal to or between -90 degrees and 90 degrees
        it is a cis isomer
    
    If the angle is between 90 and 270 degrees
        the conformer is trans
    Else
        the conformer is cis
    '''
    
    get_relavent_atoms(dipep_conformer)
    
    vec1 = get_atom_vector(dipep_conformer, '4', '3')
    vec2 = get_atom_vector(dipep_conformer, '9', '10')
    
    c = angles(vec1, vec2)

    #the function returns the angle in radians
    #converting the angle to degrees from radians
    angle= math.degrees(c)
        
    if (90 < angle < 270):
        isomer = 't'
    else: 
        isomer = 'c'
    
    return(angle, isomer)

def validate_ct_labeling(dipeptide, dipeptide_conformer):
    '''
    Validate the cis/trans isomer labeling by comparing 
    the prediction to the manual label. 
    
    If the classification is correct, the function returns True
    and prints a statement with the isomer label and name of the conformer. 
    
    If the classification is wrong, the function returns False
    and prints a statement with the name of the conformer
    '''
    
    globals()[f'{dipeptide}_df'] = pd.read_csv('data/' + dipeptide + '_neutrals/' + dipeptide + '_workup.csv', index_col = 0) 
    dpp.conformer_to_sql(dipeptide, globals()[f'{dipeptide}_df'])

    sql = f'''
    SELECT isomer
    FROM t{dipeptide}
    WHERE conformer_id = '{dipeptide_conformer}'
    '''

    manual_isomer_label = pd.read_sql(sql, conn).iloc[0]['isomer']
    calc_isomer_label = c_t_isomer(dipeptide_conformer)[1]

    if calc_isomer_label == manual_isomer_label:
        print(f"Correct classification of {manual_isomer_label} isomer, {dipeptide_conformer}")
        return(True)
    else:
        print(f"Incorrect classification of {dipeptide_conformer}")
        return(False)