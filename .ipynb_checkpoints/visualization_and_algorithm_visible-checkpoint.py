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

#dipeptide df/tables
def make_dipep_tables(dipeptide):
    dipep = dipeptide
    
    subset_dipeptide = dpp.read_workup(dipep)
    df_dipeptide = dpp.split_description(subset_dipeptide)
    
    dpp.conformer_to_sql(dipep, df_dipeptide)
    
    return(True)


def make_conf_tables(dipeptide, dipeptide_conformer):
    dipep = dipeptide
    dipep_conformer = dipeptide_conformer
    
    conformer_df = dpp.read_coords(dipep, dipep_conformer)
    dpp.coord_to_sql(dipep_conformer, conformer_df)
    
    return(True)


def conf_visualization(dipeptide_conformer):
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

#let's go ahead and isolate these points 
def get_relavent_atoms(dipep_conformer):
    sql = f'''
    SELECT *
    FROM t{dipep_conformer}
    WHERE atom_idx IN (3, 4, 9, 10)
    '''

    globals()[f'{dipep_conformer}_isomer_df'] = pd.read_sql(sql, conn)
    
    return


def get_atom_coords(dipep_conformer, origin_point, end_point):
    conf_df = globals()[f'{dipep_conformer}_isomer_df']
    
    x1 = float(conf_df.loc[conf_df['atom_idx'] == origin_point]['x_coord'])
    y1 = float(conf_df.loc[conf_df['atom_idx'] == origin_point]['y_coord'])
    z1 = float(conf_df.loc[conf_df['atom_idx'] == origin_point]['z_coord'])
    
    x2 = float(conf_df.loc[conf_df['atom_idx'] == end_point]['x_coord'])
    y2 = float(conf_df.loc[conf_df['atom_idx'] == end_point]['y_coord'])
    z2 = float(conf_df.loc[conf_df['atom_idx'] == end_point]['z_coord'])
    
    return(x1, y1, z1, x2, y2, z2)


#function to get coordinate of interest from point of interest
def get_plane_vector(dipep_conformer, origin_point = '4', end_point = '9'):

    x1, y1, z1, x2, y2, z2 = get_atom_coords(dipep_conformer, origin_point, end_point)
    vector = np.array([x2-x1, y2-y1, z2-z1])
    
    return(vector)


def projection(p,a):
    lambda_val = np.dot(p,a)/np.dot(a,a)
    return p - lambda_val * a


def get_atom_vector(dipep_conformer, origin_point, end_point):
    
    x1, y1, z1, x2, y2, z2 = get_atom_coords(dipep_conformer, origin_point, end_point)
    
    plane_vector = get_plane_vector(dipep_conformer, '4', '9')
    
    atom1 = [x1, y1, z1]
    atom2 = [x2, y2, z2]
    
    atom1_p = projection(atom1, plane_vector)
    atom2_p = projection(atom2, plane_vector)
    
    vector = np.array(np.subtract(atom2_p, atom1_p))
    
    return(vector)


#function for calculating the angle using numpy
def angles(u, v): 
  #using the arccos function from numpy
  return arccos(u.dot(v)/(norm(u)*norm(v)))


def c_t_isomer(dipep_conformer):
    
    get_relavent_atoms(dipep_conformer)
    
    vec1 = get_atom_vector(dipep_conformer, '4', '3')
    vec2 = get_atom_vector(dipep_conformer, '9', '10')
    
    c = angles(vec1, vec2)

    #the function returns the angle in radians
    #converting the angle to degrees from radians
    angle= math.degrees(c)
    
    if (angle <= 90 or angle >= 270):
        isomer = 'c'
    else: 
        isomer = 't'
    
    return(angle, isomer)