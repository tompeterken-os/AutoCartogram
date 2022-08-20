import geopandas as gpd
import pandas as pd
import numpy as np
from scipy.stats import rankdata
import copy

def _homogenise_gdf_(
    geodataframe,
    id_header,
    name_header
):
    gdf = copy.deepcopy(geodataframe)
    
    # Add a dummy variable to allow cross-joins later
    gdf['__dummy__'] = gdf.apply(lambda row: 1, axis=1)
    
    # Find the centroids
    gdf['__centroid__'] = gdf.centroid
    
    return gdf[[id_header, name_header, gdf.geometry.name, '__dummy__', '__centroid__']]


class _Index_Id_Name_Transformer_:
    # Class to provide conversion methods between geography IDs, names, and indexes.
    # Indexes are not equivalent to the GeoDataFrame index column of the input file.
    
    def __init__(
        self,
        geo_list,
        name_list
    ):
        
        # Generate lookup dicts
        g_i = {}
        g_n = {}
        i_g = {}
        i_n = {}
        n_g = {}
        n_i = {}
        
        for i,g,n in zip(range(len(geo_list)), geo_list, name_list):
            g_i[g] = i
            g_n[g] = n
            i_g[i] = g
            i_n[i] = n
            n_g[n] = g
            n_i[n] = i
            
        self._geo_to_index_ = g_i
        self._geo_to_name_ = g_n
        self._index_to_geo_ = i_g
        self._index_to_name_ = i_n
        self._name_to_geo_ = n_g
        self._name_to_index_ = n_i
        
    # Define functions to use the dicts to translate
    def id_to_ix(self, code):
        return self._geo_to_index_[code]
    def id_to_name(self, code):
        return self._geo_to_name_[code]
    def ix_to_id(self, index):
        return self._index_to_geo_[index]
    def ix_to_name(self, index):
        return self._index_to_name_[index]
    def name_to_id(self, name):
        return self._name_to_geo_[name]
    def name_to_ix(self, name):
        return self._name_to_index_[name]
    
            
        
        
def _turn_df_into_matrix_(
    row_labels, col_labels, # Must have values of form ID
    values,
    matrix_size,
    transformer, # Must be of class _Index_Id_Name_Transformer_
    default_value=0
):
    # Takes a set of row-column-value lists and creates a 2D numpy array
    
    dtype = values.dtype
    matrix = np.ones((matrix_size, matrix_size), dtype=dtype)*default_value
    
    for row_lab, col_lab, value in zip(row_labels, col_labels, values):
        matrix[transformer.id_to_ix(row_lab), transformer.id_to_ix(col_lab)] = value
        
    return matrix


def _make_neighbours_array_(
    inf, # Must be class Input_File
    use_fractional_borders
):
    # Generates a 2D matrix of weights, where the weights are equal to the fraction of column geometry's
    # border which is shared with row geometry.

    # Add an extra geometry column to keep during crossjoining
    g_copy = copy.deepcopy(inf.geodata)
    g_copy['geocopy'] = g_copy.geometry

    # Find all geometries which touch each other
    neighbours = g_copy.sjoin(g_copy, how='inner', predicate='touches')
    neighbours = neighbours.loc[neighbours[inf.id_header+'_left']!=neighbours[inf.id_header+'_right']]

    if use_fractional_borders:
        # Find the length of the border as a fraction of the left geometry's total border
        neighbours['intersection'] = neighbours.apply(
            lambda row: row['geocopy_left'].intersection(row['geocopy_right']),
            axis=1
        )
        neighbours['weight'] = neighbours.apply(
            lambda row: row['intersection'].length/row['geocopy_left'].length,
            axis=1
        )

    else:
        neighbours['weight'] = neighbours.apply(lambda row: 1.0, axis=1)

    return _turn_df_into_matrix_(
        neighbours[inf.id_header+'_left'].values,
        neighbours[inf.id_header+'_right'].values,
        neighbours['weight'].values,
        inf.number_of_geographies,
        inf.lookups
    )


def _make_coastline_array_(
    inf
):
    # Generates a 1D vector of weights, where the weights are equal to the geometry's coast as a
    # fraction of its perimeter
    
    coastline_vector = 1-np.sum(inf.neighbours_array(use_fractional_borders=True), axis=1)
    coastline_vector[coastline_vector<0.01] = 0 # To allow for rounding errors
    return coastline_vector

def _make_distance_array_(
    inf
):
    # Generates a 2D matrix of the distances between each row/column geometry
    
    geom_name = inf.geodata.geometry.name
    
    # Do a cross-join
    crossjoin = inf.geodata.merge(
        inf.geodata,
        how='inner',
        on='__dummy__',
        suffixes=('_left', '_right')
    )
    
    # Find the separation between each pair
    crossjoin = crossjoin.loc[crossjoin[inf.id_header+'_left']!=crossjoin[inf.id_header+'_right']]
    crossjoin['distance'] = crossjoin.apply(
        lambda row: row[geom_name+'_left'].distance(row[geom_name+'_right']),
        axis=1
    )
    
    # Rank the distances
    crossjoin['distance_rank'] = crossjoin.groupby(inf.id_header+'_left')['distance'].rank()
    
    return _turn_df_into_matrix_(
        crossjoin[inf.id_header+'_left'].values,
        crossjoin[inf.id_header+'_right'].values,
        crossjoin['distance_rank'].values,
        inf.number_of_geographies,
        inf.lookups
    )
    
    