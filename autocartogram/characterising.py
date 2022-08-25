import geopandas as gpd
import copy
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from shapely.geometry import MultiLineString


class Input_File:
    
    def __init__(
        self,
        geodataframe,
        id_header,
        name_header
    ):
        
        assert isinstance(geodataframe, gpd.GeoDataFrame), "Input must be a GeoDataFrame"
        assert (id_header != "__dummy__") & (name_header != "__dummy__"), "Sorry, '__dummy__' is a reserved column... please rename!"
        assert (id_header != "__centroid__") & (name_header != "__centroid__"), "Sorry, '__centroid__' is a reserved column... please rename!"
        
        self.id_header = id_header
        self.name_header = name_header
        
        self.geodata = _homogenise_gdf_(geodataframe, id_header, name_header)
        self.geo_list = self.geodata[id_header].values
        self.name_list = self.geodata[name_header].values
        self.number_of_geographies = len(geodataframe)
        
        self.lookups = _Index_Id_Name_Transformer_(self.geo_list, self.name_list)
                
        self.coastline = MultiLineString([x.exterior for x in geodataframe.dissolve()['geometry'].values[0].geoms])
        
        
    def neighbours_array(self, use_fractional_borders=True):        
        
        # We only need to run the whole process if it hasn't been done before, or if the use_fractional_borders has changed
        if '_neighbours_array_' in dir(self):
            need_to_analyse = False if self._neighbours_array_use_fractional_borders_ == use_fractional_borders else True
        else:
            need_to_analyse = True
        
        if need_to_analyse:
            self._neighbours_array_ = _make_neighbours_array_(self, use_fractional_borders)
            self._neighbours_array_use_fractional_borders_ = use_fractional_borders
            
        return self._neighbours_array_
    
    
    def coastline_array(self):        
        if '_coastline_array_' not in dir(self):
            self._coastline_array_ = _make_coastline_array_(self)        
        return self._coastline_array_
    
    
    def distance_array(self):        
        if '_distance_array_' not in dir(self):
            self._distance_array_ = _make_distance_array_(self)            
        return self._distance_array_
    
    
    def orientation_arrays(self):
        if '_orientation_arrays_' not in dir(self):
            self._orientation_arrays_ = _make_orientation_arrays_(self)
        return self._orientation_arrays_
    
    
    
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
    use_fractional_borders,
    only_non_coastal=True
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
        
        if only_non_coastal:
            neighbours['weight'] = neighbours.apply(
                lambda row: (
                    row['weight']/(
                        1-inf.coastline_array()[
                            inf.lookups.id_to_ix(
                                row[inf.id_header+'_left']
                            )
                        ]
                    )
                ),
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
    coastline_lengths = inf.geodata.intersection(inf.coastline).length.values
    return coastline_lengths/inf.geodata.length.values



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


def _make_orientation_arrays_(
    inf
):
    
    x_vals = inf.geodata['__centroid__'].apply(lambda c: c.xy[0][0]).values
    y_vals = inf.geodata['__centroid__'].apply(lambda c: c.xy[1][0]).values
    ew_ranks = rankdata(x_vals)-1
    ns_ranks = rankdata(y_vals)-1
    return (ns_ranks, ew_ranks)