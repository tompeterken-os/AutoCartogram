import autocartogram_utils as utils
import geopandas as gpd
import copy


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
        
        self.geodata = utils._homogenise_gdf_(geodataframe, id_header, name_header)
        self.geo_list = self.geodata[id_header].values
        self.name_list = self.geodata[name_header].values
        self.number_of_geographies = len(geodataframe)
        
        self.lookups = utils._Index_Id_Name_Transformer_(self.geo_list, self.name_list)
        
        
    def neighbours_array(self, use_fractional_borders=True):        
        
        # We only need to run the whole process if it hasn't been done before, or if the use_fractional_borders has changed
        if '_neighbours_array_' in dir(self):
            need_to_analyse = False if self._neighbours_array_use_fractional_borders_ == use_fractional_borders else True
        else:
            need_to_analyse = True
        
        if need_to_analyse:
            self._neighbours_array_ = utils._make_neighbours_array_(self, use_fractional_borders)
            self._neighbours_array_use_fractional_borders_ = use_fractional_borders
            
        return self._neighbours_array_
    
    
    def coastline_array(self):        
        if '_coastline_array_' not in dir(self):
            self._coastline_array_ = utils._make_coastline_array_(self)        
        return self._coastline_array_
    
    
    def distance_array(self):        
        if '_distance_array_' not in dir(self):
            self._distance_array_ = utils._make_distance_array_(self)            
        return self._distance_array_
    
#     def 
            