import numpy as np
from shapely.geometry import Point

# Create a class for gridding
class square_grid:
    
    def __init__(self, num_of_cells, scale_factor, centre_offset):
        
        self.number_of_cells = num_of_cells
        self.centre_offset = centre_offset
        
        # Create initial grid locations
        gridpoints = np.linspace(0, scale_factor*num_of_cells, num=num_of_cells)
        self.x_coords = gridpoints-centre_offset[0]
        self.y_coords = gridpoints-centre_offset[1]           
        
        
    def num_to_gridcoords(self, cell_number):
        return np.unravel_index(cell_number, (self.number_of_cells, self.number_of_cells), order='F')
        
    def gridcoords_to_num(self, x_grid, y_grid):
        return np.ravel_multi_index((x_grid, y_grid), (self.number_of_cells, self.number_of_cells), order='F')
    
    def num_to_geocoords(self, cell_number):
        (x_grid, y_grid) = self.num_to_gridcoords(cell_number)
        x_geo = self.x_coords[x_grid]
        y_geo = self.y_coords[y_grid]
        return (x_geo, y_geo)
    
    def geocoords_to_num(self, x_geo, y_geo):
        x_grid = np.argmin(abs(self.x_coords-x_geo))
        y_grid = np.argmin(abs(self.y_coords-y_geo))
        return self.gridcoords_to_num(x_grid, y_grid)
    
    def PointGeoms(self, indices):
        loclist = []
        for ind in indices:
            loclist.append(Point(self.num_to_geocoords(ind)))
        return loclist
    
#     def Borders(
    
    
# Define 