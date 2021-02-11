import flowgrid

# Construct a single grid
grid = flowgrid.CMG()
fname_in = 'data/CASE30/r1h1_CASE30.dat'
fname_out = 'data/CASE30/r1h1_CASE30.out'
perm = 'data/CASE30/PERM.inp'
poro = 'data/CASE30/PORO.inp'

grid.CORNER(fname_in, ['ZCORN', 'DI', 'DJ'])

# # Use 'external' since these props are defined in standalone text files
# # If these props were defined in .dat, grid.readProperty(...) would be used
grid.read_ext_prop(poro, 'POR')
grid.read_ext_prop(perm, 'PERMI')
grid.read_ext_prop(perm, 'PERMJ')
grid.read_ext_prop(perm, 'PERMK', mult=.1)

# Formatted as [[GEM prop keyword, user defined prop title],...]
outputProps = [['Pressure(kpa)', 'Pressure']]

# This method should be performed before well operations so timesteps can be read
# Reads props for all timesteps and adds to grid#
grid.read_outputs(fname_out, outputProps)

# # Build a Wells object from well input info
wells = grid.get_wells('data/CASE30/wells.txt')

# Read well output info, add to grid
keys = ['Inst Surface Production Rates']
subkeys = [['Water', 'Gas']]
well_output = wells.read_output(fname_out, keys, subkeys)

# Export well info as Numpy array
# When reading back in with Numpy.load, set allow_pickle=True
grid.export_wells(well_output, 'well_outputs')

# Build and export VTK grids and output prop arrays for every timestep
grid.export_grid('r1h1_CASE30_', toVTK=True, toNumpy=True)
