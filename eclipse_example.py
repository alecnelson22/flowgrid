import flowgrid

grid = flowgrid.GRDECL()
fname_in = 'MILLION_CELL_MODEL/P50_SMALL_DOMAIN_GRID.GRDECL'
fname_out = 'MILLION_CELL_MODEL/SMARTSD_17.PRT'
fname_well = 'MILLION_CELL_MODEL/SMARTSD_17.RSM'
fname_perm = 'MILLION_CELL_MODEL/P50_SMALL_DOMAIN_PROP_PERMX.GRDECL'
fname_poro = 'MILLION_CELL_MODEL/P50_SMALL_DOMAIN_PROP_PORO.GRDECL'

# Set output directory name (optional, default is 'output')
grid.set_out_dir('MILLION_CELL_MODEL/output_compressed')

# Build grid geometry
grid.loadNodes(fname_in)

# Add input property to grid
poro = grid.read_prop(fname_poro, 'PORO')
perm = grid.read_prop(fname_perm, 'PERM')
grid.export_prop('PORO', poro)
grid.export_prop('PERM', perm)

# Get desired well info from .RSM file
well_output = grid.read_well_output(fname_well, ['WBHP', 'WWPR', 'WGPR', 'WGIR'])
grid.export_prop('WELLS', well_output)

# Pack property keywords to read from .PRT into list of lists
# This saves significant time for large grids as only one pass is required through .PRT file
# Inner lists should contain keywords that enable line denoting property section to be uniquely identified
# Order props as they appear in .PRT file
# grid.read_outputs(fname_out, [['PRESSURE', 'Max'], ['SGAS', 'Max']], toVTK=False, toNumpy=False, toCompressed=True)
grid.read_outputs(fname_out, [['PRESSURE', 'Max'], ['SGAS', 'Max'], ['SWAT', 'Max']], toVTK=True, toNumpy=True, toCompressed=False)
