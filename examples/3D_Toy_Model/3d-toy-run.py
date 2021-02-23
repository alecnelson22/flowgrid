import sys
import os
sys.path.insert(0, os.sep.join(['..', '..', 'regrid', 'flowgrid']))
import flowgrid
import time
import numpy as np

start_time = time.time()


# f_permi = './data/input/permi-k1.inc'
# f_por = './data/input/por-k1.inc'
# f_dat = './data/input/k1r1-h.dat'
# f_out = './data/output/k1r1-h.out'

# # Build CMG grid object
# grid = flowgrid.CMG()
# grid.CART(f_dat)

# grid.read_ext_prop(f_por,'por-k1')
# grid.read_ext_prop(f_permi,'permi-k1')

# grid.read_outputs(f_out, [
#     ['Gas Saturation', 'sat_k1r1-h'],
#     ['Pressure (psia)', 'pres_k1r1-h']
#     ])

# # grid.read_outputs(f_out, [
# #     ['Current Porosity', 'por_k1r1-h']
# #     ])

# # Build a Wells object from well input info
# wells = grid.get_wells(f_dat)

# # Read well output info, add to grid
# keys = ['Inst Surface Production Rates', 'Inst Surface Injection Rates']
# subkeys = [['Water', 'Gas'], ['Water', 'Gas']]
# well_output = wells.read_output(f_out, keys, subkeys)

# # Export well info as Numpy array
# # When reading back in with Numpy.load, set allow_pickle=True
# grid.export_wells(wells.wells, 'well_inputs')
# grid.export_wells(well_output, 'well_outputs')

# # Build and export VTK grids and output prop arrays for every timestep
# grid.export_grid('k1r1-h', toVTK=False, toNumpy=True)


input_dir = './data/input/'
output_dir = './data/output/'

perms = ['k1','k2','k3']
props = ['por','permi']
rates = {'k1':['r1','r1-5','r2','r2-5','r3','r3-5','r4','r5','r6'],
        'k2':['r1','r2','r3','r4','r5','r6','r7','r8','r9'],
        'k3':['r1','r2','r3','r4','r5','r6','r7','r8','r9']}
outputs = ['sat','pres']

# write input properties
for perm in perms:
    for rate in rates[perm]:

        # Build CMG grid object
        grid = flowgrid.CMG()

        f_name = perm + rate + '-h'
        f_dat = input_dir + f_name + '.dat'
        f_out = output_dir + f_name + '.out'

        print('Reading '+f_name)

        grid.CART(f_dat)

        for prop in props:
            f_inc = input_dir + prop + '-' + perm + '.inc'
            grid.read_ext_prop(f_inc, prop+'-'+perm)

        grid.read_outputs(f_out, [
            ['Gas Saturation', 'sat_'+f_name],
            ['Pressure (psia)', 'pres_'+f_name]
            ])

        # Build a Wells object from well input info
        wells = grid.get_wells(f_dat)

        # Optional: allows you to visualize wells in Paraview
        wells.build_cylinders(f_name, zscale=5)

        # Read well output info, add to grid
        keys = ['Inst Surface Production Rates', 'Inst Surface Injection Rates']
        subkeys = [['Water', 'Gas'], ['Water', 'Gas']]
        well_output = wells.read_output(f_out, keys, subkeys)

        # Export well info as Numpy array
        # When reading back in with Numpy.load, set allow_pickle=True
        grid.export_wells(wells.wells, 'well_inputs_'+f_name)
        grid.export_wells(well_output, 'well_outputs_'+f_name)

        # Build and export VTK grids and output prop arrays for every timestep
        grid.export_grid(f_name, toVTK=False, toNumpy=True)

print("--- %s minutes ---" % ((time.time() - start_time)/60))

