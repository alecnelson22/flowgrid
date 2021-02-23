import numpy as np

print('pres_k1r1-h')
pres_out = np.load('output/pres_k1r1-h/pres_k1r1-h_2.npz')
print(pres_out['arr_0'].shape)
print(pres_out['arr_0'])

print('por-k1')
por_in = np.load('./output/por-k1/por-k1_0.npz')
print(por_in['arr_0'].shape)
print(por_in['arr_0'])

print('well input')
well_input = np.load('./output/well_inputs_k1r1-h/well_inputs_k1r1-h.npz', allow_pickle=True)
print(well_input['arr_0'].shape)
print(well_input['arr_0'])

#print('well output')
well_output = np.load('./output/well_outputs_k1r1-h/well_outputs_k1r1-h.npz', allow_pickle=True)
print(well_output['arr_0']) # it's a dictionary
