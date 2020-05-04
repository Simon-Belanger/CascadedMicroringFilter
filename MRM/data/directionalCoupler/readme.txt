# coupler = {'1500': None, '1548': None, '1600': None}

# f = open(os.getcwd() + '/fdtdSimul/coupler1500.txt','r')
# gap, t, k = [], [], []
# for line in f.read().split('\n')[0:-1]:
#     gap.append(complex(line.split('\t')[0].replace('i','j')))
#     k.append(complex(line.split('\t')[1].replace('i','j'))) 
#     t.append(complex(line.split('\t')[2].replace('i','j')))
# f.close()
# gap, k, t = np.asarray(gap), np.asarray(k), np.asarray(t)
# coupler['1500'] = {'gap': np.asarray(gap), 'k': np.asarray(k), 't': np.asarray(t)}

# f = open(os.getcwd() + '/fdtdSimul/coupler1548.txt','r')
# gap, t, k = [], [], []
# for line in f.read().split('\n')[0:-1]:
#     gap.append(complex(line.split('\t')[0].replace('i','j')))
#     k.append(complex(line.split('\t')[1].replace('i','j'))) 
#     t.append(complex(line.split('\t')[2].replace('i','j')))
# f.close()
# gap, k, t = np.asarray(gap), np.asarray(k), np.asarray(t)
# coupler['1548'] = {'gap': np.asarray(gap), 'k': np.asarray(k), 't': np.asarray(t)}

# f = open(os.getcwd() + '/fdtdSimul/coupler1600.txt','r')
# gap, t, k = [], [], []
# for line in f.read().split('\n')[0:-1]:
#     gap.append(complex(line.split('\t')[0].replace('i','j')))
#     k.append(complex(line.split('\t')[1].replace('i','j'))) 
#     t.append(complex(line.split('\t')[2].replace('i','j')))
# f.close()
# gap, k, t = np.asarray(gap), np.asarray(k), np.asarray(t)
# coupler['1600'] = {'gap': np.asarray(gap), 'k': np.asarray(k), 't': np.asarray(t)}
# saveData(coupler, os.getcwd()+ '/data/directionalCoupler/coupler1')