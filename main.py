import numpy as np

f = open('Commands.txt','w')
omega = np.linspace(1.0e-3, 2.0e-3, 11)

for o in omega:
	print(o);
	cmdline = 'nohup '+'./simDER' + ' option.txt ' + '--' + ' rodRadius ' +  '%.4f' % o + '\n'
	f.write(cmdline)

f.close()
