import numpy as np

f = open('Commands.txt','w')
omega = np.linspace(0, 1, 6)

for o in omega:
	print(o);
	cmdline = 'nohup '+'./simDER' + ' option.txt ' + '--' + ' friction ' +  '%.2f' % o + '\n'
	f.write(cmdline)

f.close()
