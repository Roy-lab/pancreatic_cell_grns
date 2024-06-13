import sys

def readNames(inname):
	names = []
	f = open(inname,'r')
	for l in f:
		names.append(l.strip())
	f.close()
	return names

def readVal(inname):
	data = {}
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		data[parts[0]] = float(parts[1])
	f.close()
	return data

def readDeg(inname):
	data = {}
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		data[parts[0]] = int(parts[1])
	f.close()
	return data

def readNet(inname):
	net = {}
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		if parts[2] == 'NaN':
			parts[2] = '0.0'
		net[(parts[0],parts[1])] = parts[2]
	f.close()
	return net

def readAllVals(inname):
	allvals = []
	f = open(inname,'r')
	for l in f:
		vals = readVal(l.strip())
		allvals.append(vals)
	f.close()
	return allvals

def readAllNet(inname):
	allnets = []
	f = open(inname,'r')
	cnt = 0
	for l in f:
		cnt += 1
		if cnt == 1:
			conf = readNet(l.strip())
		else:
			net  =  readNet(l.strip())
			allnets.append(net)
	f.close()
	return (conf,allnets)

def writeAtt(outname,deg,mod,allvals,allnames):
	names = set([])
	for g in deg:
		names.add(g)
	for g in mod:
		names.add(g)
	for i in range(len(allnames)):
		f = open('%s_%s_att.txt' % (outname,allnames[i]),'w')
		f.write('Name\tType\tMod\tDeg\tExp\n')
		for g in names:
			if g in deg:
				t = 'TF'
			else:
				t = 'G'
			if g in mod:
				m = 'MOD'
			else:
				m = 'OUT'
			f.write('%s\t%s\t%s\t%d\t%f\n' % (g,t,m,deg.get(g,0),allvals[i][g]))
		f.close()

def writeNet(outname,conf,allnets,allnames,suff):
	f = open('%s_conf.txt' % (outname),'w')
	f.write('Source\tTarget\tConf\n')
	for (tf,tg) in conf:
		f.write('%s\t%s\t%s\n' % (tf,tg,conf[(tf,tg)]))
	f.close()
	for i in range(len(allnames)):
		f = open('%s_%s_%s.txt' % (outname,allnames[i],suff),'w')
		f.write('Name\tColor\n')
		for (tf,tg) in conf:
			f.write('%s (interacts with) %s\t%s\n' % (tf,tg,allnets[i][(tf,tg)]))
		f.close()

if __name__ == '__main__':
	(conf,allreg) = readAllNet(sys.argv[1])
	(conf,allcc)  = readAllNet(sys.argv[2])
	allvals = readAllVals(sys.argv[3])
	deg = readDeg(sys.argv[4])
	mod = readDeg(sys.argv[5])
	allnames = readNames(sys.argv[6])
	writeAtt(sys.argv[7],deg,mod,allvals,allnames)
	writeNet(sys.argv[7],conf,allreg,allnames,'reg')
	writeNet(sys.argv[7],conf,allcc,allnames,'cc')
