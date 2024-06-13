import sys

def readData(inname):
	data = {}
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		data[parts[0]] = parts[1:]
	f.close()
	return data

def readNet(inname):
	net = {}
	names = set([])
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		tg = parts[1]
		tf = parts[0]
		t  = net.get(tg,set([]))
		t.add(tf)
		net[tg] = t
		names.add(tf)
		names.add(tg)
	f.close()
	return (net,names)

def writeAdj(outname,names,net):
	nmap = {}
	for i in range(len(names)):
		nmap[names[i]] = i+1
	f = open(outname,'w')
	for tg in net:
		t = net[tg]
		for tf in t:
			f.write('%d\t%d\t1\n' % (nmap[tg],nmap[tf]))
	f.close()

def writeData(outname,names,data):
	f = open(outname,'w')
	for n in names:
		f.write('%s\t%s\n' % (n,'\t'.join(data[n])))
	f.close()

if __name__ == '__main__':
	data1 = readData(sys.argv[1])
	(net,names) = readNet(sys.argv[2])
	names = list(names)
	names.sort()
	writeAdj(sys.argv[3],names,net)
	writeData(sys.argv[4],names,data1)
