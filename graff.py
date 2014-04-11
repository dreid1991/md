

def graffme(plt, polyPts, atomPts):
	colors = ['red', 'green', 'blue', 'brown']
	for idx, poly in enumerate(polyPts):
		plt.fill(poly[0], poly[1], fill = False)
	plt.scatter(atomPts[0], atomPts[1])
	plt.show()



def displacement(plt, dispTuples):
	print 'HERE'

	times = [tup[0] for tup in dispTuples]
	disps = [tup[1] for tup in dispTuples]
	print dispTuples
	plt.plot(times, disps)
	plt.show()


def printRe(leRe):
	print 'ere'
	print leRe

def isNone(x):
	print 'func'
	return int(not bool(x))

