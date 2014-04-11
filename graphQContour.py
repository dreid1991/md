

def graphQContour(plt, xs, ys, zs, levels, xSize, ySize, maxLevel, maxVal, label, atomXs, atomYs):
	fig = plt.figure(figsize=(xSize, ySize))
	myPlt = fig.add_subplot(111)
	if not len(levels):
		levels = None
	contourStuff = myPlt.contourf(xs, ys, zs, levels = levels)
	axisLabelMax = 0.
	if maxLevel != -1:
		axisLevelMax = maxLevel
	else:
		axisLevelMax = maxVal
	ticks = []
	for i in range(int(axisLevelMax * 10)):
		ticks.append(i / 10.)
	cbar = fig.colorbar(contourStuff, ticks=ticks)
	cbar.ax.set_ylabel(label)
	if len(atomXs):
		plt.scatter(atomXs, atomYs, s=5)
	plt.show()
	#print 'xs '
	#print xs
	#print 'ys '
	#print len(ys)
	#print 'zs'
	#print len(zs)
	#print len(zs[0])
	#plt.show()
	#print "here!"
	#return 5

def addThings(a):
	for x in a:
		print x
