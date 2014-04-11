import math
import random
import sys

fileOut = open(sys.argv[1], 'w')
numAtoms = 1000
numSteps = 200
distStep = 1.
xlo = 0.
xhi = 30.
ylo = 0.
yhi = 30.

xspan = xhi-xlo
yspan = yhi-ylo
#x = 15.
#y = 15.
paths = [[] for i in range(numAtoms)]
for path in paths:
	x = random.random() * xspan + xlo
	y = random.random() * yspan + ylo
	for s in range(numSteps):

		relX = ((x-xlo) % xspan) / xspan
		relY = ((y-ylo) % yspan) / yspan
		path.append((relX, relY))
		theta = 2 * math.pi * random.random()
		x += math.cos(theta) * distStep
		y += math.sin(theta) * distStep

for i in range(numSteps):
	fileOut.write('ITEM: TIMESTEP\n')
	fileOut.write(str(i) + '\n')
	fileOut.write('ITEM: NUMER OF ATOMS\n')
	fileOut.write(str(numAtoms))
	fileOut.write('\n')
	fileOut.write('Item: Box Bounds pp ff pp\n')
	fileOut.write(str(xlo) + ' ' + str(xhi) + '\n')
	fileOut.write(str(ylo) + ' ' + str(yhi) + '\n')
	fileOut.write('-3 3\n')
	fileOut.write('ITEM: Atoms id type xs ys zs\n')
	for j, path in enumerate(paths):
		fileOut.write(str(j) + ' 1 ' + str(path[i][0]) + ' ' + str(path[i][1]) + ' 0.5\n')


