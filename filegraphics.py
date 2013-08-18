#!/usr/bin/env python

from sys import argv
from visual import scene, sphere, curve, rate
from json import loads

def main():
	if len(argv) < 2:
		raise Exception('>>> ERROR! Please supply a data file name <<<')
	dataFile = open(argv[1], 'r')
	# scene basics
	scene.center = (0,0,0)
	scene.width = scene.height = 1000
	scene.range = (20.0, 20.0, 20.0)
	# get data dimensions
	line = dataFile.readline()
	coordinates = loads(line)
	#  set up the balls
	colours = [ (1.0, 1.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.7, 0.7, 0.7), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (0.3, 0.3, 0.3) ]
	spheres = []
	ball = sphere(pos = (0.0, 0.0, 0.0), radius = 2.0, color = colours[3])
	spheres.append(ball)
	ball = sphere(pos = (coordinates['x'], coordinates['y'], coordinates['z']), radius = 0.1, color = colours[2])
	ball.trail = curve(color = ball.color)
	spheres.append(ball)
	while line:
		rate(60)
		coordinates = loads(line)
		ball = spheres[1]
		position = (coordinates['x'], coordinates['y'], coordinates['z'])
		ball.pos = position
		ball.trail.append(pos = position)
		line = dataFile.readline()

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
