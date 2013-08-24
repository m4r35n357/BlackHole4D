#!/usr/bin/env python

from sys import argv
from math import sqrt
from visual import scene, sphere, curve, rate, ellipsoid
from json import loads

def main():
	if len(argv) < 3:
		raise Exception('>>> ERROR! Please supply a data file name and a parameter file name <<<')
	dataFile = open(argv[1], 'r')
	parameterFile = open(argv[2], 'r')
	# scene basics
	scene.center = (0,0,0)
	scene.width = scene.height = 1000
	scene.range = (20.0, 20.0, 20.0)
	# get parameters
	parameters = loads(parameterFile.read())
	a = parameters['a']
	horizon = parameters['M'] * (1.0 + sqrt(1.0 - a * a));
	# get data dimensions
	line = dataFile.readline()
	coordinates = loads(line)
	#  set up the balls
	colours = [ (1.0, 1.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.7, 0.7, 0.7), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (1.0, 1.0, 0.0) ]
	spheres = []
	ball = sphere(pos = (0.0, 0.0, 0.0), radius = horizon, color = colours[3], opacity = 1.0)
	spheres.append(ball)
	ball = sphere(pos = (coordinates['x'], coordinates['y'], coordinates['z']), radius = 0.2, color = colours[2])
	ball.trail = curve(color = ball.color)
	spheres.append(ball)
	myell = ellipsoid(pos = (0.0, 0.0, 0.0), length = 4.0, height = 4.0, width = 2.0 * horizon, color = colours[3], opacity = 0.2) 
	while line:
		rate(60)
		coordinates = loads(line)
		ball = spheres[1]
		error = coordinates['H']
		if error < -120.0:
			ball.color = colours[2]
		elif error < -90.0:
			ball.color = colours[5]
		elif error < -60.0:
			ball.color = colours[8]
		else:
			ball.color = colours[1]
		ball.trail.color = ball.color
		position = (coordinates['x'], coordinates['y'], coordinates['z'])
		ball.pos = position
		ball.trail.append(pos = position)
		line = dataFile.readline()

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
