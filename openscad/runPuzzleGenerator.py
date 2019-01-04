#!/usr/bin/env python
import os
import random
import subprocess

totalexamples=5
numparts=0
#thefiles=['ambercuppoisson_100K.stl']
#thefiles=['tocutUSETHISFORPUZZLE_500K_poisson.stl']
thefiles=['elephantreconstructedsmoothclean125KWith7mmthickness.stl']

for file in thefiles:
	print(file)
	# create 2 example of puzzle fragmented spheres
	for x in range(0, totalexamples):
		generatedpiece=0;
		while(generatedpiece==0):
			# do between 7 and 14 parts
			#numparts = random.randint(7,14)
			numparts = random.randint(5,10)
			print('Puzzle '+str(x)+' with '+str(numparts)+' parts')
			breakSphereprogram='/Users/photoscan/workspace/breakSphere/source/breakSphere -numparts '+str(numparts)
			command_run = subprocess.call(breakSphereprogram,shell=True)
			if command_run == 0:
				print "Its worked!!"
				generatedpiece=1
			else:
				print "Try running again"


		#for each fragment piece, intersect with the 3D model
		resultfolder='result/'+file+'_puzzle'+str(x)
		os.system('mkdir '+resultfolder)	
		for piece in range (0, numparts):
			print("I will break the piece: "+str(piece))
			os.system('cp fragmentsphere/fragment_'+str(piece)+'.stl '+resultfolder+'/')
			openscadcmd = 'openscad -D \'file="geometry/'+file+'"\' -D num='+str(piece)+' -o '+resultfolder+'/puzzlepiece_'+str(piece)+'of'+str(numparts)+'.stl puzzlepieces.scad'
			print(openscadcmd)
			os.system(openscadcmd) 
