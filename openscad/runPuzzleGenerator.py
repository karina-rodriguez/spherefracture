#!/usr/bin/env python
import os
import random
import subprocess

totalexamples=2
numparts=0
settings=[[6, 6, 0.3, 0.1, 0.85, 0.0,'standard'],[6, 6, 0.3, 0.1, 0.85, 1.0,'standard_smooth'],[7, 1, 0.1, 0.1, 0.5, 0.0,'low_niter'],[4, 4, 0.5, 0.1, 0.5, 1.0,'low_amplitude'],[4, 4, 0.5, 0.6, 0.5, 1.0,'high_amplitude']]
thefiles=['saltdeanpot_30K.stl']
#thefiles=['ambercup_smooth','saltdeanpot_smooth','elephant_smooth']

filenum=1
for file in thefiles:
	print(file)
	# create 2 example of puzzle fragmented spheres
	#for x in range(0, totalexamples):
	for puzzleset in settings:
		generatedpiece=0;
		while(generatedpiece==0):
			# do between 7 and 14 parts
			#numparts = random.randint(7,14)
			numparts = 10
			seedLine = 22
			seedOrientation = 7
			print('Puzzle '+str(puzzleset[6])+' with '+str(numparts)+' parts')
			breakSphereprogram='/Users/photoscan/workspace/breakSphere/source/breakSphere -numparts '+str(numparts)+' -seedLine '+str(seedLine)+' -seedOrientation '+str(seedOrientation)+' -m '+str(puzzleset[0])+' -niter '+str(puzzleset[1])+' -jitter '+str(puzzleset[2])+' -amplitude '+str(puzzleset[3])+' -decay '+str(puzzleset[4])+' -smoothing '+str(puzzleset[5])
			command_run = subprocess.call(breakSphereprogram,shell=True)
			if command_run == 0:
				print "Its worked!!"
				generatedpiece=1
			else:
				print "Try running again"


		#for each fragment piece, intersect with the 3D model
		resultfolder='result/'+file+'_puzzle_'+str(puzzleset[6])
		os.system('mkdir '+resultfolder)	
		for piece in range (0, numparts):
			print("I will break the piece: "+str(piece))
			os.system('cp fragmentsphere/fragment_'+str(piece)+'.stl '+resultfolder+'/')
			openscadcmd = 'openscad -D \'file="geometry/'+file+'"\' -D num='+str(piece)+' -D geo='+str(filenum)+' -o '+resultfolder+'/puzzlepiece_'+str(piece)+'of'+str(numparts)+'.stl puzzlepieces.scad'
			print(openscadcmd)
			os.system(openscadcmd) 
	filenum=filenum+1
