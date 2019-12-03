#!/usr/bin/env python
import os
import random
import subprocess

numparts=0
settings=[[6, 6, 0.3, 0.1, 0.5, 1.0,'standard_smooth']]
#settings=[[6, 6, 0.3, 0.2, 0.2, 1,'high_amplitude_1']]
#thefiles=['saltdeanpot_30K.stl']
#thefiles=['ambercup_smooth','saltdeanpot_smooth','elephant_smooth']
#thefiles=['sphere']
seedLines = [239183,77364]
thefiles=['vaseFrank.stl']


filenum=1
for theseed in seedLines:
	for file in thefiles:
		print(file)
		# create 2 example of puzzle fragmented spheres
		for puzzleset in settings:
			generatedpiece=0;
			while(generatedpiece==0):
				# do between 7 and 14 parts
				#numparts = random.randint(7,14)
				#numparts = 5
				numparts = 6
				seedLine = theseed
				print('Puzzle '+str(puzzleset[6])+' with '+str(numparts)+' parts')
				breakSphereprogram='/Users/kre/workspace/research/breakSphere/source/breakSphere -numparts '+str(numparts)+' -seedLine '+str(seedLine)+' -m '+str(puzzleset[0])+' -niter '+str(puzzleset[1])+' -jitter '+str(puzzleset[2])+' -amplitude '+str(puzzleset[3])+' -decay '+str(puzzleset[4])+' -smoothing '+str(puzzleset[5])
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
				openscadcmd = '/Applications/MacPorts/OpenSCAD.app/Contents/MacOS/OpenSCAD -D \'file="geometry/'+file+'"\' -D num='+str(piece)+' -D geo='+str(filenum)+' -o '+resultfolder+'/puzzlepiece_'+str(piece)+'of'+str(numparts)+'.stl puzzlepiecesFrank.scad'
				print(openscadcmd)
				os.system(openscadcmd) 
		filenum=filenum+1
