# -*- coding: utf-8 -*-
# sh Python
# Single commodity primal dual algorithm ..
import sys, subprocess,datetime,time ,os, shutil,math,random


instance_name = str(sys.argv[1]); # this the argoment of python I give ls


arg_base = 'mpic++ -DNDEBUG -DIL_STD -c main_async.cpp -g' 
subprocess.call(arg_base,shell=True)
##os.sytem(arg_base1)
#raw_input("Checking 1")
arg_base = 'mpic++ -g -o main_async main_async.o -lilocplex -lconcert -lcplex -lm -lpthread' 
subprocess.call(arg_base,shell=True)
#raw_input("Checking 2")
				

if instance_name[0] == 'S':
	#Problem Info
	ORIGIN_DESTIN  =  [40]#14,40,80]
	#SCENARIOS      = [10,20,60,90]
	NODES	       = [16]#16,30]#16,
	
	for odkey in range(len(ORIGIN_DESTIN)):

		if ORIGIN_DESTIN[odkey] == 14:
			SCENARIOS      = [10]#,20]
		else:
			SCENARIOS      = [1000]#20,60,90]
			
		for nodekey in range(len(NODES)):
		
			# copy of instance file to data folder
			ins_name = 'SND'+str(ORIGIN_DESTIN[odkey])+'C'+str(NODES[nodekey])+'N'+'.dat'
			print ins_name

			fullPathins1 = "./instances-S/" + ins_name
			subprocess.call("cp" + " " + fullPathins1 + " " + "./",shell=True )
			os.rename("./" +ins_name,'instance.dat')
			subprocess.call("mv" + " " + './instance.dat' + " " + './Data',shell=True )		
			#raw_input("Checking 3")
			for sckey in range(len(SCENARIOS)):
				# copy of scenario file to data folder
				ins_name = 'scenario'+str(ORIGIN_DESTIN[odkey])+'C'+str(SCENARIOS[sckey])+'S'+'.dat'
				print ins_name

				fullPathins2 = "./instances-S/" + ins_name
				subprocess.call("cp" + " " + fullPathins2 + " " + "./",shell=True )
				os.rename("./" +ins_name,'scenario.dat')
				subprocess.call("mv" + " " + './scenario.dat' + " " + './Data',shell=True )		
				#raw_input("Checking 3")
				


				n_proc = 5

				dd = sys.argv[1]; # this the argoment of python I give S or R

				print 'mpirun --bind-to none -np ' + str(n_proc) \
				+' bash -c "ulimit -s unlimited &&'+' ./main_async '+dd[0]+' '+str(NODES[nodekey])+' '+str(ORIGIN_DESTIN[odkey])\
				+' '+str(SCENARIOS[sckey]) + ' ' + fullPathins1 + ' ' + fullPathins2+'"'
				
				rr= subprocess.call('mpirun --bind-to none -np ' + str(n_proc) \
				+' bash -c "ulimit -s unlimited &&'+' ./main_async '+dd[0]+' '+str(NODES[nodekey])+' '+str(ORIGIN_DESTIN[odkey])\
				+' '+str(SCENARIOS[sckey])+' '+'./Data/instance.dat'+' '+'./Data/scenario.dat'+' '+'"',shell=True)


					

				print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
				
			

#########----------------------------------------------------------------------------------------------------------------------------####

elif instance_name[0] == 'R':
	
	#Problem Info
	Ins_R_name         = [4]#4,6,7,8,9,10,11]#[4,5,6,7,8,9,10,11]#,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
	FIXCOST_CAP_RATIOS = [9]#1,3,5,7,9]   # there are 5 instances and I use only one (no: 5)
	CORRS    	   = [0.2]#[0,0.2,0.8]
	SCENARIOS 	   = [64]#16,32,64]
	t = -1
	for inskey in range(len(Ins_R_name)):
		if Ins_R_name[inskey] == 4:
			  n_nodes = 10
			  n_arcs  = 60
			  n_od    = 10
		elif Ins_R_name[inskey] == 5:
			  n_nodes = 10
			  n_arcs  = 60
			  n_od    = 25
		elif Ins_R_name[inskey] == 6:
			  n_nodes = 10
			  n_arcs  = 60
			  n_od    = 50
		elif Ins_R_name[inskey] == 7:
			  n_nodes = 10
			  n_arcs  = 82
			  n_od    = 10
		elif Ins_R_name[inskey] == 8:
			  n_nodes = 10
			  n_arcs  = 83
			  n_od    = 25
		elif Ins_R_name[inskey] == 9:
			  n_nodes = 10
			  n_arcs  = 83
			  n_od    = 50
		elif Ins_R_name[inskey] == 10:
			  n_nodes = 20
			  n_arcs  = 120
			  n_od    = 40
		elif Ins_R_name[inskey] == 11:
			  n_nodes = 20
			  n_arcs  = 120
			  n_od    = 100
		elif Ins_R_name[inskey] == 12:
			  n_nodes = 20
			  n_arcs  = 120
			  n_od    = 200
		elif Ins_R_name[inskey] == 13:
			  n_nodes = 20
			  n_arcs  = 220
			  n_od    = 40
		elif Ins_R_name[inskey] == 14:
			  n_nodes = 20
			  n_arcs  = 220
			  n_od    = 100
		elif Ins_R_name[inskey] == 15:
			  n_nodes = 20
			  n_arcs  = 220
			  n_od    = 200
		elif Ins_R_name[inskey] == 16:
			  n_nodes = 20
			  n_arcs  = 314
			  n_od    = 40
		elif Ins_R_name[inskey] == 17:
			  n_nodes = 20
			  n_arcs  = 318
			  n_od    = 100
		elif Ins_R_name[inskey] == 18:
			  n_nodes = 20
			  n_arcs  = 315
			  n_od    = 200
		elif Ins_R_name[inskey] == 19:
			  n_nodes = 25
			  n_arcs  = 100
			  n_od    = 10
		elif Ins_R_name[inskey] == 20:
			  n_nodes = 25
			  n_arcs  = 400
			  n_od    = 10
		elif Ins_R_name[inskey] == 21:
			  n_nodes = 25
			  n_arcs  = 100
			  n_od    = 30
		elif Ins_R_name[inskey] == 22:
			  n_nodes = 100
			  n_arcs  = 400
			  n_od    = 30
		elif Ins_R_name[inskey] == 23:
			  n_nodes = 20
			  n_arcs  = 230
			  n_od    = 40
		elif Ins_R_name[inskey] == 24:
			  n_nodes = 20
			  n_arcs  = 294
			  n_od    = 40
		elif Ins_R_name[inskey] == 25:
			  n_nodes = 30
			  n_arcs  = 519
			  n_od    = 100
		elif Ins_R_name[inskey] == 26:
			  n_nodes = 30
			  n_arcs  = 680
			  n_od    = 100
		elif Ins_R_name[inskey] == 27:
			  n_nodes = 20
			  n_arcs  = 230
			  n_od    = 200
		elif Ins_R_name[inskey] == 28:
			  n_nodes = 20
			  n_arcs  = 294
			  n_od    = 200
		elif Ins_R_name[inskey] == 29:
			  n_nodes = 30
			  n_arcs  = 520
			  n_od    = 300
		else:
			  n_nodes = 30
			  n_arcs  = 685
			  n_od    = 300

		for fcrkey in range(len(FIXCOST_CAP_RATIOS)):
			# copy of instance file to data folder
			ins_name = 'r0'+str(Ins_R_name[inskey]) + '.' + str(FIXCOST_CAP_RATIOS[fcrkey])+'.dow'
			if Ins_R_name[inskey] >= 10:
			      ins_name = 'r'+str(Ins_R_name[inskey]) + '.' + str(FIXCOST_CAP_RATIOS[fcrkey])+'.dow'
			
			print "Network file: " + ins_name
			fullPathins = "./instances-R/" + ins_name
			subprocess.call("cp" + " " + fullPathins + " " + "./",shell=True )
			os.rename("./" +ins_name,'instance.dat')
			subprocess.call("mv" + " " + './instance.dat' + " " + './Data',shell=True )
			
			for corrkey in range(len(CORRS)):
				for sckey in range(len(SCENARIOS)):
					# copy of scenario file to data folder
					ins_name = 'r0'+str(Ins_R_name[inskey]) +'-'+str(CORRS[corrkey])+'-'+str(SCENARIOS[sckey])

					if Ins_R_name[inskey] >= 10:
						ins_name = 'r'+str(Ins_R_name[inskey]) +'-'+str(CORRS[corrkey])+'-'+str(SCENARIOS[sckey])
					
					print "Scenario file: " + ins_name
					fullPathins1 = "./instances-R/scenarios/" + ins_name
					subprocess.call("cp" + " " + fullPathins1 + " " + "./",shell=True )
					os.rename("./" +ins_name,'scenario.dat')
					subprocess.call("mv" + " " + './scenario.dat' + " " + './Data',shell=True )		
					#raw_input("Checking 3")
					
					
					
					
					n_proc = 5
					#--------------------------------------------------------------------------------------------
					print 'mpirun --bind-to none -np ' + str(n_proc) \
					+' main_async '+instance_name[0]+' '+str(n_nodes)+' '+str(n_od)\
					+' '+str(SCENARIOS[sckey])+' '+ fullPathins +' '+ fullPathins1 
					
					
					dd = sys.argv[1]; # this the argoment of python I give S or R
					
					rr= subprocess.call('mpirun --bind-to none -np ' + str(n_proc) \
					+' main_async '+instance_name[0]+' '+str(n_nodes)+' '+str(n_od)\
					+' '+str(SCENARIOS[sckey])+' '+'./Data/instance.dat'+' '+'./Data/scenario.dat',shell=True)

					print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
					t =t+1
					
else:
	print "\nThe argment that you are entered should be S or R"
		
		
		
		
