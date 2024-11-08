# plotDBthreads.py
# This here is a python script to load up the DBSCAN regressors, and then to plot them nicely and save them too! :)

# Here is the regressor file:
# subject='sub10448'
# folder='/Users/kevinaquino/projects/GSR_data/UCLA_data_niftis/fmriprep/'
# confound_file=folder+confounds

# regressors = 'sub10448_dbscan_liberal_regressors.csv'
# confounds  = 'sub-10448_task-rest_bold_confounds.tsv'
# folder_save = '/Users/kevinaquino/'

def main(raw_args=None):
	# Parse in inputs
	from argparse import ArgumentParser
	parser = ArgumentParser(epilog="plotDBthreads.py -- A function to generate a plot with all the regressors inculding FD. Kevin Aquino 2018 BMH")
	parser.add_argument("-s", dest="subject",
		help="The subject", metavar="subject_label")
	parser.add_argument("-l", dest="label",
		help="Label for the title and used for", metavar="label_")
	parser.add_argument("-d", dest="folder",
		help="folder for to save the carpet plot", metavar="saving_dir")
	parser.add_argument("-reg", dest="regressors",
		help="The regressors that DBSCAN outputs -- assumes the format that clusterCorrect.py outputs", metavar="regressors.csv")
	parser.add_argument("-cf", dest="confounds",
		help="The location of the confounds that is automatically generated by fmriprep", metavar="XXX_task-rest_bold_confounds.tsv")

	# Here we are parsing the arguments
	args = parser.parse_args(raw_args)

	# Parsing the arguments
	subject 	= args.subject
	label 		= args.label
	folder 		= args.folder
	regressors 	= args.regressors
	confounds 	= args.confounds
	# Now here have to plot the regressors:



	import csv	
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib import gridspec

	with open(regressors, 'r') as csvfile:
		spamreader = csv.reader(csvfile, delimiter=',')
		# row_count = sum(1 for row_dummy in spamreader)	
		total = []
		counter = 0;
		# counter = 
		for num,row in enumerate(spamreader):
			if(num>0):
				# print row			
				val2=np.array(row)			
				vals = val2.astype(float)			
				total = np.append(total,vals[1:len(vals)])
				counter=counter+1
				nframes=len(row)-1
				
	total2=np.reshape(total,(nframes,counter),order="F")


	confound_num = 0
	with open(confounds,'r') as tsvin:
		tsvin = csv.reader(tsvin, delimiter='\t')
		total_confounds = []	
		# counter = 
		for num,row in enumerate(tsvin):
			if(num==0):
				headings=row
			else:
				# print row			
				vals=np.array(row)			
				# import pdb; pdb.set_trace()
				# vals = val2.astype(float)			
				total_confounds = np.append(total_confounds,vals)		
				# confound_num = confound_num+1

	# import pdb; pdb.set_trace()
	total_confounds2=np.reshape(total_confounds,(len(headings),nframes),order="F")


	# Here start the plotting
	fig = plt.figure(figsize=(10, 3.5), dpi=300) 

	gs = gridspec.GridSpec(counter+2, 1) 

	# Setupcolor vector:
	cols = np.array(([1,0,0],[0,1,0],[0,0,1],[0.5,0.5,0],[0,0.5,0.5],[0,0.5,0.5]))

	# ====== Plot FD ========
	ax=plt.subplot(gs[0])
	ax.get_xaxis().set_visible(False)
	FD = total_confounds2[6,:]
	FD[0] = 0
	FD = FD.astype(float)
	plt.plot(FD)
	plt.autoscale(enable=True, axis='x', tight=True)
	ax.set_ylabel('FD')
	ax.yaxis.set_label_position("right")

	# Should also plot DVARS
	ax=plt.subplot(gs[1])
	ax.get_xaxis().set_visible(False)
	DVARS = total_confounds2[3,:]
	DVARS[0] = 1
	DVARS = DVARS.astype(float)
	plt.plot(DVARS)
	plt.autoscale(enable=True, axis='x', tight=True)
	ax.set_ylabel('DVARS')
	ax.yaxis.set_label_position("right")
	# =====-----------========

	# Plot the regressors #

	#for k in range(2,counter+2):
	#	ax=plt.subplot(gs[k])
	#	plt.plot(total2[:,k-2],c=cols[k-2,:])
	#	plt.autoscale(enable=True, axis='x', tight=True)
	#	ax.set_ylabel('REG:'+str(k-1))
	#	ax.yaxis.set_label_position("right")
		# plt.legend(loc=2)
	#	if(k<counter+2):
	#		ax.get_xaxis().set_visible(False)

	saveName=folder+subject+'_regressors.png'
	title=subject+' '+label
	plt.suptitle(title)
	plt.savefig(saveName)
	plt.clf()

	# Now also show some correlation plots just to make sure nothing is really correlated to FD?
	# for k in range(0,counter):
	# 	plt.scatter(FD,total2[:,k],c=cols[k,:],label='REG:'+str(k))	
	# plt.legend(loc=2)
	# saveName=folder_save+subject+'_regressors_cor.png'
	# plt.savefig(saveName)
	# plt.clf()


