import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from utils.utils import  sumVect, cancers, genders, rows,rows_genes
from collections import Counter
def step17(algo):
	print ("step 17 start")
	input_file = algo[:-4] + '/data/'+algo[:-4] +'_patients.tsv'
	algo_name = algo[:-4]
	input_file_genes = algo[:-4] +'/data_genes/'+algo[:-4] + '_patients_genes.tsv'

	from textwrap import wrap
	try:
		os.mkdir(algo_name)
	except:
		pass
	try:
		os.mkdir(algo_name+"/cumulative histograms")
	except:
		pass
	try:
		os.mkdir(algo_name+"/data_genes")
	except:
		pass
	try:
		os.mkdir(algo_name+"/genes plots")
	except:
		pass
	try:
		os.mkdir(algo_name+"/genes plots/stage")
	except:
		pass
	header = "Tumor stage"+"\t"+"\t".join(rows)+"\n"
	header_genes = "Tumor stage"+"\t"+"\t".join(rows_genes)+"\n"
	stages1 = ['I']#'I or II NOS'
	stages2 = ['II']
	stages3 = ['III']
	stages4 = ['IV']
	colors_genes = {0:'lightcoral',1:'red',2:'maroon',3:'cyan',4:'c',5:'teal',6:'orange',7:'darkorange',8:'violet',9:'purple',10:'lightgreen'}
	total_events = {}
	male_events ={}
	female_events = {}
	male_events_genes ={}
	female_events_genes = {}
	total_events_genes = {}

	with open(input_file,'r') as inpt:
		inpt.readline().split("\t")
		
		for line in inpt:
			s = line.split("\t")
			if s[4] in stages1:
				key = '1'
			if s[4] in stages2:
				key = '2'
			if s[4] in stages3:
				 key = '3'
			if s[4] in stages4:
				key = '4'
			if key not in total_events.keys():
				total_events[key] = [0 for i in range(11)]
				total_events[key].append(set())
			total_events[key][-1].add(s[0]) #add patient
			for i in range(11):
				total_events[key][i] += int(s[5+i])
			if s[2] == "MALE":
				if key not in male_events.keys():
					male_events[key] = [0 for i in range(11)]
					male_events[key].append(set())
				male_events[key][-1].add(s[0]) #add patient
				for i in range(11):
					male_events[key][i] += int(s[5+i])
			###females
			if s[2] == "FEMALE":
				if key not in female_events.keys():
					female_events[key] = [0 for i in range(11)]
					female_events[key].append(set())
				female_events[key][-1].add(s[0]) #add patient
				for i in range(11):
					female_events[key][i] += int(s[5+i])

					
	with open(input_file_genes,'r') as inpt:
		inpt.readline().split("\t")
		
		for line in inpt:
			s = line.split("\t")
			if s[4] in stages1:
				key = '1'
			if s[4] in stages2:
				key = '2'
			if s[4] in stages3:
				 key = '3'
			if s[4] in stages4:
				key = '4'
			if key not in total_events_genes.keys():
				total_events_genes[key] = [[] for i in range(11)]
				
			for i in range(10):
				dd = s[5+i]
				if i == 9:
					dd = s[5+i][:-1]
				total_events_genes[key][i].extend(dd.split(','))
				total_events_genes[key][10].extend(dd.split(','))
			if s[2] == "MALE":
				if key not in male_events_genes.keys():
					male_events_genes[key] = [[] for i in range(11)]
					
				for i in range(10):
					dd = s[5+i]
					if i == 9:
						dd = s[5+i][:-1]
					male_events_genes[key][i].extend(dd.split(','))
					male_events_genes[key][10].extend(dd.split(','))
			###females
			if s[2] == "FEMALE":
				if key not in female_events_genes.keys():
					female_events_genes[key] = [[] for i in range(11)]
					

				for i in range(10):
					dd = s[5+i]
					if i == 9:
						dd = s[5+i][:-1]
					female_events_genes[key][i].extend(dd.split(','))
					female_events_genes[key][10].extend(dd.split(','))
					
	keys = ['1','2','3','4']
	keys_ = {'1':'I','2':'II','3':'III','4':'IV'}
	with open(algo_name+'/data_genes/'+algo_name+"_distribution_stages_genes.tsv",'w') as out:
		out.write(header_genes)
		genders_plot_genes = {}
		for c in keys:
			out.write(c+"\t")
			a = total_events_genes[c]#[:-1]
			a = [Counter(i).most_common() for i in a]
			line = []
			genders_plot_genes[c] = a
			for j in a:
				line.append(",".join([i[0]+":"+str(i[1])  for i in j if i[0] not in ["-",'']]))
			line = "\t".join(line)+"\n"
			out.write(line)
	x = keys
	data = genders_plot_genes

	for i in x:#gender/stage...
		fig, ax =  plt.subplots(4, 3, squeeze=False)
		
		for j in range(11):#тип мутации
			k = 0
			title = rows_genes[j]
			#fig, ax = plt.subplots()
			xx,yy = j//3,j%3
			if j==8:
				xx,yy = j//3+1,0
			elif j>8:
				xx,yy = j//3,j%3+1
			for driver in data[i][j]:#top 10 genes
				if k <10 and driver[0] not in ['','-']:
					ax[xx][yy].barh(driver[0], int(driver[1]),color = colors_genes[j],linewidth = 0.1)
					k+=1
			if k!=0:
				for ik in range(k,9):
					ax[xx][yy].barh(str(ik+100), 0,color = colors_genes[j],linewidth = 0.1)
					plt.setp(ax[xx][yy].get_yticklabels()[-1], visible=False)
			ax[xx][yy].set_xlabel('Number of events', fontdict = {'fontsize':18})
			ax[xx][yy].set_title(title,pad = 18, fontdict = {'fontsize':20})
			fig.set_figwidth(36)
			fig.set_figheight(36)
			fig.set_facecolor('floralwhite')
			ax[xx][yy].set_facecolor('seashell')
			ax[xx][yy].invert_yaxis()
			ax[xx][yy].tick_params(axis="y", labelsize=16)
			ax[xx][yy].tick_params(axis="x", labelsize=13)
			if k == 0:
				ax[xx][yy].set_yticks([])
				ax[xx][yy].set_yticks([], minor=True)
		fig.delaxes(ax[-2][-1])
		plt.savefig(algo_name+"/"+"/genes plots/stage/"+algo_name+"_"+"distribution_stage_"+keys_[i]+".pdf",dpi = 600)
	with open(algo_name+'/data_genes/'+algo_name+"_distribution_stages_males_genes.tsv",'w') as out:
		out.write(header_genes)
		genders_plot_genes = {}
		for c in keys:
			out.write(c+"\t")
			a = male_events_genes[c]#[:-1]
			a = [Counter(i).most_common() for i in a]
			line = []
			genders_plot_genes[c] = a
			for j in a:
				line.append(",".join([i[0]+":"+str(i[1])  for i in j if i[0] not in ["-",'']]))
			line = "\t".join(line)+"\n"
			out.write(line)
	x = keys
	data = genders_plot_genes

	for i in x:#gender/stage...
		fig, ax =  plt.subplots(4, 3, squeeze=False)
		
		for j in range(11):#тип мутации
			k = 0
			title = rows_genes[j]
			#fig, ax = plt.subplots()
			xx,yy = j//3,j%3
			if j==8:
				xx,yy = j//3+1,0
			elif j>8:
				xx,yy = j//3,j%3+1
			for driver in data[i][j]:#top 10 genes
				if k <10 and driver[0] not in ['','-']:
					ax[xx][yy].barh(driver[0], int(driver[1]),color = colors_genes[j],linewidth = 0.1)
					k+=1
			if k!=0:
				for ik in range(k,9):
					ax[xx][yy].barh(str(ik+100), 0,color = colors_genes[j],linewidth = 0.1)
					plt.setp(ax[xx][yy].get_yticklabels()[-1], visible=False)
			ax[xx][yy].set_xlabel('Number of events', fontdict = {'fontsize':18})
			ax[xx][yy].set_title(title,pad = 18, fontdict = {'fontsize':20})
			fig.set_figwidth(36)
			fig.set_figheight(36)
			fig.set_facecolor('floralwhite')
			ax[xx][yy].set_facecolor('seashell')
			ax[xx][yy].invert_yaxis()
			ax[xx][yy].tick_params(axis="y", labelsize=16)
			ax[xx][yy].tick_params(axis="x", labelsize=13)
			if k == 0:
				ax[xx][yy].set_yticks([])
				ax[xx][yy].set_yticks([], minor=True)
		fig.delaxes(ax[-2][-1])
		plt.savefig(algo_name+"/"+"/genes plots/stage/"+algo_name+"_"+"distribution_stage_males_"+keys_[i]+".pdf",dpi = 600)
	with open(algo_name+'/data_genes/'+algo_name+"_distribution_stages_females_genes.tsv",'w') as out:
		out.write(header_genes)
		genders_plot_genes = {}
		for c in keys:
			out.write(c+"\t")
			a = female_events_genes[c]#[:-1]
			a = [Counter(i).most_common() for i in a]
			genders_plot_genes[c] = a
			line = []
			for j in a:
				line.append(",".join([i[0]+":"+str(i[1])  for i in j if i[0] not in ["-",'']]))
			line = "\t".join(line)+"\n"
			out.write(line)
			
	x = keys
	data = genders_plot_genes

	for i in x:#gender/stage...
		fig, ax =  plt.subplots(4, 3, squeeze=False)
		
		for j in range(11):#тип мутации
			k = 0
			title = rows_genes[j]
			#fig, ax = plt.subplots()
			xx,yy = j//3,j%3
			if j==8:
				xx,yy = j//3+1,0
			elif j>8:
				xx,yy = j//3,j%3+1
			for driver in data[i][j]:#top 10 genes
				if k <10 and driver[0] not in ['','-']:
					ax[xx][yy].barh(driver[0], int(driver[1]),color = colors_genes[j],linewidth = 0.1)
					k+=1
			if k!=0:
				for ik in range(k,9):
					ax[xx][yy].barh(str(ik+100), 0,color = colors_genes[j],linewidth = 0.1)
					plt.setp(ax[xx][yy].get_yticklabels()[-1], visible=False)
			ax[xx][yy].set_xlabel('Number of events', fontdict = {'fontsize':18})
			ax[xx][yy].set_title(title,pad = 18, fontdict = {'fontsize':20})
			fig.set_figwidth(36)
			fig.set_figheight(36)
			fig.set_facecolor('floralwhite')
			ax[xx][yy].set_facecolor('seashell')
			ax[xx][yy].invert_yaxis()
			ax[xx][yy].tick_params(axis="y", labelsize=16)
			ax[xx][yy].tick_params(axis="x", labelsize=13)
			if k == 0:
				ax[xx][yy].set_yticks([])
				ax[xx][yy].set_yticks([], minor=True)
		fig.delaxes(ax[-2][-1])
		plt.savefig(algo_name+"/"+"/genes plots/stage/"+algo_name+"_"+"distribution_stage_females_"+keys_[i]+".pdf",dpi = 600)
			
			
					
	keys = ['1','2','3','4']
	with open(algo_name+'/data/'+algo_name+"_distribution_stages.tsv",'w') as out:
		out.write(header)
		stages_plot = {}
		for c in keys:
			out.write(c+"\t")
			num = len(total_events[c][-1])
			a = total_events[c][:-1]
			a = list(map(lambda x:x/num,a))
			stages_plot[c] = a
			line = "\t".join(list(map(lambda x:str(x),a)))+"\n"
			out.write(line)

	x = keys
	data = stages_plot
	fig, ax = plt.subplots()
	title = 'Driver event distribution by cancer stage'
	for i in x:
		ax.bar(i, data[i][0],color = 'lightcoral')
		bott = data[i][0]
		ax.bar(i, data[i][1], bottom = bott,color = 'red')
		bott+= data[i][1]
		ax.bar(i, data[i][2], bottom = bott,color = 'maroon')
		bott+= data[i][2]
		ax.bar(i, data[i][3], bottom = bott,color = 'cyan')
		bott+= data[i][3]
		ax.bar(i, data[i][4], bottom = bott,color = 'c')
		bott+= data[i][4]
		ax.bar(i, data[i][5], bottom = bott,color = 'teal')
		bott+= data[i][5]
		ax.bar(i, data[i][6], bottom = bott,color = 'orange')
		bott+= data[i][6]
		ax.bar(i, data[i][7], bottom = bott,color = 'darkorange')
		bott+= data[i][7]
		ax.bar(i, data[i][8], bottom = bott,color = 'violet')
		bott+= data[i][8]
		ax.bar(i, data[i][9], bottom = bott,color = 'purple')
	ax.set_ylabel('Average number driver of events', fontdict = {'fontsize':30})
	ax.set_xlabel('Cancer stage', fontdict = {'fontsize':30})
	ax.set_title(title, pad = 18,fontdict = {'fontsize':30})
	plt.xticks(fontsize = 19)
	plt.yticks(fontsize = 19)
	fig.set_figwidth(10)
	fig.set_figheight(16)
	fig.set_facecolor('floralwhite')
	ax.set_facecolor('seashell')
	plt.savefig(algo_name+"/"+"/cumulative histograms/"+algo_name+"_"+"distribution_stages"+".pdf",dpi = 600)

	with open(algo_name+'/data/'+algo_name+"_distribution_stages_males.tsv",'w') as out:
		out.write(header)
		males_plot = {}
		for c in keys:
			out.write(c+"\t")
			num = len(male_events[c][-1])
			a = male_events[c][:-1]
			a = list(map(lambda x:x/num,a))
			males_plot[c] = a
			line = "\t".join(list(map(lambda x:str(x),a)))+"\n"
			out.write(line)
	x = keys
	data = males_plot
	fig, ax = plt.subplots()
	title = 'Driver event distribution by cancer stage in males'
	for i in x:
		ax.bar(i, data[i][0],color = 'lightcoral')
		bott = data[i][0]
		ax.bar(i, data[i][1], bottom = bott,color = 'red')
		bott+= data[i][1]
		ax.bar(i, data[i][2], bottom = bott,color = 'maroon')
		bott+= data[i][2]
		ax.bar(i, data[i][3], bottom = bott,color = 'cyan')
		bott+= data[i][3]
		ax.bar(i, data[i][4], bottom = bott,color = 'c')
		bott+= data[i][4]
		ax.bar(i, data[i][5], bottom = bott,color = 'teal')
		bott+= data[i][5]
		ax.bar(i, data[i][6], bottom = bott,color = 'orange')
		bott+= data[i][6]
		ax.bar(i, data[i][7], bottom = bott,color = 'darkorange')
		bott+= data[i][7]
		ax.bar(i, data[i][8], bottom = bott,color = 'violet')
		bott+= data[i][8]
		ax.bar(i, data[i][9], bottom = bott,color = 'purple')

	ax.set_ylabel('Average number of driver events', fontdict = {'fontsize':30})
	ax.set_xlabel('Cancer stage', fontdict = {'fontsize':30})
	ax.set_title("\n".join(wrap(title,35)), pad = 18,fontdict = {'fontsize':30})
	plt.xticks(fontsize = 19)
	plt.yticks(fontsize = 19)
	fig.set_figwidth(10)
	fig.set_figheight(16)
	fig.set_facecolor('floralwhite')
	ax.set_facecolor('seashell')
	plt.savefig(algo_name+"/"+"/cumulative histograms/"+algo_name+"_"+"distribution_stages_males"+".pdf",dpi = 600)

	with open(algo_name+'/data/'+algo_name+"_distribution_stages_females.tsv",'w') as out:
		females_plot = {}
		out.write(header)
		for c in keys:
			out.write(c+"\t")
			num = len(female_events[c][-1])
			a = female_events[c][:-1]
			a = list(map(lambda x:x/num,a))
			females_plot[c] = a
			line = "\t".join(list(map(lambda x:str(x),a)))+"\n"
			out.write(line)
	x = keys
	data = females_plot
	fig, ax = plt.subplots()
	title = 'Driver event distribution by cancer stage in females'
	for i in x:
		ax.bar(i, data[i][0],color = 'lightcoral')
		bott = data[i][0]
		ax.bar(i, data[i][1], bottom = bott,color = 'red')
		bott+= data[i][1]
		ax.bar(i, data[i][2], bottom = bott,color = 'maroon')
		bott+= data[i][2]
		ax.bar(i, data[i][3], bottom = bott,color = 'cyan')
		bott+= data[i][3]
		ax.bar(i, data[i][4], bottom = bott,color = 'c')
		bott+= data[i][4]
		ax.bar(i, data[i][5], bottom = bott,color = 'teal')
		bott+= data[i][5]
		ax.bar(i, data[i][6], bottom = bott,color = 'orange')
		bott+= data[i][6]
		ax.bar(i, data[i][7], bottom = bott,color = 'darkorange')
		bott+= data[i][7]
		ax.bar(i, data[i][8], bottom = bott,color = 'violet')
		bott+= data[i][8]
		ax.bar(i, data[i][9], bottom = bott,color = 'purple')
	ax.set_ylabel('Average number of driver events', fontdict = {'fontsize':30})
	ax.set_xlabel('Cancer stage', fontdict = {'fontsize':30})
	ax.set_title("\n".join(wrap(title,35)),pad = 18, fontdict = {'fontsize':30})
	plt.xticks(fontsize = 19)
	plt.yticks(fontsize = 19)
	fig.set_figwidth(10)   
	fig.set_figheight(16)  
	fig.set_facecolor('floralwhite')
	ax.set_facecolor('seashell')
	plt.savefig(algo_name+"/"+"/cumulative histograms/"+algo_name+"_"+"distribution_stages_females"+".pdf",dpi = 600)
	print ("step 17 start")