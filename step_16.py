from tqdm import tqdm
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from utils.utils import med, sumVect, cancers, genders,rows,rows_genes
from collections import Counter
def step16(algo):
	print ("step 16 start")
	input_file = algo[:-4] +'/data/'+algo[:-4] + '_patients.tsv'
	input_file_genes = algo[:-4] +'/data_genes/'+algo[:-4] + '_patients_genes.tsv'
	algo_name = algo[:-4]
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
		os.mkdir(algo_name+"/genes plots/gender")
	except:
		pass
	header = "Gender"+"\t"+"\t".join(rows)+"\n"
	header_genes = "Gender"+"\t"+"\t".join(rows_genes)+"\n"
	male_events ={}
	female_events = {}
	total_events = {}

	male_events_genes ={}
	female_events_genes = {}
	total_events_genes = {}
	
	with open(input_file,'r') as inpt:
		inpt.readline().split("\t")
		for line in inpt:
			s = line.split("\t")
			if s[2] not in total_events.keys():
				total_events[s[2]] = [0 for i in range(11)]
				total_events[s[2]].append(set())
			total_events[s[2]][-1].add(s[0]) #add patient
			for i in range(11):
				total_events[s[2]][i] += int(s[5+i])
#######################################################
	with open(input_file_genes,'r') as inpt:
		inpt.readline().split("\t")
		for line in inpt:
			s = line.split("\t")
			if s[2] not in total_events_genes.keys():
				total_events_genes[s[2]] = [[] for i in range(11)]
				
			#total_events_genes[s[2]][-1].add(s[0]) #add patient
			for i in range(10):
				dd = s[5+i]
				if i == 9:
					dd = s[5+i][:-1]
				total_events_genes[s[2]][i].extend(dd.split(','))
				total_events_genes[s[2]][10].extend(dd.split(','))
	with open(algo_name+'/data_genes/'+algo_name+"_distribution_gender_genes.tsv",'w') as out:
		out.write(header_genes)
		genders_plot_genes = {}
		for c in total_events_genes.keys():
			out.write(c+"\t")
			
			a = total_events_genes[c]#[:-1]
			a = [Counter(i).most_common() for i in a]
			genders_plot_genes[c] = a
			line = []
			for j in a:
				line.append(",".join([i[0]+":"+str(i[1])  for i in j if i[0] not in ["-",'']]))
			line = "\t".join(line)+"\n"
			out.write(line)
			
	with open(algo_name+'/data/'+algo_name+"_distribution_gender.tsv",'w') as out:
		out.write(header)
		genders_plot = {}
		for c in total_events.keys():
			out.write(c+"\t")
			num = len(total_events[c][-1])
			a = total_events[c][:-1]
			a = list(map(lambda x:x/num,a))
			genders_plot[c] = a
			line = "\t".join(list(map(lambda x:str(x),a)))+"\n"
			out.write(line)
	colors_genes = {0:'lightcoral',1:'red',2:'maroon',3:'cyan',4:'c',5:'teal',6:'orange',7:'darkorange',8:'violet',9:'purple',10:'lightgreen'}
			
	x = genders
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
		plt.savefig(algo_name+"/"+"/genes plots/gender/"+algo_name+"_"+"distribution_gender_"+i+".pdf",dpi = 600)
			
#############################3
	x = genders
	data = genders_plot
	fig, ax = plt.subplots()
	title = 'Driver event distribution by gender'
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
	ax.set_title(title,pad = 18, fontdict = {'fontsize':30})
	plt.yticks(fontsize = 19)
	plt.xticks(fontsize = 19)
	fig.set_figwidth(10)
	fig.set_figheight(16)
	fig.set_facecolor('floralwhite')
	ax.set_facecolor('seashell')
	plt.savefig(algo_name+"/"+"/cumulative histograms/"+algo_name+"_"+"distribution_gender"+".pdf",dpi = 600)
	print ("step 16 end")