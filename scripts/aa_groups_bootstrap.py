###IMPORT
import sys
import numpy as np
import random
###FILES
#Add human reference proteome fasta from Uniprot to PROTEOME
PROTEOME=""
BOUNDARY="cc_prediction_boundaries.csv"
MUTATION="humsavar_mutations.csv"
CLUSTER="protein_clusters_cdhit.csv"
###FIXED_VAR
DATASET="human"#DEFINE: Deepcoil, Marcoil, Ncoils, Paircoil or human!
try:
	DATASET=sys.argv[1]
	if DATASET=="Deepcoil" or DATASET=="Marcoil" or DATASET=="Ncoils" or DATASET=="Paircoil" or DATASET=="human":
		pass
	else:
		print("DEFINE: Deepcoil, Marcoil, Ncoils, Parcoil or human")
		sys.exit()
except IndexError:
	print("DEFINE: Deepcoil, Marcoil, Ncoils, Paircoil or human")
	sys.exit()
groups={}
groups['H']='positive'
groups['K']='positive'
groups['R']='positive'
groups['D']='negative'
groups['E']='negative'
groups['A']='hydrophobic'
groups['I']='hydrophobic'
groups['L']='hydrophobic'
groups['V']='hydrophobic'
groups['M']='hydrophobic'
groups['C']='other'
groups['F']='other'
groups['G']='other'
groups['N']='other'
groups['P']='other'
groups['Q']='other'
groups['S']='other'
groups['T']='other'
groups['W']='other'
groups['Y']='other'
###READ FILES
inf=open(PROTEOME)
s=""
ID=""
seqs={}
while 1:
	line=inf.readline()
	if line=="":
		seqs[ID]=s
		break
	if line[0]==">":
		if s!="":
			seqs[ID]=s
		s=""
		ID=line.split("|")[1]
	else:
		s=s+line.strip() 
inf.close()


inf=open(BOUNDARY)
line=inf.readline()#header
boundary={}
s=0
while 1:
	line=inf.readline()
	if line=="":
		break
	line=line.split(";")
	AC=line[0]
	method=line[1]
	start=int(line[2])
	end=int(line[3].strip())
	try:
		boundary[AC][method].append([start,end])
	except KeyError:
		try:
			boundary[AC][method]=[[start,end]]
		except KeyError:
			boundary[AC]={method:[[start,end]]}
inf.close()

inf=open(MUTATION)
line=inf.readline()#header
mutation={}
while 1:
	line=inf.readline()
	if line=="":
		break
	line=line.split(";")
	AC=line[0]
	position=int(line[1])
	wt=line[2]
	mt=line[3]
	phenotype=line[4].strip()
	try:
		mutation[AC][position][phenotype].append([wt,mt])
	except KeyError:
		try:
			mutation[AC][position][phenotype]=[[wt,mt]]
		except KeyError:
			try:
				mutation[AC][position]={phenotype:[[wt,mt]]}
			except KeyError:
				mutation[AC]={position:{phenotype:[[wt,mt]]}}
inf.close()
s=0
s2=0
inf=open(CLUSTER)
line=inf.readline()#header
valid={}
while 1:
	line=inf.readline()
	if line=="":
		break
	line=line.split(";")
	represent=line[1]
	valid[represent]="nr"

inf.close()

###FUN
def cc_group(method,dataset,phenotype_,boundary_,mutation_,valid_,seqs_):
	mx={}
	for key in groups:
		mx[groups[key]]={}
		for key2 in groups:
			mx[groups[key]][groups[key2]]=0
	if method!="human":
		for AC in seqs_:
			if random.random()<0.8:
				accept=True
				minimx={}	
				for key in groups:
					minimx[groups[key]]={}
					for key2 in groups:
						minimx[groups[key]][groups[key2]]=0			
				try:
					if valid_[AC]==dataset:	
						for i in range (0, len(boundary_[AC][method])):			
							if len(seqs_[AC])>boundary_[AC][method][i][1]-1:
								for position in mutation_[AC]:
									if boundary_[AC][method][i][0]<=position<=boundary_[AC][method][i][1]:		
										for j in range (0, len(mutation_[AC][position][phenotype_])):	
											if seqs_[AC][position-1]==mutation_[AC][position][phenotype_][j][0]:
												minimx[groups[mutation_[AC][position][phenotype_][j][0]]][groups[mutation_[AC][position][phenotype_][j][1]]]+=1
											else:
												accept=False
							else:
								accept=False
				except KeyError:
					pass
				if accept==True:
					for key in mx:
						for key2 in mx[key]:
							mx[key][key2]+=minimx[key][key2]	
		return mx
		
	else:
		for AC in seqs_:
			if random.random()<0.8:
				accept=True
				minimx={}	
				for key in groups:
					minimx[groups[key]]={}
					for key2 in groups:
						minimx[groups[key]][groups[key2]]=0			
				try:
					if valid_[AC]==dataset:	
				
						for position in mutation_[AC]:												
							if position>len(seqs_[AC])-1:
								accept=False		
							else:	
								try:
									for j in range (0, len(mutation_[AC][position][phenotype_])):	
										
										if seqs_[AC][position-1]==mutation_[AC][position][phenotype_][j][0]:
											minimx[groups[mutation_[AC][position][phenotype_][j][0]]][groups[mutation_[AC][position][phenotype_][j][1]]]+=1
										else:
											accept=False
								except KeyError:
									pass
					else:
						accept=False
				except KeyError:
					pass
				if accept==True:
					for key in mx:
						for key2 in mx[key]:
							mx[key][key2]+=minimx[key][key2]	
		return mx
###MAIN
print(DATASET+" DM;hydrophobic;positive;negative;other;"+DATASET+" DM;hydrophobic;positive;negative;other")
changemx={}
for key in groups:
	changemx[groups[key]]={}
	for key2 in groups:
		changemx[groups[key]][groups[key2]]=[]
for i in range (0,2):
	changemx_=cc_group(DATASET,"nr","Disease",boundary,mutation,valid,seqs)
	for key in changemx_:		
		for key2 in changemx_[key]:
			changemx[key][key2].append(changemx_[key][key2])
print("hydrophobic;"+str(np.mean(changemx['hydrophobic']['hydrophobic']))+";"+str(np.mean(changemx['hydrophobic']['positive']))+";"+str(np.mean(changemx['hydrophobic']['negative']))+";"+str(np.mean(changemx['hydrophobic']['other']))+";hydrophobic;"+str(np.std(changemx['hydrophobic']['hydrophobic']))+";"+str(np.std(changemx['hydrophobic']['positive']))+";"+str(np.std(changemx['hydrophobic']['negative']))+";"+str(np.std(changemx['hydrophobic']['other'])))
print("positive;"+str(np.mean(changemx['positive']['hydrophobic']))+";"+str(np.mean(changemx['positive']['positive']))+";"+str(np.mean(changemx['positive']['negative']))+";"+str(np.mean(changemx['positive']['other']))+";positive;"+str(np.std(changemx['positive']['hydrophobic']))+";"+str(np.std(changemx['positive']['positive']))+";"+str(np.std(changemx['positive']['negative']))+";"+str(np.std(changemx['positive']['other'])))
print("negative;"+str(np.mean(changemx['negative']['hydrophobic']))+";"+str(np.mean(changemx['negative']['positive']))+";"+str(np.mean(changemx['negative']['negative']))+";"+str(np.mean(changemx['negative']['other']))+";negative;"+str(np.std(changemx['negative']['hydrophobic']))+";"+str(np.std(changemx['negative']['positive']))+";"+str(np.std(changemx['negative']['negative']))+";"+str(np.std(changemx['negative']['other'])))
print("other;"+str(np.mean(changemx['other']['hydrophobic']))+";"+str(np.mean(changemx['other']['positive']))+";"+str(np.mean(changemx['other']['negative']))+";"+str(np.mean(changemx['other']['other']))+";other;"+str(np.std(changemx['other']['hydrophobic']))+";"+str(np.std(changemx['other']['positive']))+";"+str(np.std(changemx['other']['negative']))+";"+str(np.std(changemx['other']['other'])))
print(DATASET+" PM;hydrophobic;positive;negative;other;"+DATASET+" PM;hydrophobic;positive;negative;other")
changemx={}
for key in groups:
	changemx[groups[key]]={}
	for key2 in groups:
		changemx[groups[key]][groups[key2]]=[]
for i in range (0,2):
	changemx_=cc_group(DATASET,"nr","Polymorphism",boundary,mutation,valid,seqs)
	for key in changemx_:		
		for key2 in changemx_[key]:
			changemx[key][key2].append(changemx_[key][key2])	
print("hydrophobic;"+str(np.mean(changemx['hydrophobic']['hydrophobic']))+";"+str(np.mean(changemx['hydrophobic']['positive']))+";"+str(np.mean(changemx['hydrophobic']['negative']))+";"+str(np.mean(changemx['hydrophobic']['other']))+";hydrophobic;"+str(np.std(changemx['hydrophobic']['hydrophobic']))+";"+str(np.std(changemx['hydrophobic']['positive']))+";"+str(np.std(changemx['hydrophobic']['negative']))+";"+str(np.std(changemx['hydrophobic']['other'])))
print("positive;"+str(np.mean(changemx['positive']['hydrophobic']))+";"+str(np.mean(changemx['positive']['positive']))+";"+str(np.mean(changemx['positive']['negative']))+";"+str(np.mean(changemx['positive']['other']))+";positive;"+str(np.std(changemx['positive']['hydrophobic']))+";"+str(np.std(changemx['positive']['positive']))+";"+str(np.std(changemx['positive']['negative']))+";"+str(np.std(changemx['positive']['other'])))
print("negative;"+str(np.mean(changemx['negative']['hydrophobic']))+";"+str(np.mean(changemx['negative']['positive']))+";"+str(np.mean(changemx['negative']['negative']))+";"+str(np.mean(changemx['negative']['other']))+";negative;"+str(np.std(changemx['negative']['hydrophobic']))+";"+str(np.std(changemx['negative']['positive']))+";"+str(np.std(changemx['negative']['negative']))+";"+str(np.std(changemx['negative']['other'])))
print("other;"+str(np.mean(changemx['other']['hydrophobic']))+";"+str(np.mean(changemx['other']['positive']))+";"+str(np.mean(changemx['other']['negative']))+";"+str(np.mean(changemx['other']['other']))+";other;"+str(np.std(changemx['other']['hydrophobic']))+";"+str(np.std(changemx['other']['positive']))+";"+str(np.std(changemx['other']['negative']))+";"+str(np.std(changemx['other']['other'])))

