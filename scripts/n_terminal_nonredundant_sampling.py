###IMPORT
from scipy.stats import chi2_contingency
import numpy as np
import pandas as pd
import random
###FILES
#Add human reference proteome fasta from Uniprot to PROTEOME
PROTEOME=""
BOUNDARY="cc_prediction_boundaries.csv"
MUTATION="humsavar_mutations.csv"
CLUSTER="protein_clusters_cdhit.csv"

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
def cc_residue(method,dataset,boundary_,mutation_,valid_,seqs_):
	cc_=0
	DM_=0
	PM_=0
	cc_7=0
	DM_7=0
	PM_7=0	
	for AC in seqs_:
		accept=True
		cc_p=0
		DM_p=0
		PM_p=0
		cc_p7=0
		DM_p7=0
		PM_p7=0		
		if random.random()<0.8:
			try:
				if valid_[AC]==dataset:	
					for i in range (0, len(boundary_[AC][method])):			
						if len(seqs_[AC])>boundary_[AC][method][i][1]-1:
							cc_p+=(boundary_[AC][method][i][1]-boundary_[AC][method][i][0]+1)
							cc_p7+=min((boundary_[AC][method][i][1]-boundary_[AC][method][i][0]+1),7)
							for position in mutation_[AC]:
								if boundary_[AC][method][i][0]<=position<=boundary_[AC][method][i][1]:
									for phenotype in mutation_[AC][position]:			
										for j in range (0, len(mutation_[AC][position][phenotype])):	
											if seqs_[AC][position-1]==mutation_[AC][position][phenotype][j][0]:
												if phenotype=="Polymorphism":
													PM_p+=1
												else: 
													DM_p+=1
												if position-1<=boundary_[AC][method][i][0]+7:											
													if phenotype=="Polymorphism":
														PM_p7+=1
													else: 
														DM_p7+=1
											else:
												accept=False
						else:
							accept=False
			except KeyError:
				pass
		if accept==True:
			cc_+=cc_p
			DM_+=DM_p
			PM_+=PM_p
			cc_7+=cc_p7
			DM_7+=DM_p7
			PM_7+=PM_p7			
	return [cc_-cc_7,DM_-DM_7,PM_-PM_7,cc_7,DM_7,PM_7]
			
###MAIN
print("NON REDUNDANT DATASET ALL DATA")
dc_dm7freq=[]
dc_pm7freq=[]
dc_dmfreq=[]
dc_pmfreq=[]

pc_dm7freq=[]
pc_pm7freq=[]
pc_dmfreq=[]
pc_pmfreq=[]

mc_dm7freq=[]
mc_pm7freq=[]
mc_dmfreq=[]
mc_pmfreq=[]

nc_dm7freq=[]
nc_pm7freq=[]
nc_dmfreq=[]
nc_pmfreq=[]
print(";NON REDUNDANT DATASET RANDOM SAMPLING;;;;;;;;;;;;;;;;;;")
print(";;Mean Rel. Freq (DM);Std Rel. Freq (DM);Mean Rel. Freq (PM);Std Rel. Freq (PM);;;;;;;;;;;;;;")

for sampling in range (0,100):

	[cc,DM,PM,cc7,DM7,PM7]=cc_residue("Deepcoil","nr",boundary,mutation,valid,seqs)
	dc_dm7freq.append(DM7/cc7)
	dc_pm7freq.append(PM7/cc7)
	dc_dmfreq.append(DM/cc)
	dc_pmfreq.append(PM/cc)	
	[cc,DM,PM,cc7,DM7,PM7]=cc_residue("Marcoil","nr",boundary,mutation,valid,seqs)
	mc_dm7freq.append(DM7/cc7)
	mc_pm7freq.append(PM7/cc7)
	mc_dmfreq.append(DM/cc)
	mc_pmfreq.append(PM/cc)		
	[cc,DM,PM,cc7,DM7,PM7]=cc_residue("Ncoils","nr",boundary,mutation,valid,seqs)
	nc_dm7freq.append(DM7/cc7)
	nc_pm7freq.append(PM7/cc7)
	nc_dmfreq.append(DM/cc)
	nc_pmfreq.append(PM/cc)		
	[cc,DM,PM,cc7,DM7,PM7]=cc_residue("Paircoil","nr",boundary,mutation,valid,seqs)
	pc_dm7freq.append(DM7/cc7)
	pc_pm7freq.append(PM7/cc7)
	pc_dmfreq.append(DM/cc)
	pc_pmfreq.append(PM/cc)		

total_dm7m=[]
total_dm7s=[]
total_pm7m=[]
total_pm7s=[]

total_dmm=[]
total_dms=[]
total_pmm=[]
total_pms=[]

dm7m=np.mean(dc_dm7freq)
dm7s=np.std(dc_dm7freq)
pm7m=np.mean(dc_pm7freq)
pm7s=np.std(dc_pm7freq)

total_dm7m.append(dm7m)
total_dm7s.append(dm7s)
total_pm7m.append(pm7m)
total_pm7s.append(pm7s)

dmm=np.mean(dc_dmfreq)
dms=np.std(dc_dmfreq)
pmm=np.mean(dc_pmfreq)
pms=np.std(dc_pmfreq)

total_dmm.append(dmm)
total_dms.append(dms)
total_pmm.append(pmm)
total_pms.append(pms)

print("Deepcoil;1-7 residues;"+str(dm7m)+";"+str(dm7s)+";"+str(pm7m)+";"+str(pm7s))
print(";other;"+str(dmm)+";"+str(dms)+";"+str(pmm)+";"+str(pms)+"")

dm7m=np.mean(mc_dm7freq)
dm7s=np.std(mc_dm7freq)
pm7m=np.mean(mc_pm7freq)
pm7s=np.std(mc_pm7freq)

total_dm7m.append(dm7m)
total_dm7s.append(dm7s)
total_pm7m.append(pm7m)
total_pm7s.append(pm7s)

dmm=np.mean(mc_dmfreq)
dms=np.std(mc_dmfreq)
pmm=np.mean(mc_pmfreq)
pms=np.std(mc_pmfreq)

total_dmm.append(dmm)
total_dms.append(dms)
total_pmm.append(pmm)
total_pms.append(pms)

print("Marcoil;1-7 residues;"+str(dm7m)+";"+str(dm7s)+";"+str(pm7m)+";"+str(pm7s))
print(";other;"+str(dmm)+";"+str(dms)+";"+str(pmm)+";"+str(pms)+"")


dm7m=np.mean(nc_dm7freq)
dm7s=np.std(nc_dm7freq)
pm7m=np.mean(nc_pm7freq)
pm7s=np.std(nc_pm7freq)

total_dm7m.append(dm7m)
total_dm7s.append(dm7s)
total_pm7m.append(pm7m)
total_pm7s.append(pm7s)

dmm=np.mean(nc_dmfreq)
dms=np.std(nc_dmfreq)
pmm=np.mean(nc_pmfreq)
pms=np.std(nc_pmfreq)

total_dmm.append(dmm)
total_dms.append(dms)
total_pmm.append(pmm)
total_pms.append(pms)

print("Ncoils;1-7 residues;"+str(dm7m)+";"+str(dm7s)+";"+str(pm7m)+";"+str(pm7s))
print(";other;"+str(dmm)+";"+str(dms)+";"+str(pmm)+";"+str(pms)+"")


dm7m=np.mean(pc_dm7freq)
dm7s=np.std(pc_dm7freq)
pm7m=np.mean(pc_pm7freq)
pm7s=np.std(pc_pm7freq)

total_dm7m.append(dm7m)
total_dm7s.append(dm7s)
total_pm7m.append(pm7m)
total_pm7s.append(pm7s)

dmm=np.mean(pc_dmfreq)
dms=np.std(pc_dmfreq)
pmm=np.mean(pc_pmfreq)
pms=np.std(pc_pmfreq)

total_dmm.append(dmm)
total_dms.append(dms)
total_pmm.append(pmm)
total_pms.append(pms)

print("Parcoil;1-7 residues;"+str(dm7m)+";"+str(dm7s)+";"+str(pm7m)+";"+str(pm7s))
print(";other;"+str(dmm)+";"+str(dms)+";"+str(pmm)+";"+str(pms)+"")
print("Mean;1-7 residues;"+str(np.mean(total_dm7m))+";"+str(np.mean(total_dm7s))+";"+str(np.mean(total_pm7m))+";"+str(np.mean(total_pm7s)))
print(";other;"+str(np.mean(total_dmm))+";"+str(np.mean(total_dms))+";"+str(np.mean(total_pmm))+";"+str(np.mean(total_pms)))


