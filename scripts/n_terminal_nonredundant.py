###IMPORT
from scipy.stats import chi2_contingency
import numpy as np
import pandas as pd
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
print(";;;;;;;;;;;;DM;PM;Odds ratio")
[cc,DM,PM,cc7,DM7,PM7]=cc_residue("Deepcoil","nr",boundary,mutation,valid,seqs)
odds_a=[]
dm_a=[]
pm_a=[]
dm7_a=[]
pm7_a=[]
mut={"1-7":{"DM":DM7,"PM":PM7},"rest":{"DM":DM,"PM":PM}}
mut2=pd.DataFrame.from_dict(mut)
print("Deepcoil;1-7;"+str(cc7)+";"+str(DM7)+";"+str(PM7)+";"+str(chi2_contingency(mut2)[1])+";"+str((DM7/DM)/(PM7/PM))+";"+str(DM7/cc7)+";"+str(PM7/cc7))
print("Deepcoil;other;"+str(cc)+";"+str(DM)+";"+str(PM)+";;;"+str(DM/cc)+";"+str(PM/cc))
odds_a.append((DM7/DM)/(PM7/PM))
dm_a.append(DM/cc)
pm_a.append(PM/cc)
dm7_a.append(DM7/cc7)
pm7_a.append(PM7/cc7)
[cc,DM,PM,cc7,DM7,PM7]=cc_residue("Marcoil","nr",boundary,mutation,valid,seqs)
mut={"1-7":{"DM":DM7,"PM":PM7},"rest":{"DM":DM,"PM":PM}}
mut2=pd.DataFrame.from_dict(mut)
print("Marcoil;1-7;"+str(cc7)+";"+str(DM7)+";"+str(PM7)+";"+str(chi2_contingency(mut2)[1])+";"+str((DM7/DM)/(PM7/PM))+";"+str(DM7/cc7)+";"+str(PM7/cc7))
print("Marcoil;other;"+str(cc)+";"+str(DM)+";"+str(PM)+";;;"+str(DM/cc)+";"+str(PM/cc))
mut={"1-7":{"DM":DM7,"PM":PM7},"rest":{"DM":DM,"PM":PM}}
mut2=pd.DataFrame.from_dict(mut)
odds_a.append((DM7/DM)/(PM7/PM))
dm_a.append(DM/cc)
pm_a.append(PM/cc)
dm7_a.append(DM7/cc7)
pm7_a.append(PM7/cc7)
[cc,DM,PM,cc7,DM7,PM7]=cc_residue("Ncoils","nr",boundary,mutation,valid,seqs)
mut={"1-7":{"DM":DM7,"PM":PM7},"rest":{"DM":DM,"PM":PM}}
mut2=pd.DataFrame.from_dict(mut)
print("Ncoils;1-7;"+str(cc7)+";"+str(DM7)+";"+str(PM7)+";"+str(chi2_contingency(mut2)[1])+";"+str((DM7/DM)/(PM7/PM))+";"+str(DM7/cc7)+";"+str(PM7/cc7))
print("Ncoils;other;"+str(cc)+";"+str(DM)+";"+str(PM)+";;;"+str(DM/cc)+";"+str(PM/cc))
mut={"1-7":{"DM":DM7,"PM":PM7},"rest":{"DM":DM,"PM":PM}}
mut2=pd.DataFrame.from_dict(mut)
odds_a.append((DM7/DM)/(PM7/PM))
dm_a.append(DM/cc)
pm_a.append(PM/cc)
dm7_a.append(DM7/cc7)
pm7_a.append(PM7/cc7)
[cc,DM,PM,cc7,DM7,PM7]=cc_residue("Paircoil","nr",boundary,mutation,valid,seqs)
mut={"1-7":{"DM":DM7,"PM":PM7},"rest":{"DM":DM,"PM":PM}}
mut2=pd.DataFrame.from_dict(mut)
print("Parcoil;1-7;"+str(cc7)+";"+str(DM7)+";"+str(PM7)+";"+str(chi2_contingency(mut2)[1])+";"+str((DM7/DM)/(PM7/PM))+";"+str(DM7/cc7)+";"+str(PM7/cc7))
print("Parcoil;other;"+str(cc)+";"+str(DM)+";"+str(PM)+";;;"+str(DM/cc)+";"+str(PM/cc))
mut={"1-7":{"DM":DM7,"PM":PM7},"rest":{"DM":DM,"PM":PM}}
mut2=pd.DataFrame.from_dict(mut)
odds_a.append((DM7/DM)/(PM7/PM))
dm_a.append(DM/cc)
pm_a.append(PM/cc)
dm7_a.append(DM7/cc7)
pm7_a.append(PM7/cc7)
print("Mean values;1-7;;;;;"+str(np.mean(odds_a))+";"+str(np.mean(dm7_a))+";"+str(np.mean(pm7_a)))
print(";other;;;;;;"+str(np.mean(dm_a))+";"+str(np.mean(pm_a)))



