###IMPORT
import sys
###FILES
#Add human reference proteome fasta from Uniprot to PROTEOME
PROTEOME=""
BOUNDARY="cc_prediction_boundaries.csv"
MUTATION="humsavar_mutations.csv"
CLUSTER="protein_clusters_cdhit.csv"
###FIXED_VAR
DATASET="human"#DEFINE: Deepcoil, Marcoil, Ncoils, Paircoil or human!
try:
	FILE1=sys.argv[1]
	DATASET1=sys.argv[2]
	FILE2=sys.argv[3]
	DATASET2=sys.argv[4]		
	try:
		pr=open(FILE1)
		pr.close()
		pr=open(FILE2)
		pr.close()		
	except IOError:
		print("DEFINE dataset 1 from aa_groups_bootstrap")
		print("DEFINE dataset 2 from aa_groups_bootstrap")
		print("e.g. python3 aa_groups_significance.py Deepcoil.txt DM human.txt DM")
		sys.exit()
	bool1=False	
	if DATASET1=="DM" or DATASET1=="PM":
		bool1=True
	if DATASET2=="DM" or DATASET2=="PM": 
		bool1=True
	if bool1==False:	
		print("DEFINE dataset 1 from aa_groups_bootstrap")
		print("DEFINE dataset 2 from aa_groups_bootstrap")
		print("e.g. python3 aa_groups_significance.py Deepcoil.txt DM human.txt DM")	
except IndexError:
	print("DEFINE dataset 1 from aa_groups_bootstrap")
	print("DEFINE dataset 2 from aa_groups_bootstrap")
	print("e.g. python3 aa_groups_significance.py Deepcoil.txt DM human.txt DM")
	sys.exit()

inf=open(FILE1)

mx1=[]
mx1_d=[]
sum1=0
while 1:
	line=inf.readline()
	if line=="":
		break
	if len(line.split())>1:
		if line.split()[1][:2]==DATASET1:
			
			bool1=True
			line=inf.readline()
			line=line.split(";")
			mx1.append([])
			for i in range (1, 5):
				mx1[0].append(float(line[i].strip()))
				sum1+=mx1[len(mx1)-1][len(mx1[len(mx1)-1])-1]
			mx1_d.append([])
			for i in range (6, 10):
				mx1_d[0].append(float(line[i].strip()))		
			line=inf.readline()
			line=line.split(";")			
			mx1.append([])
			for i in range (1, 5):
				mx1[1].append(float(line[i].strip()))
				sum1+=mx1[len(mx1)-1][len(mx1[len(mx1)-1])-1]
			mx1_d.append([])
			for i in range (6, 10):
				mx1_d[1].append(float(line[i].strip()))				
			line=inf.readline()
			line=line.split(";")			
			mx1.append([])
			for i in range (1, 5):
				mx1[2].append(float(line[i].strip()))
				sum1+=mx1[len(mx1)-1][len(mx1[len(mx1)-1])-1]
			mx1_d.append([])
			for i in range (6, 10):
				mx1_d[2].append(float(line[i].strip()))				
			line=inf.readline()
			line=line.split(";")			
			mx1.append([])
			for i in range (1, 5):
				mx1[3].append(float(line[i].strip()))	
				sum1+=mx1[len(mx1)-1][len(mx1[len(mx1)-1])-1]
			mx1_d.append([])
			for i in range (6, 10):
				mx1_d[3].append(float(line[i].strip()))				
			break
	
inf.close()
for i in range (0, len(mx1)):
	for j in range (0, len(mx1)):
		mx1[i][j]=mx1[i][j]/sum1
for i in range (0, len(mx1_d)):
	for j in range (0, len(mx1_d)):
		mx1_d[i][j]=mx1_d[i][j]/sum1
inf=open(FILE2)

mx2=[]
mx2_d=[]
sum2=0
while 1:
	line=inf.readline()
	if line=="":
		break
	if len(line.split())>1:
		if line.split()[1][:2]==DATASET2:
			
			bool1=True
			line=inf.readline()
			line=line.split(";")
			mx2.append([])
			for i in range (1, 5):
				mx2[0].append(float(line[i].strip()))
				sum2+=mx2[len(mx2)-1][len(mx2[len(mx2)-1])-1]
			mx2_d.append([])
			for i in range (6, 10):
				mx2_d[0].append(float(line[i].strip()))		
			line=inf.readline()
			line=line.split(";")			
			mx2.append([])
			for i in range (1, 5):
				mx2[1].append(float(line[i].strip()))
				sum2+=mx2[len(mx2)-1][len(mx2[len(mx2)-1])-1]
			mx2_d.append([])
			for i in range (6, 10):
				mx2_d[1].append(float(line[i].strip()))				
			line=inf.readline()
			line=line.split(";")			
			mx2.append([])
			for i in range (1, 5):
				mx2[2].append(float(line[i].strip()))
				sum2+=mx2[len(mx2)-1][len(mx2[len(mx2)-1])-1]
			mx2_d.append([])
			for i in range (6, 10):
				mx2_d[2].append(float(line[i].strip()))				
			line=inf.readline()
			line=line.split(";")			
			mx2.append([])
			for i in range (1, 5):
				mx2[3].append(float(line[i].strip()))	
				sum2+=mx2[len(mx2)-1][len(mx2[len(mx2)-1])-1]
			mx2_d.append([])
			for i in range (6, 10):
				mx2_d[3].append(float(line[i].strip()))			
			break
	
inf.close()
for i in range (0, len(mx2)):
	for j in range (0, len(mx2)):
		mx2[i][j]=mx2[i][j]/sum2
for i in range (0, len(mx2_d)):
	for j in range (0, len(mx2_d)):
		mx2_d[i][j]=mx2_d[i][j]/sum2

print(FILE1+" "+DATASET1+" vs "+FILE2+" "+DATASET2+";hydrophobic;positive;negative;other")
s="hydrophobic"
for i in range (0, len(mx1[0])):
	if mx1[0][i]+mx1_d[0][i]*3<mx2[0][i]-mx2_d[0][i]*3 or mx1[0][i]-mx1_d[0][i]*3>mx2[0][i]+mx2_d[0][i]*3:
		s=s+";***"
	elif mx1[0][i]+mx1_d[0][i]*3<mx2[0][i]-mx2_d[0][i]*2 or mx1[0][i]-mx1_d[0][i]*3>mx2[0][i]+mx2_d[0][i]*2:
		s=s+";**"
	elif mx1[0][i]+mx1_d[0][i]*3<mx2[0][i]-mx2_d[0][i] or mx1[0][i]-mx1_d[0][i]*3>mx2[0][i]+mx2_d[0][i]:
		s=s+";*"
	else:
		s=s+";0"	
print(s)	
s="positive"
for i in range (0, len(mx1[1])):
	if mx1[1][i]+mx1_d[1][i]*3<mx2[1][i]-mx2_d[1][i]*3 or mx1[1][i]-mx1_d[1][i]*3>mx2[1][i]+mx2_d[1][i]*3:
		s=s+";***"
	elif mx1[1][i]+mx1_d[1][i]*3<mx2[1][i]-mx2_d[1][i]*2 or mx1[1][i]-mx1_d[1][i]*3>mx2[1][i]+mx2_d[1][i]*2:
		s=s+";**"
	elif mx1[1][i]+mx1_d[1][i]*3<mx2[1][i]-mx2_d[1][i] or mx1[1][i]-mx1_d[1][i]*3>mx2[1][i]+mx2_d[1][i]:
		s=s+";*"
	else:
		s=s+";0"	
print(s)	
s="negative"
for i in range (0, len(mx1[2])):
	if mx1[2][i]+mx1_d[2][i]*3<mx2[2][i]-mx2_d[2][i]*3 or mx1[2][i]-mx1_d[2][i]*3>mx2[2][i]+mx2_d[2][i]*3:
		s=s+";***"
	elif mx1[2][i]+mx1_d[2][i]*3<mx2[2][i]-mx2_d[2][i]*2 or mx1[2][i]-mx1_d[2][i]*3>mx2[2][i]+mx2_d[2][i]*2:
		s=s+";**"
	elif mx1[2][i]+mx1_d[2][i]*3<mx2[2][i]-mx2_d[2][i] or mx1[2][i]-mx1_d[2][i]*3>mx2[2][i]+mx2_d[2][i]:
		s=s+";*"
	else:
		s=s+";0"	
print(s)	
s="other"
for i in range (0, len(mx1[3])):
	if mx1[3][i]+mx1_d[3][i]*3<mx2[3][i]-mx2_d[3][i]*3 or mx1[3][i]-mx1_d[3][i]*3>mx2[3][i]+mx2_d[3][i]*3:
		s=s+";***"
	elif mx1[3][i]+mx1_d[3][i]*3<mx2[3][i]-mx2_d[3][i]*2 or mx1[3][i]-mx1_d[3][i]*3>mx2[3][i]+mx2_d[3][i]*2:
		s=s+";**"
	elif mx1[3][i]+mx1_d[3][i]*3<mx2[3][i]-mx2_d[3][i] or mx1[3][i]-mx1_d[3][i]*3>mx2[3][i]+mx2_d[3][i]:
		s=s+";*"
	else:
		s=s+";0"	
print(s)				
