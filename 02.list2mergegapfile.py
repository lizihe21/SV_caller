import click
import os
import sys
import re
from collections import Counter
import numpy as np
import pandas as pd
@click.command()
@click.option('--maffile', default="maffile", help='maffile')
@click.option('--out', default="outfile", help='outputfile')
@click.option('--outgroup', default="groupfile", help='outgroupfile')
def argin(maffile, out, outgroup):
	"""pypy 00.change2list.py --maffile merge.dog.human.filter.sort.v1.maf.list --speciesfile merge.dog.human.filter.sort.v1.maf.test --outfile testlistfile"""
	mainscript(maffile, out, outgroup)
def readspeciesandmaf(lf): # changetopandasfile format
	matrixall=pd.read_csv(lf,sep="\t",dtype=str)
	return matrixall
def readbedregion(pdmaf,outgroup):
	df = pdmaf
	H=[]
	outgrouplist=[]
	out=[]
	listoutsel=[]
	AA=["-","X"]
	with open(outgroup,'r') as f:
		for i in f:
			outgrouplist.append(i.strip())
	ingrouplist=[]
	headname_1=pdmaf.columns.values[2:-1]
	rawspecies_2=pdmaf.columns.values[-1]
	for i in headname_1:
		if i in outgrouplist:
			continue
		ingrouplist.append(i)
	species1_sm=ingrouplist[0]
	species2_sm=ingrouplist[-1]
	try:
		df2=df.replace("X","-")
		df3=df2[outgrouplist]
		tmpoutgroup1=[]
		for row in df3.itertuples():
			tmpoutgroup1.append(''.join(list(row)[1:]))
		df2.loc[:,"outgroupmerge"]=tmpoutgroup1
		mathchenlen="-"*len(outgrouplist)
		listoutsel=df2[df2["outgroupmerge"] == mathchenlen].index.tolist()
	except IndexError:
		sys.exit()
	else:
		location=[]
		for i in listoutsel:
			location.append([i,i+1])
		intervals = sorted(location,key=lambda x: x[0])
		merged = []
		length = len(df)
		for interval in intervals:
			if not merged or merged[-1][-1] < interval[0]:
				merged.append(interval)
			else:
				merged[-1][-1] = max(merged[-1][-1], interval[-1])
		for i in merged:
			start,end=getlistback(df,i,length)
			start_out2,end_out2=getlistbackout2(df,i,length,rawspecies_2)
			if start == -1 or end == -1:
				continue
			if i[1] - i[0] > 4:
				dfingroup=df.loc[i[0]:i[1]-1,ingrouplist]
				dfeslecttmp1=df.loc[i[0]:i[1]-1,[species1_sm]]
				dfeslecttmp2=df.loc[i[0]:i[1]-1,[species2_sm]]
				if checkpass(dfeslecttmp1) == 1 or checkpass(dfeslecttmp2) == 1:
					continue
				identityin=calculidentitynew(dfingroup)
				dfinup=df.loc[i[0]-50:i[0]-1,ingrouplist]
				dfindown=df.loc[i[1]:50+i[1]-1,ingrouplist]
				dfoutup=df.loc[i[0]-50:i[0]-1,outgrouplist]
				dfoutdown=df.loc[i[1]:50+i[1]-1,outgrouplist]
				identityinup=calculidentitynew(dfinup)
				identityindown=calculidentitynew(dfindown)
				identityoutup=calculidentitynew(dfoutup)
				identityoutdown=calculidentitynew(dfoutdown)
				length_set=len(dfingroup)
				tmpout=[df.loc[i[0],"chr"],start,end,start_out2,end_out2,identityin,identityinup,identityindown,identityoutup,identityoutdown,length_set]
				out.append(tmpout)
	return out
def checkpass(A):
	flagtmp1=0
	for row in A.itertuples():
		if list(row)[1:] == ["-"]:
			flagtmp1 = 1
		else:
			flagtmp1 = 0
			break
	if flagtmp1 == 0:
		return 0
	else: 
		return 1 
def getlistback(A,B,C):
	df = A
	start,end=-1,-1
	for i in range(20):
		if B[0]-i <=0 or B[0]-i >=C-2:
			return [start,end]
		if df.loc[B[0]-i,"pos"] != "-":
			start=df.loc[B[0]-i,"pos"]
			break
	for i in range(20):
		if B[1]+i >=C-2 or B[1]+i <=0:
			return [start,end]
		if df.loc[B[1]+i,"pos"] != "-":
			end = df.loc[B[1]+i,"pos"]
			break
	return [start,end]
def getlistbackout2(A,B,C,D):
	df = A
	start,end=-1,-1
	for i in range(40):
		if B[0]-i+1 <=0 or B[0]-i+1 >=C-2:
			return [start,end]
		if df.loc[B[0]-i+1,D] != "Nonesp" and "-" not in df.loc[B[0]-i+1,D]:
			start=df.loc[B[0]-i+1,D]
			break
	for i in range(40):
		if B[1]+i >=C-2 or B[1]+i <=0:
			return [start,end]
		if df.loc[B[1]+i,D] != "Nonesp" and df.loc[B[1]+i,D].split("|")[1] != "-":
			end = df.loc[B[1]+i,D]
			break
	return [start,end]
def calculidentitynew(dflist):
	tmp=[]
	length=len(dflist.columns)
	count = 0
	out = 0
	for row in dflist.itertuples():
		tmplist=list(row)[1:]
		i=tmplist
		most=float(Counter(i).most_common(1)[0][1])*1.0/len(i)
		base=Counter(i).most_common(1)[0][0]
		if str(base) == "-" or str(base) == "X" or str(base) == "N":
			most = 0
		count =  count + 1
		out = out + most
	identityout=float(out*1.0/count)
	return identityout
def getout(filename,outbed):
	ax = open(filename,'w')
	for i in outbed:
		ax.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10]))
	ax.close()
def mainscript(maffile, out, outgroup): # main script
	pandasformastfile=readspeciesandmaf(maffile)
	outfile=readbedregion(pandasformastfile,outgroup)
	if len(outfile) == 0:
		pass
	else:
		getout(out,outfile)
if __name__ == '__main__':
	argvs=argin()
