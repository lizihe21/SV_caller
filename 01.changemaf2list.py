import click
import sys
@click.command()
@click.option('--maffile', default="caiji", help='mfa file')
@click.option('--speciesfile', default="caiji", help='species file')
@click.option('--outfile', default="caiji", help='out list file')
@click.option('--sp1addlc', default="caiji", help='othersp list file')
def argin(maffile, speciesfile, outfile,sp1addlc):
	"""python $0 --maffile maf --speciesfile speciesfile --outfile out."""
	partopen(maffile, speciesfile, outfile,sp1addlc)
def partopen(partfile,speciesfile,outfile,sp1addlc):
	specieslist=[]
	with open(speciesfile,'r') as f:
		for i in f:
			specieslist.append(str(i.strip()))
	dictfile={}
	with open(partfile,'r') as f:
		for i in f:
			if i.startswith("s"):
				refspecies=i.split()[1].split(".")[0]
				chrref=i.split()[1].split(".")[1]
				break
	count=0
	with open(partfile,'r') as f:
		for i in f:
			if i.startswith("s"):
				dictfile[count]=[i.split()[1].split(".")[0],i.split()[1].split(".")[1],int(i.split()[2])+1,i.split()[4],int(i.split()[5]),i.split()[-1].strip()]
				count = count + 1 #species chr starttmp fangxiang length seq
	startblock,endblock=0,0
	ax = open(outfile,'w')
	speciestmp = "chr" + "\t" + "pos" + "\t"
	for i in specieslist:
		speciestmp = speciestmp + i + "\t"
	ax.write("%s%s_loc\n"%(speciestmp,sp1addlc))
	for i in range(count):
		filepart=[]
		if i ==0:
			continue
		if refspecies == dictfile[i][0]:
			startblock=endblock
			endblock=i
			for j in range(startblock,endblock):
				filepart.append(dictfile[j])
			ax1 = partinanchors(filepart,refspecies,outfile,specieslist,sp1addlc)
			for i in ax1:
				ax.write(i)
	filepart=[]
	for j in range(endblock,count):
		filepart.append(dictfile[j])
		ax1 = partinanchors(filepart,refspecies,outfile,specieslist,sp1addlc)
	for i in ax1:
		ax.write(i)
	ax.close()
def partinanchors(filepart,refspecies,outfile,specieslist,sp1addlc):
	#species chr starttmp seq 
	#species chr starttmp fangxiang length seq
	start=filepart[0][2]
	location=[]
	string=filepart[0][5]
	species=[]
	for i in filepart:
		species.append(i[0])
	count=0
	flag = []
	for i in string:
		if i == "-":
			location.append("-")
		else:
			location.append(str(start+count))
			count = count + 1
	flag2=writeout(filepart,location,species,specieslist,sp1addlc)
	return flag2
def writeout(filepart,location,species,specieslist,sp1addlc):
	flag = []
	tmp=""
	CHR=filepart[0][1]
	countout=0
	for i in range(len(location)):
		tmp = CHR + "\t"
		tmp = tmp + str(location[i])+"\t"
		for j in specieslist:
			if j in species:
				tmp = tmp + filepart[species.index(j)][5][i].upper()+"\t"
				continue
			else:
				tmp = tmp + "X"+"\t"
		if sp1addlc in species:
			if filepart[species.index(sp1addlc)][3]=="+":
				if filepart[species.index(sp1addlc)][5][i].upper() == "-":
					tmp = tmp + filepart[species.index(sp1addlc)][1] + str("|")+"-"
				else:
					tmp = tmp + filepart[species.index(sp1addlc)][1] + str("|")+str(int(filepart[species.index(sp1addlc)][2])+countout)
					countout = countout + 1
			else:
				if filepart[species.index(sp1addlc)][5][i].upper() == "-":
					tmp = tmp + filepart[species.index(sp1addlc)][1] + str("|")+"-"
				else:#species chr starttmp fangxiang length seq
					tmp = tmp + filepart[species.index(sp1addlc)][1] + str("|")+ str(int(filepart[species.index(sp1addlc)][4])-int(filepart[species.index(sp1addlc)][2])-countout)
					countout = countout + 1
		else:
			tmp=tmp+"Nonesp"
		tmp = tmp + "\n"
		flag.append(tmp)
	return flag
if __name__ == '__main__':
	argvs=argin()
