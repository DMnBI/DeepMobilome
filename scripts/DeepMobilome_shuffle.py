import argparse
import pandas as pd
import random

def init():
	parser = argparse.ArgumentParser(
	"""
	Read plk files and Shuffle data and Write new plk files

	"""
	)
	parser.add_argument("-i", "--inFile", dest="inFile",
						metavar="FILE", required=True, help="input plk list")
	parser.add_argument("-in","--inDir",dest = "inDir",metavar = "DIR",help="input plk directory")
	parser.add_argument("-o", "--outPre", dest="outPre", help="Output prefix file")
	parser.add_argument("-b","--batch",dest = "batch",type = int ,help ="batchSize")
	parser.add_argument("-s","--sfl",dest= "sfl",type=int, help ="shuffle file number")
	
	args = parser.parse_args()

	return args




def read_data(fileList,batchsize,out,sfl):

	
	num =0 
	data = pd.DataFrame()
	outNum =0 
	while num < len(fileList):
		if num%sfl==0:
			
			if num!=0:

				dataLen = len(data)
#				print("In:",dataLen)
				data = data.sample(frac=1).reset_index(drop=True)
				while dataLen>batchsize:
				
					data.iloc[:batchsize].to_pickle(out+'_'+str(outNum)+".plk")
					
					data = data.iloc[batchsize:]
					dataLen = len(data)
					outNum+=1

			print(fileList[num])

			data = pd.read_pickle(fileList[num])

		else:
			print(fileList[num])	
			tdf = pd.read_pickle(fileList[num])
			data = data.append(tdf, ignore_index=True)
		num+=1
	dataLen = len(data)
	data = data.sample(frac=1).reset_index(drop=True)
	while dataLen>batchsize:
#		print("Out:",dataLen)	
		data.iloc[:batchsize].to_pickle(out+'_'+str(outNum)+".plk")
		
		data = data.iloc[batchsize:]
		dataLen = len(data)
		outNum+=1


		#data=pd.DataFrame()
		#empty=0






def read_file(indir,infile):
	infileList =[]
	f = open(infile)
	for line in f:
		line = line.rstrip('\r\n')
		infileList.append(indir+'/'+line)

	f.close()
	random.shuffle(infileList)
	return infileList

def main():
	args = init()

	random.seed(10)
	batch_size = args.batch

	fileList = read_file(args.inDir,args.inFile)
	read_data(fileList,batch_size,args.outPre,args.sfl)

if __name__ == "__main__":
	main()
