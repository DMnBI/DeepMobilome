
from __future__ import absolute_import, division, print_function, unicode_literals
import argparse
import pandas as pd
import tensorflow as tf
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Conv1D
from keras.layers import MaxPooling1D
from keras.layers import AveragePooling1D
from keras.layers import Activation
from keras.layers import Dropout
from keras.layers import Flatten
from keras.models import Model
from keras.layers import Input
from keras.layers import BatchNormalization
from keras.layers import concatenate
from keras.callbacks import CSVLogger
from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint
from sklearn.model_selection import KFold
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
import re
from io import StringIO
import pickle
import numpy as np

from keras import backend as K
from tensorflow.keras.utils import Sequence
from matplotlib import pyplot as plt
import math
import logging, sys, subprocess,os

MAXLENGTH=8000

#NEED TO PRINT ITS PARAMETERS
#Draw Accuracy and Loss graph

class Dataloader(Sequence):
	def __init__(self,dataset,batch_size,shuffle=False):
		self.dataset = dataset
		self.batch_size = batch_size
		self.shuffle=shuffle
		self.on_epoch_end()

	def __getitem__(self,idx):
		indices = self.indices[idx*self.batch_size:(idx+1)*self.batch_size]
		filedata = list(np.array(self.dataset)[indices])
		batch_data = read_data(filedata)
		batch_x, batch_y = extend_data(batch_data)
		return np.array(batch_x),np.array(batch_y)

	def on_epoch_end(self):
		self.indices = np.arange(len(self.dataset))
		if self.shuffle ==True:
			np.random.shuffle(self.indices)
	def __len__(self):
		return math.floor(len(self.dataset) / self.batch_size)

def makeLog(logFile):
	logging.basicConfig(
		level=logging.DEBUG,
		format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
		filename=logFile,
		filemode='w',
		encoding ="utf-8"
	)
	
	return logging
def read_data(infile):

	num =0
	data = pd.DataFrame()
	for line in infile:
		if num ==0:
			data = pd.read_pickle(line)
			num=1
		else:
				
			tdf = pd.read_pickle(line)
			data = data.append(tdf, ignore_index=True)
	data = data.loc[:, ~data.columns.str.contains('^Unnamed')]
	return data


def draw_graph(history):
	train_history = history.history["loss"]
	validation_history = history.history["val_loss"]
	fig = plt.figure(figsize=(8, 8))
	plt.title("Loss History")
	plt.xlabel("EPOCH")
	plt.ylabel("LOSS Function")
	plt.plot(train_history, "red")
	plt.plot(validation_history, 'blue')
	fig.savefig("train_history.png")

	train_history = history.history["accuracy"]
	validation_history = history.history["val_accuracy"]
	fig = plt.figure(figsize=(8, 8))
	plt.title("Accuracy History")
	plt.xlabel("EPOCH")
	plt.ylabel("Accuracy")
	plt.plot(train_history, "red")
	plt.plot(validation_history, 'blue')
	fig.savefig("accuracy_history.png")


class StreamToLogger(object):
	"""
	Fake file-like stream object that redirects writes to a logger instance.
	"""
	def __init__(self, logger, log_level=logging.INFO):
		self.logger = logger
		self.log_level = log_level
		self.linebuf = ''
		def write(self, buf):
			for line in buf.rstrip().splitlines():
				self.logger.log(self.log_level, line.rstrip())

def extend_data(df):
	dataset = df.copy()
	if not dataset.empty:
		dataset = dataset.drop(columns=["Depth1", "Depth2","Depth3", "Answer","Sample"])
		count = dataset.join(pd.DataFrame(df["Depth3"].tolist(), columns=['D3_'+str(i) for i in range(MAXLENGTH)], index=df.index))

		depth1 = dataset.join(pd.DataFrame(df["Depth1"].tolist(), columns=['D1_'+str(i) for i in range(MAXLENGTH)], index=df.index))
		depth2 = dataset.join(pd.DataFrame(df["Depth2"].tolist(), columns=['D2'+str(i) for i in range(MAXLENGTH)], index=df.index))


		ans =  df.loc[:, ["Answer"]]
	else:
		depth1 =pd.DataFrame()
		depth2 =pd.DataFrame()
		count =pd.DataFrame()
		ans = pd.DataFrame()

	count = count.to_numpy()
	depth1 =depth1.to_numpy()
	
	depth2 = depth2.to_numpy()
	ans = ans.to_numpy()

	rgb = np.array([depth1, depth2, count])
	rgb = np.transpose(rgb)
	rgb = np.moveaxis(rgb,1,0)
	return	rgb, ans





#get list of plk
def init():
	parser = argparse.ArgumentParser(
	"""
	Read plk list to open plk files and read the data.
	Categorize the reads
	Make Depth calculation
	Train CNN model

	Input: plk List
	Output: Trained CNN model 
	{outPre}_model.h5, {outPre}_weights.h5
	{outPre}_loss.log, {outPre}.log

	"""
	)
	parser.add_argument("-i", "--inFile", dest="inFile",
						metavar="FILE", required=True, help="input plk list")
	parser.add_argument("-in","--inDir",dest = "inDir",metavar = "DIR",help="input plk directory")
	parser.add_argument("-o", "--outPre", dest="outPre", help="Output Prefix ")
	
	args = parser.parse_args()

	return args


def create_cnn(width, depth,mylog,fnnLayer = (128), filters=(32, 64, 128), regress=False):
	inputShape = (width, depth)
	chanDim = -1
	inputs = Input(shape=inputShape)

	filters = (32,64,128)
	fnnLayer= (128,8)


	mylog.info(f"CNN filter layers {filters}")
	mylog.info(f"fnn Layer {fnnLayer}")

	for (i, f) in enumerate(filters):
		if i == 0:
			x = inputs
		x = Conv1D(f, 3, padding="same")(x)
		x = Activation("relu")(x)
		x = BatchNormalization(axis=chanDim)(x)
		x = AveragePooling1D(pool_size=2)(x)
	x = Flatten()(x)

	for (i,f) in enumerate(fnnLayer):
		x = Dense(f)(x)
		x = Activation("relu")(x)
		x = BatchNormalization(axis=chanDim)(x)
		x = Dropout(0.5)(x)


	x = Dense(1, activation="sigmoid")(x)
	model = Model(inputs, x)
#	model.summary()
	return model

def get_fileList(inDir,inFile):
	fileList = []
	f= open(inFile)
	for line in f:
		line = line.rstrip('\r\n')
		fileList.append(inDir+'/'+line)
	f.close()
	return fileList

def main():

	gpus = tf.config.experimental.list_physical_devices('GPU')
	if gpus:
		try:
			tf.config.experimental.set_virtual_device_configuration(gpus[0], [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=1024*15)])
		except RuntimeError as e:
			print(e)

	args = init()
	epoch =5000
	batch_size = 1 #single file 256 cases
	random_state =27

	dataset2 = get_fileList(args.inDir,args.inFile)
	trainFilelist,testFilelist = train_test_split(dataset2,random_state=random_state,test_size=0.1)
	trainFilelist,valFilelist = train_test_split(trainFilelist,random_state=random_state, test_size=0.1)

	train_loader = Dataloader(trainFilelist,batch_size,shuffle=True)
	val_loader = Dataloader(valFilelist,batch_size)
	test_loader = Dataloader(testFilelist,batch_size)

#checkfile = f"{args.outPre}_checkpoint-{epoch:02d}.h5"

#	ch_checkpoint = ModelCheckpoint(checkfile,monitor='val_loss',verbose=1,save_best_only = False,save_weights_only=True,mode='auto',save_freq='epoch')
	
	del dataset2 

	mylog = makeLog(f"{args.outPre}.log")
	mylog.info("Current Directory : "+ str(os.getcwd()).rstrip('\n'))
	mylog.info ("Used Command : python "+' '.join(sys.argv))
	mylog.info (f"MAX epoch : {epoch}")
	mylog.info(f"Batch size : {batch_size}")
	mylog.info(f"Random State: {random_state}")



	csv_logger = CSVLogger(f"{args.outPre}_loss.log",append=True,separator=';')

	model = create_cnn(MAXLENGTH, 3,mylog, regress=False)

	callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=5)
	model.compile(loss="binary_crossentropy", optimizer= "adam", metrics = ["accuracy"])


	history = model.fit(train_loader,validation_data = val_loader ,epochs = epoch,callbacks=[csv_logger, callback])
	print('------------------------------------------------------------------------')
	test_loss,test_acc =	model.evaluate(test_loader)
	print(f"test loss, test acc: {test_loss}, {test_acc}")
	mylog.info(f"test loss, test acc: {test_loss}. {test_acc}")


	model.save(args.outPre+'_model.h5')
	model.save_weights(args.outPre+'_weights.h5')
	del model



	K.clear_session()

	mylog.info("END")

if __name__ == "__main__":
	main()
