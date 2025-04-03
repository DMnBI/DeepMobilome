#!/usr/bin/env python3

import sys
import argparse
import subprocess
import math
import statistics
import logging
import re
import pandas as pd
import tensorflow as tf
import numpy as np

from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Conv1D
from keras.layers import AveragePooling1D
from keras.layers import Activation
from keras.layers import Dropout
from keras.layers import Flatten
from keras.models import Model
from keras.layers import Input
from keras.layers import BatchNormalization
from keras.layers import concatenate
from sklearn.model_selection import KFold
from sklearn.metrics import classification_report
from tensorflow.keras.utils import plot_model

from DeepMobilome_common import makeLog
from DeepMobilome_common import ref_dic_realdata
from DeepMobilome_common import sam_parse
from DeepMobilome_common import save_read

from collections import defaultdict
from numpy import loadtxt
from io import StringIO

#Upgrade
#Point - single sam file -> single output file


MAX_LENGTH = 8000
mol = ['A','C','G','T','N']
Sample = []
Sum = []
Depth1 =[]
Depth2 =[]
Depth3 = []
Count = []
Ans = []


def init_argv():
    parser = argparse.ArgumentParser(description="""
                    Predict Target

                    Input: Model weight, SAM file
                    Output: Predicted result

                    """ 
                ,formatter_class= argparse.RawTextHelpFormatter)
    parser.add_argument("-i","--inFile",dest="inFile",metavar="FILE",required=True, help="sam File : sample reads aligned to target sequence")
    parser.add_argument("-o","--outFile",dest="outFile",metavar="FILE",required=True,help="predicted result / samplename,Answer,Predict,AvgDepth / csv file")
    parser.add_argument("-l","--logFile",dest="logFile",metavar="FILE",required=True,help="log file")
    parser.add_argument("-r","--ref",dest="ref",metavar="FILE",required=True, help="target sequence file / fna file")
    parser.add_argument("-c","--count",dest = "count",required=True, help ="read count / samplename,read count / csv file",)
    parser.add_argument("-a","--ans",dest= "ans",metavar ="FILE", help = "answer for the given case / samplename,answer / csv file",required=True)
    parser.add_argument("-s","--save", dest="save",help="model weight", required=True)


    args = parser.parse_args()
    return args

def extend_data(df):
    dataset = df.copy()

    dataset = dataset.drop(
        columns=["Depth1", "Depth2","Depth3", "Answer","Sample"])

    count = dataset.join(pd.DataFrame(df["Depth3"].tolist(), columns=[
                         'D3_'+str(i) for i in range(len(df["Depth3"][0]))], index=df.index))
    

    depth1 = dataset.join(pd.DataFrame(df["Depth1"].tolist(), columns=[
                          'D1_'+str(i) for i in range(len(df["Depth1"][0]))], index=df.index))
    depth2 = dataset.join(pd.DataFrame(df["Depth2"].tolist(), columns=[
                          'D2'+str(i) for i in range(len(df["Depth2"][0]))], index=df.index))

    for i in range(len(df["Depth1"][0]),MAX_LENGTH):
        count["D3_"+str(i)]=-5
        depth1["D1_"+str(i)]=-5
        depth2["D2_"+str(i)]=-5
    count = count.to_numpy()
    depth1 =depth1.to_numpy()
    depth2 = depth2.to_numpy()

    ans = df.loc[:, ["Answer"]].to_numpy()

    rgb = np.array([depth1, depth2, count])
    rgb = np.transpose(rgb)
    rgb = np.moveaxis(rgb,1,0)
    return    rgb, ans

def convert_to_0(predict):

    return list(map(lambda x: 0 if x < 0.5 else 1, predict))




def create_cnn(width, depth, filters=(32, 64, 128), regress=False):
    inputShape = (width, depth)
    chanDim = -1
    fnnLayer = (128,8)

    inputs = Input(shape=inputShape)
    
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
#    model.summary()
#    plot_model(model,to_file="model.png",show_shapes=True)    
    return model


def calculate_depth_coverage(key,refDic,countRead,refRead,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth,refLengthDic):
    mappedlength = 0
    mappeddepth = 0
    sample_dom = 0 
    t_entrophy =0
    linklengthfe =0 
    linklengthmid =0 
    linkdepthfe =0
    linklengthmid =0

    entrophy = 0
    sum_read =0
    mismatchpoint = 0
    steeppoint =0
    for k,v in refDic.items():
        refitems = k.split('_')
        reflength=len(refDic[k])
        realLength = refLengthDic[k]
        err =500
        refmid = reflength -(err*2)
        reffe = err*2
        mappedlength=0
        mappeddepth=0
        linklengthfe =0
        linklengthmid =0
        linkdepthfe =0
        linkdepthmid =0

        sample_dom =0
        t_entrophy=0
        mismatchpoint =0
        steeppoint =0
        depth_list = []

        tdepth1=[]
        tdepth2=[]
        tcount=[]
        tdepth3=[]
        for li,lj in linkmatrix[k].items():
            if li !=reflength:
                if len(readposcount[k][li] & readposcount[k][li+1])==0 and len(readposcount[k][li])!=0 and len(readposcount[k][li+1])!=0:
                    mismatchpoint+=1
                    if abs(len(readposcount[k][li])-len(readposcount[k][li+1]))>=10:
                        steeppoint+=1
            if (li < err or li>= reflength-err ) and lj !=0:
                linklengthfe +=1
                linkdepthfe +=lj
            elif lj!=0:
                linklengthmid +=1
                linkdepthmid+=lj
        for i,j in v.items():
            if j >0:
                mappedlength+=1
                mappeddepth +=j
                depth_list.append(j)
                entrophy=0
                for dna in mol:
                    if  sampleDNA[k][i][dna]!= 0 and j>=0:
                        entrophy+= -(float(sampleDNA[k][i][dna])/(j+1))* math.log(float(sampleDNA[k][i][dna])/(j+1),2)
                t_entrophy +=entrophy
            if refDNA[k][i] != max(sampleDNA[k][i],key = sampleDNA[k][i].get) and sampleDNA[k][i][max(sampleDNA[k][i],key = sampleDNA[k][i].get)]!=0:
                sample_dom +=1
        

        sum_read +=refRead[k]
        if len(depth_list)==0:
            median_depth =0
        else:
            median_depth = statistics.median(depth_list)
        Sample.append(key+"_"+k)
        countkey = key.split('_')[0] 
        #countRead key needs to be same with the sam list name


        Sum.append((float(countRead[countkey]),float(mappedlength)/float(realLength),float(mappeddepth)/float(realLength),((float(mappeddepth)/float(realLength))/float(countRead[countkey]))*float(10**8),(float(refRead[k])*float(10**3*10**6))/(float(countRead[countkey])*float(realLength)),float(float(median_depth)/float(countRead[countkey]))*float(10**9),float(linklengthfe)/float(reffe),float(linkdepthfe)/float(reffe),float(linklengthmid)/float(refmid),float(linkdepthmid)/float(refmid),float(linklengthfe+linklengthmid)/float(realLength),float(linkdepthfe+linkdepthmid)/float(realLength),mismatchpoint,steeppoint,mappedlength))
        for refnum in range(1,MAX_LENGTH+1):
            if refRead[k]>0:
                if v[refnum]>=0:

                    tdepth1.append(float(v[refnum])/float(refRead[k])*float(realLength))
                    tdepth3.append(float(insertdepth[k][refnum])/float(refRead[k])*float(realLength))
                    if refnum!= reflength:
                        tdepth2.append(float(linkmatrix[k][refnum])/float(refRead[k])*float(realLength))
                        tcount.append(len(readposcount[k][refnum] & readposcount[k][refnum+1]))
                    else:
                        tdepth2.append(float(linkmatrix[k][refnum])/float(refRead[k])*float(realLength))
                        tcount.append(len(readposcount[k][refnum] & readposcount[k][refnum-1]))
                else:
                    tdepth1.append(float(v[refnum]))
                    tdepth3.append(float(insertdepth[k][refnum]))
                    if refnum!= reflength:
                        tdepth2.append(float(linkmatrix[k][refnum]))
                        tcount.append(len(readposcount[k][refnum] & readposcount[k][refnum+1]))
                    else:
                        tdepth2.append(float(linkmatrix[k][refnum]))
                        tcount.append(len(readposcount[k][refnum] & readposcount[k][refnum-1]))


            else: 

                tdepth1.append(0)
                tdepth3.append(0)
                if refnum!= reflength:
                    tdepth2.append(0)
                    tcount.append(0)
                else:
                    tdepth2.append(0)
                    tcount.append(0)

    Depth1.append(tuple(tdepth1))
    Depth2.append(tuple(tdepth2))
    Depth3.append(tuple(tdepth3))
    Count.append(tuple(tcount))
    print(f"[Complete] \"{key}\" Aligned read # {sum_read}")

def read_ans(ansfile):
    ansDic = {}
    f = open(ansfile)
    num =0
    for line in f:
        line = line.rstrip('\r\n')
        if num!=0:
            items = line.split(',')
            ansDic[items[0]] = int(items[1])
        num+=1
    return ansDic

def dataTransformation(inFile, refFile, countFile, ansFile, nm, sr):

    countRead = save_read(countFile)
    ansDic = read_ans(ansFile)

#fList[num] = inFile
    key = inFile.split('.sam')[0]
    key = key.split('/')[-1]
    name = key.split('_')[0]
    if name in ansDic:
        refDic,refRead,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth,refLengthDic= ref_dic_realdata(refFile,MAX_LENGTH)
        refDic,refRead,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth = sam_parse(inFile,"./",nm,refDic,refRead,sr,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth)
        calculate_depth_coverage(key,refDic,countRead,refRead,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth,refLengthDic)
        Ans.append(ansDic[name]) 
    else:
        print (f"[Error] \"{key}\" has no Answer")
        return


    data ={"Sample":Sample,
        "Depth1":Depth2,
        "Depth2":Depth1,
        "Depth3":Depth3,
        "Answer":Ans}

    df = pd.DataFrame(data)

    return df

def makePrediction(dataset2,modelweight,outFile,depth):

    rgb, ans = extend_data(dataset2)

    o = open(outFile, 'w')
#    r = open(f"{args.outPre}_report.txt", 'w')
    
    model = create_cnn(rgb.shape[1], 3, regress=False)
    model.compile(loss="binary_crossentropy", optimizer= "adam", metrics = ["accuracy"])
    model.load_weights(modelweight)

    test_ans = model.predict(rgb)

    o.write("Sample,Answer,Predicted,AvgDepth\n")
    for i in range(len(test_ans)):
        avgDepth =[]
        for dd in dataset2.iloc[i]["Depth2"]:
            if dd >0:
                avgDepth.append(dd)
        if len(avgDepth)>0:
            avgDepth = statistics.median(avgDepth)
        else:
            avgDepth=0
        if avgDepth  >= depth:
            o.write(dataset2.iloc[i]["Sample"] +','+str(dataset2.iloc[i]["Answer"])+','+str(convert_to_0(test_ans[i])[0])+','+str(round(avgDepth,2))+'\n')
        elif 1 < avgDepth < depth and convert_to_0(test_ans[i])[0] ==1:
            o.write(dataset2.iloc[i]["Sample"] +','+str(dataset2.iloc[i]["Answer"])+',NA,'+str(round(avgDepth,2))+'\n')
        else:
            o.write(dataset2.iloc[i]["Sample"] +','+str(dataset2.iloc[i]["Answer"])+',0,'+str(round(avgDepth,2))+'\n')


    o.close()


if __name__ == "__main__":
    args = init_argv()
    args.nm = 0.9
    args.sr = 0.5
    args.depth = 5


    mylog = makeLog(f"{args.logFile}")
    mylog.info("Current Directory : "+subprocess.check_output("pwd",shell=True).decode('utf-8'))
    print("[DeepMobilome data transformation Start]")
    data=dataTransformation(args.inFile, args.ref, args.count, args.ans, args.nm, args.sr)

    print("[DeepMobilome data transformation End]")

    if not data.empty:

        print("[DeepMobilome prediction Start]")
        makePrediction(data, args.save, args.outFile, args.depth)

        print("[DeepMobilome prediction End]")


            

