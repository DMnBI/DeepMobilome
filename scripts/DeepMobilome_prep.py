import sys, os, time, argparse,subprocess, re, math,statistics
import logging
from collections import defaultdict
import pandas as pd
import pickle

from DeepMobilome_common import makeLog
from DeepMobilome_common import ref_dic_realdata
from DeepMobilome_common import sam_parse
from DeepMobilome_common import save_read

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
                    Get fragment coverage by SAM file
                    Refine == only when both reads meet the criteria
                    Only paired reads with -> <- structures are used"""
                     
                ,formatter_class= argparse.RawTextHelpFormatter)
    parser.add_argument("-i","--inFile",dest="inFile",metavar="FILE",required=True, help="sam File list")
    parser.add_argument("-in","--inDir",dest="inDir",metavar="DIR", required=True,help="sam File directory")
    parser.add_argument("-o","--outPre",dest="outPre",metavar="FILE",required=True,help="output File prefix {outPre}.log, {outPre}.plk, {outPre}list.txt")
    parser.add_argument("-r","--ref",dest="ref",metavar="FILE",required=True, help="bowtie reference file")
    parser.add_argument("-nm","--nm",dest="nm",default = 0.9, help = "mapped read nm ratio default : 0.9", type =float)
    parser.add_argument("-s","--sr",dest = "sr",default = 0.5, help = "soft clipped ratio defalut : 0.5", type = float)
    parser.add_argument("-c","--count",dest = "count",required=True, help ="read count list",)
    parser.add_argument("-a","--ans",dest= "ans",metavar ="FILE", help = "Answer of the given case, samplename and answer, csv file",required=True)
    args = parser.parse_args()
    return args

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

if __name__ == "__main__":
    parser = init_argv()
    f = open(parser.inFile)
    fList = f.readlines()
    countRead = save_read(parser.count)
    ansDic = read_ans(parser.ans)

    mylog = makeLog(f"{parser.outPre}.log")
    mylog.info("Current Directory : "+subprocess.check_output("pwd",shell=True).decode('utf-8'))
    print("[DeepMobilome_prep Start]")
    print ("Used Command : python "+' '.join(sys.argv))

    for num in range(0,len(fList)):
        key = fList[num].split('.sam')[0]
        name = key.split('_')[0]
        if name in ansDic:
            refDic,refRead,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth,refLengthDic= ref_dic_realdata(parser.ref,MAX_LENGTH)
            refDic,refRead,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth = sam_parse(fList[num].rstrip('\r\n'),parser.inDir,parser.nm,refDic,refRead,parser.sr,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth)
            calculate_depth_coverage(key,refDic,countRead,refRead,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth,refLengthDic)
            Ans.append(ansDic[name]) 
        else:
            print (f"[Error] \"{key}\" has no Answer")

    data ={"Sample":Sample,
        "Depth1":Depth2,
        "Depth2":Depth1,
        "Depth3":Depth3,
        "Answer":Ans}

    df = pd.DataFrame(data)
    with open(f"{parser.outPre}.plk",'wb') as f:
        pickle.dump(df,f)
    mylog.info("[DeepMobilome_prep Complete]")
    f.close()
