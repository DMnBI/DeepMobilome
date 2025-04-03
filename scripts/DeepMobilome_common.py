import sys, os, time, argparse,subprocess, re, math,statistics
import logging
from collections import defaultdict
import pandas as pd
import pickle

def makeLog(logFile):
    logging.basicConfig(
            level= logging.DEBUG,
            format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
            filename=logFile,
            filemode ='w'
    )
    stdout_logger = logging.getLogger('STDOUT')
    sl = StreamToLogger(stdout_logger,logging.INFO)
    sys.stdout = sl
    stderr_logger = logging.getLogger('STDERR')
    sl = StreamToLogger(stderr_logger, logging.ERROR)
    sys.stderr = sl

    return logging                   


class StreamToLogger(object):
    def __init__(self, logger, log_level=logging.INFO):
        self.logger= logger
        self.log_level = log_level
        self.linebuf = ''
    def write(self,buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())
    def flush(self):
        pass

def accumulate(num_list):
    sum_list = []
    for num,value in enumerate(num_list):
        if num!=0:
            sum_list.append(value+sum_list[num-1])
        else:
            sum_list.append(value)
    return sum_list

def sam_parse(key,inDir,nm,refDic,refRead,sr,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth):
    """Parse SAM file and save data to the dictionary 

    Args:
        key (str): file key name
        nm (float): NM minimum similarity (default 0.9)
        refDic (dict)
        refRead (dict)            
        sr (float): softclip mimimum ration (default 0.5)
        refDNA (dict)
        sampleDNA (dict)
        linkmatrix (dict)           
        readposcount (dict)
        insertdepth (dict)
    
    Returns:
        refDic
        refRead
        refDNA
        sampleDNA
        linkmatrix
        readposcount
        insertdepth
    """    
    open_sam = open(inDir+'/'+key)
    all_sam = {}

    hasfus=0
    for line in open_sam:
        readfus=0
        line = line.rstrip('\r\n')
        column =line.split('\t')
        #TRUE :Not header, yes alignment tags, proper SAM file
        if not line.startswith('@') and len(column)>=11 and column[7].find(':')<0 and column[1].find('_')<0:
            if ".1" in column[0][-2:] or ".2" in column[0][-2:]:
                query = column[0][:-2]
            else:
                query  = column[0]

            flag   = int(column[1])
            CIGAR  = column[5]
            qual = column[10:]
            refPos = int(column[3])
            refName = column[2]
            seq = column[9]
            rFlag = column[6]
            rrefpos = int(column[7])
            #ipos 1st read position, apos 2nd read position 
            ipos =min(refPos,rrefpos)
            apos = max(refPos,rrefpos)
            isize = abs(int(column[8]))
            reverseflag =0 
            properflag = 0
            
            if flag & 0x10 ==0:
                reverseflag +=1
            if flag &0x20 ==0:
                reverseflag -=1
            if flag &0x2 !=0:
                properflag =1
            #TRUE: mapped, mate mapped to the same referense, insertion size limit, -> <- format     
            if flag & 0x4 == 0 and '=' in rFlag and 50<=isize<=600  and reverseflag !=0 and properflag ==1:
                charCIGAR = re.findall('\D+',CIGAR)
                numCIGAR = re.findall('\d+',CIGAR)
                accum_num = accumulate(numCIGAR)
                #single sortclip allowed 
                if CIGAR.count('M')>=1 and CIGAR.count('S')<=1:
                    mPos = getIndexPosition(charCIGAR,'M')
                    mLength =0 
                    for m in mPos:
                        mLength += int(numCIGAR[m])
                    readLength = len(seq)
                    if readLength * float(sr)< mLength:
                        for i in qual:
                            if 'NM' in i:
                                NM_num = int(i.split(':')[2])
                                if NM_num<=(1-float(nm))*mLength:
                                    if not query in all_sam:
                                        all_sam[query] = []
                                    all_sam[query].append(line)
    readnum =0
    #read pair dictionary k: read name v: pair samfile line   
    for k, v in all_sam.items():
        readnum+=1
        pair = len(v)
        pflag =0
        nextflag =0 
        sflag =0 
        linkflag =0
        if pair ==2:
            l1 = v[0].split('\t')[2]
            l2 = v[1].split('\t')[2]
            p1 = int(v[0].split('\t')[3])
            p2 = int(v[1].split('\t')[3])
            pr1 = len(v[0].split('\t')[9])
            pr2 = len(v[1].split('\t')[9])
            c1 = v[0].split('\t')[5]
            cc1 = re.findall('\D+',c1)
            cd1 = re.findall('\d+',c1)
            c2 = v[1].split('\t')[5]
            cc2 = re.findall('\D+',c2)
            cd2 = re.findall('\d+',c2)
            l1isize = abs(int(v[0].split('\t')[8]))
            if p1<p2:
                s1,s2 = p1,p2
                sr1,sr2 = pr1,pr2
            else:
                s1,s2 = p2,p1
                sr1,sr2 = pr2,pr1
            if l1 ==l2:
                pflag =1
            #only one 'M' in CIGAR
            if len(cc1) ==1  and len(cc2) ==1 and pflag ==1:
               if 50<=l1isize<=600 :
                   linkflag =1

        if linkflag ==1:
            for linkpos in range(s1,s2+sr2):
                linkmatrix[l1][linkpos]+=1
                readposcount[l1][linkpos].add(readnum)
            for linkpos in range(s1+sr1,s2):
                insertdepth[l1][linkpos]+=1
        if pflag ==1:
            for nflag,line in enumerate(v):
                column = line.split('\t')
                query  = column[0]
                flag   = int(column[1])
                CIGAR  = column[5]
                qual = column[10:]
                refPos = int(column[3])
                refName = column[2]
                seq = column[9]
                rFlag = column[6]
                rrefpos = int(column[7])
                ipos =min(refPos,rrefpos)
                apos = max(refPos,rrefpos)
                if flag & 0x4 == 0 and '=' in rFlag:
                    charCIGAR = re.findall('\D+',CIGAR)
                    numCIGAR = re.findall('\d+',CIGAR)
                    accum_num = accumulate(numCIGAR)
                    if CIGAR.count('M')>=1:
                        mPos = getIndexPosition(charCIGAR,'M')
                        mLength =0 
                        for m in mPos:
                            mLength += int(numCIGAR[m])
                        readLength = len(seq)
                        if readLength * float(sr)< mLength:
                            for i in qual:
                                if 'NM' in i:
                                    NM_num = int(i.split(':')[2])
                                    if NM_num<(1-float(nm))*mLength:
                                        refRead[refName]+=1
                                        delNum =0 
                                        currentpos =0
                                        for index,cigar in enumerate(charCIGAR):
                                            if cigar == 'M' or cigar=='D':
                                                for num in range(0,int(numCIGAR[index])):
                                                    if num+refPos+delNum+currentpos in refDic[refName]:
                                                        refDic[refName][currentpos+num+refPos+delNum]+=1
                                                    else:
                                                        print (refDic[refName])
                                                        print ("Wrong Position" + line +','+str(num+refPos+delNum+currentpos)+'\n')
                                                                    
                                                currentpos += int(numCIGAR[index])+1
                                            elif cigar =='I':
                                                delNum -=int(numCIGAR[index])
                                                 

                                    break
           
    open_sam.close()
    return refDic,refRead,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth


def ref_dic_realdata(refFile,MAX_LENGTH):
    """ Make the baseline reference dictionary by reference length and sequence - Read data
    
    Args: 
        refFile (str): The reference file (.fna/.fa) name 
        MAX_LENGTH (int): maximum number of pattern length

    Returns: 
        All empty dictionary initialized by reference length

        refDic (dict): save read depth per basepair
        refRead (dict): save read count
        refDNA (dict): save read sequence
        sampleDNA (dict): save read sequence by ACGT count
        linkmatrix (dict): save read depth per basepair
        readposcount (dict): save read name
        insertdepth (dict): save read insert size depth                      
    """


    mol = ['A','C','G','T','N']
    r = open(refFile)
    refDic ={} # save only matched points
    refRead ={}
    refDNA = {}
    sampleDNA ={}
    linkmatrix = {}# save only match read
    readposcount = {}# save read number to find steep position
    insertdepth={} #count insertion depth 
    refLengthDic = {}

    total =0
    name = ""
    seq = ""
    for line in r:
        line = line.rstrip('\r\n') 
        if line[0]==">" and name =="":
            name = line[1:].split(' ')[0]
            refDic[name]={}
            linkmatrix[name]={}
            insertdepth[name]={}
            readposcount[name]={}
            refRead[name]=0
            refDNA[name] = {}
            sampleDNA[name]={}
            refLengthDic[name] = 0
        elif line[0]==">" and seq!="":
            refLengthDic[name]=total
            for i in range(1,total+1):
                refDic[name][i]=0
                linkmatrix[name][i]=0
                insertdepth[name][i]=0
                readposcount[name][i]=set()

                refDNA[name][i]=seq[i-1]
                sampleDNA[name][i]={}
                for j in mol:
                    sampleDNA[name][i][j]=0
                sampleDNA[name][i][seq[i-1]]+=1
            if total < MAX_LENGTH:
                for i in range(total,MAX_LENGTH+1):
                    refDic[name][i]=-5
                    linkmatrix[name][i]=-5
                    insertdepth[name][i]=-5
                    readposcount[name][i]=set()

                    refDNA[name][i]='N'
                    sampleDNA[name][i]={}
                    for j in mol:
                        sampleDNA[name][i][j]=0
                    sampleDNA[name][i]['N']+=1

            name =line[1:].split(' ')[0]
            refDic[name]={}
            linkmatrix[name]={}
            insertdepth[name]={}
            readposcount[name]={}
            refLengthDic[name]=0

            refRead[name]=0
            refDNA[name]={}
            sampleDNA[name] = {}
            total =0
            seq = ""
        else :
            total +=len(line)
            seq += line
    if seq!="":
        refLengthDic[name]=total
        for i in range(1,total+1):
            refDic[name][i]=0
            linkmatrix[name][i]=0
            insertdepth[name][i]=0
            readposcount[name][i]=set()
            refDNA[name][i]=seq[i-1]
            sampleDNA[name][i]={}
            for j in mol:
                sampleDNA[name][i][j]=0
            sampleDNA[name][i][seq[i-1]]+=1
        if total < MAX_LENGTH:
            for i in range(total,MAX_LENGTH+1):
                refDic[name][i]=-5
                linkmatrix[name][i]=-5
                insertdepth[name][i]=-5
                readposcount[name][i]=set()

                refDNA[name][i]='N'
                sampleDNA[name][i]={}
                for j in mol:
                    sampleDNA[name][i][j]=0
                sampleDNA[name][i]['N']+=1


    return refDic,refRead,refDNA,sampleDNA,linkmatrix,readposcount,insertdepth,refLengthDic


def getIndexPosition(listOfElements, element):
    indexPosList = []
    indexPos = 0
    while True:
        try :
            indexPos = listOfElements.index(element, indexPos)
            indexPosList.append(indexPos)
            indexPos +=1 
        except ValueError as e:
            break
    return indexPosList


def save_read(countFile):
    tmp = {}
    num =0
    with open(countFile,'r') as f:
        for line in f:
            if num!=0:
                line = line.rstrip('\r\n')
                items =line.split(',')
                if not "Read" in line:
                    if 'fastq' in items[0]:
                        tmp[items[0].split('.fastq')[0].split('_')[0]] = str(float(items[1]))
                    elif 'fq' in items[0]:
                        
                        tmp[items[0].split('.fq')[0].split('_')[0]] = str(float(items[1]))

            num+=1
    return tmp

