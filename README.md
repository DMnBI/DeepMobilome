## DeepMobilome: predicting mobile genetic elements using sequencing reads of microbiomes

This repository includes the implementation of ' DeepMobilome: predicting mobile genetic elements using sequencing reads of microbiomes '. Please cite our paper if you use our models. Fill free to report any issue for maintenance of our model.

## Citation
If you have used DeepMobilome in your research, please cite the following publication:
Youna Cho, Mina Rho

## Overview
Here is the brief overview of DeepMobilome’s prediction workflow.

  <img width="451" alt="image" src="https://github.com/user-attachments/assets/9d01c2ee-2c85-46cc-b03e-13352dcf3ec8" />

DeepMobilome predicts target MGEs through three phases: input generation, data transformation, and model prediction (Fig. 1). a) In the input generation phase, sample reads are aligned to the target MGE sequences using a read aligner such as Bowtie2. b) The read alignment information is processed in the data transformation phase. c) In the model prediction phase, DeepMobilome predicts the presence of the target MGEs in the sample based on the representation learned from the read mapping information.

## 1.Setup
We strongly recommend you to use python virtual environment with Anaconda/Miniconda. This model can run on GPU/TPU machine, but this is not mandatory. The details of the machine/environment we used are as follows:
+ NVIDIA Tesla V100 DGXS 32 GB graphic memory
+ CUDA 11.0
+ python 3.9.12
+ tensorflow=2.4.1

### 1.1 Build environment
```
git clone https://(gitclone site)
cd environment
conda env create -f DeepMobilome.yml
conda activate DeepMobilome
export PATH=/path/to/DeepMobilome/Script/:$PATH
```

## 2. Predict 

## 2.1 Running DeepMobilome

you can find details of required/optional parameters for each function with -h option.
```
deepmobilome -h
```
```
usage: deepmobilome [-h] -i FILE -o FILE -l FILE -r FILE -c COUNT -a FILE -s SAVE

                    Predict Target

                    Input: Model weight, SAM file
                    Output: Predicted result

                    

optional arguments:
  -h, --help            show this help message and exit
  -i FILE, --inFile FILE
                        sam File : sample reads aligned to target sequence
  -o FILE, --outFile FILE
                        predicted result / samplename,Answer,Predict,AvgDepth / csv file
  -l FILE, --logFile FILE
                        log file
  -r FILE, --ref FILE   target sequence file / fna file
  -c COUNT, --count COUNT
                        read count / samplename,read count / csv file
  -a FILE, --ans FILE   answer for the given case / samplename,answer / csv file
  -s SAVE, --save SAVE  model weight
```
If you don't know the answer, just put either 0 or 1 for every given sample.


### 2.2 Tutorial
Predict existence of the target sequence in sample using the trained model.
Sample1 and sample2 are human data from [Cho, Y., et al., Deciphering Resistome in Patients With Chronic Obstructive Pulmonary Diseases and Clostridioides difficile Infections. Front Microbiol, 2022. 13: p. 919907.]
Target sequence is 'IS26/AAC(6')-Ie-APH(2'’)-Ia/N-acetyltransferase/IS26'
fastq file of sample1 and sample2 is uploaded at xxx(dropbox link)

```
cd tutorial/

```

### 2.2.1 Positive sample
Predict positive sample as positive

```
deepmobilome \
--inFile positive_sample/sample1_target1.sam \
-o output/sample1Target1.csv -l output/sample1Target1.log \
--ref target/Target1seq.fna --count positive_sample/sample1_read_detail.txt \
--ans positive_sample/sample1target1answer.csv --save ../model/DeepMobilome_weight.h5
```

```
[Prediction example]
Sample,Answer,Predicted,AvgDepth
sample1_target1_Target1,1,1,148.44
```


### 2.2.2 Negative sample
Predict negative sample as negative

```
deepmobilome \
--inFile negative_sample/sample2_target1.sam \
-o output/sample2Target1.csv -l output/sample2Target1.log \
--ref target/Target1seq.fna --count negative_sample/sample2_read_detail.txt \
--ans negative_sample/sample2target1answer.csv --save ../model/DeepMobilome_weight.h5
```

```
[Prediction example]
Sample,Answer,Predicted,AvgDepth
sample2_target1_Target1,0,0,161.63
```
