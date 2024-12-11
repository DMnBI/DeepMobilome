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
```bash
git clone https://(gitclone site)
cd env
conda env create -f DeepMobilome.yml
conda activate DeepMobilome
export PATH=/path/to/DeepMobilome/Script/:$PATH
```

## 2. Predict 
Predict existence of the target sequence in sample using the trained model.
Sample1 and sample2 are human data from 

### 2.1 Positive sample


