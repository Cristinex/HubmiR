# SiCmiR

Predicting miRNA Expression Level in Single-Cell Clusters for hub-miRNA Discovery

## Usage

## Dependency

The following packages are required for running SiCmiR:

```
R

Seurat
```

```
python
	
sys, torch, torch.nn, numpy, pandas, sklearn.model_selection
```

## Command

- Step1: Process the input data `location/for/Rscript /location/for/package/SiCmiR/script/input_processing.R args[1]  args[2] args[3] `.
		
    - ALL Parameters should be quoted by "" or ''.
	
    - parameter `args[1]` select which type of pooling of data:
	  "celltypeavg" for average of annotated cell type.
      "pooledavg" for random pooling of cells in the same cell type.
	- parameter `args[2]`  file directory for input profile. Format see the `Demo_input.csv` in files folder. 
	  This directory should contain ONLY the input files. End with "/".
	- parameter `args[3]` output directory. End with "/".
	


- Step2: Input the processed data into SiCmiR for miRNA expression profile prediction `location/for/python /location/for/package/SiCmiR/script/DNN.py args[4] `.
		
	- ALL Parameters should be quoted by "" or ''.
	- parameter `args[4]` Same directory as `args[3]` where the processed file(s) was(were) stored. End with "/".
