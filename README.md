# MBMethPred
Medulloblastoma Subgroups Prediction from Methylation Data using AI 

## Install MBMethPred package
```{r}
remotes::install_github("sharifrahmanie/MBMethPred")
```
```{r}
require(MBMethPred)
```

## Input file for prediction

The `ReadMethylFile` is a function for reading DNA methylation files and use them as new data for prediction by every model. The input for this function should be either CSV or TSV file format.
 

### Usage
```{r}
set.seed(1234)
fac <- ncol(Data1)
NewData <- sample(data.frame(t(Data1[,-fac])),10)
NewData <- cbind(rownames(NewData), NewData)
colnames(NewData)[1] <- "ID"
write.csv(NewData, "NewData.csv", quote = FALSE, row.names = FALSE)
methyl <- ReadMethylFile(File = "NewData.csv")
file.remove("NewData.csv")
```

This function has only one argument, the File. While the first column is CpG methylation probs, starting with cg and followed by a number, other columns are samples with methylation values. All columns should be named.

## Box plot

The `BoxPlot` function draws a box plot out of the DNA methylation dataset or other data frames.

### Usage

```{r}

data <- Data2[1:20,]
data <- cbind(rownames(data), data)
colnames(data)[1] <- "ID"
BoxPlot(File = data, Projname = NULL)
```

This function has two arguments as follow:

* `File`    A data frame with the first column as ID. 
* `Projname` A string to name the plot.

## t-SNE 3D plot
The `TSNEPlot` function draws a 3D t-SNE plot for DNA methylation dataset using the K-means clustering technique. This function has two arguments `File` (any matrices) and `NCluster` ( number of clusters for K-Means clustering). 

### Usage 

```{r}
data <- data.frame(t(Data2[1:100,]))
data <- cbind(rownames(data), data)
colnames(data)[1] <- "ID"
TSNEPlot(File = data, NCluster = 4)
```

An R window will appear with a 3D projection of the t-SNE result. The plot object can be saved with the next line of code (uncomment).

```{r}
# rgl.snapshot('tsne3d.png', fmt = 'png')
```


## Input file for similarity network fusion (SNF)

Using `ReadSNFData` function, one can read files (any matrices with CSV or TSV format) and feed them into the similarity network fusion (SNF) function (from the SNFtools package).

### Usage

```{r}
data(Data2) # Gene expression 
Data2 <- cbind(rownames(Data2), Data2)
colnames(Data2)[1] <- "ID"
write.csv(Data2, "Data2.csv", row.names = FALSE)
Data2 <- ReadSNFData(File = "Data2.csv")
file.remove("Data2.csv")
```


## Similarity network fusion (SNF)

The `SimilarityNetworkFusion` is a function to perform SNF function (from SNFtool package) and output clusters.

### Usage 

```{r}
data(RLabels) # Real labels
data(Data2) # Methylation
data(Data3) # Gene expression
snf <- SimilarityNetworkFusion(Files = list(Data2, Data3),
                               NNeighbors  = 13,
                               Sigma = 0.75,
                               NClusters = 4,
                               CLabels = c("Group4", "SHH", "WNT", "Group3"),
                               RLabels = RLabels,
                               Niterations = 60)
snf
```

This function has several arguments as follow:

* `Files`   A list of data frames created using the ReadSNFData function or matrices.
* `NNeighbors`    The number of nearest neighbors.
* `Sigma`     The variance for local model.
* `NClusters`   The number of clusters.
* `CLabels`   A string vector to name the clusters. Optional.
* `RLabels`   The actual label of samples to calculate the Normalized Mutual Information (NMI) score. Optional.
* `Niterations`   The number of iterations for the diffusion process.

## Support vector machine model

The `SupportVectorMachineModel` is a function to train a support vector machine model to classify medulloblastoma subgroups using the DNA methylation dataset (Illumina Infinium HumanMethylation450). Prediction is followed by training if new data is provided.

Model metrics, including accuracy, precision, sensitivity F1-Score, specificity, and AUC_average can be calculated for the test dataset using the `ModelMetrics` function, which calculates the average of the above parameters from the result of the `ConfusionMatrix` function. 

The prediction result on new data can be accessed through the `NewDataPredictionResult` function, which calculates every prediction's mode across the number of cross-validation folds. 

### Usage 

```{r}
set.seed(1234)
fac <- ncol(Data1)
NewData <- sample(data.frame(t(Data1[,-fac])),10)
NewData <- cbind(rownames(NewData), NewData)
colnames(NewData)[1] <- "ID"

svm <- SupportVectorMachineModel(SplitRatio = 0.8, 
                                 CV = 10, 
                                 NCores = 1, 
                                 NewData = NewData)
ModelMetrics(Model = svm)
NewDataPredictionResult(Model = svm)
```

This function has the following arguments:

* `SplitRatio`    Train and test split ratio. A value greater or equal to zero and less than one.
* `CV`    The number of folds for cross-validation. It should be greater than one.
* `NCores`    The number of cores for parallel computing.
*  `NewData`    A methylation data for prediction. 


## K nearest neighbor model
The `KNearestNeighborModel` is a function to train a K nearest neighbor model to classify medulloblastoma subgroups using the DNA methylation dataset.

### Usage 

```{r}
set.seed(1234)
fac <- ncol(Data1)
NewData <- sample(data.frame(t(Data1[,-fac])),10)
NewData <- cbind(rownames(NewData), NewData)
colnames(NewData)[1] <- "ID"

knn <- KNearestNeighborModel(SplitRatio = 0.8, 
                             CV = 10, 
                             K = 3, 
                             NCores = 1, 
                             NewData = NewData)
ModelMetrics(Model = knn)
NewDataPredictionResult(Model = knn)
```

This function has the following arguments:

* `SplitRatio`    Train and test split ratio. A value greater or equal to zero and less than one.
* `CV`    The number of folds for cross-validation. It should be greater than one.
* `K`   The number of nearest neighbors.
* `NCores`    The number of cores for parallel computing.
*  `NewData`    A methylation data for prediction. 

## Random forest model
The `RandomForestModel` is a function to train a random forest model to classify medulloblastoma subgroups using the DNA methylation dataset.

### Usage

```{r}
set.seed(1234)
fac <- ncol(Data1)
NewData <- sample(data.frame(t(Data1[,-fac])),10)
NewData <- cbind(rownames(NewData), NewData)
colnames(NewData)[1] <- "ID"

rf <- RandomForestModel(SplitRatio = 0.8, 
                        CV = 10, 
                        NTree = 100, 
                        NCores = 1, 
                        NewData = NewData)
ModelMetrics(Model = rf)
NewDataPredictionResult(Model = rf)
```



This function has the following arguments:

* `SplitRatio`    Train and test split ratio. A value greater or equal to zero and less than one.
* `CV`    The number of folds for cross-validation. It should be greater than one.
* `NTree`   The number of trees to be grown.
* `NCores`    The number of cores for parallel computing.
*  `NewData`  A methylation data for prediction.  

## XGBoost model
The `XGBoostModel` is a A function to train an XGBoost model to classify medulloblastoma subgroups using the DNA methylation dataset.

### Usage

```{r}
set.seed(1234)
fac <- ncol(Data1)
NewData <- sample(data.frame(t(Data1[,-fac])),10)
NewData <- cbind(rownames(NewData), NewData)
colnames(NewData)[1] <- "ID"

xgboost <- XGBoostModel(SplitRatio = 0.8, 
                        CV = 10, 
                        NCores = 1, 
                        NewData = NewData)
ModelMetrics(Model = xgboost)
NewDataPredictionResult(Model = xgboost)
```

This function has the following arguments:

* `SplitRatio`    Train and test split ratio. A value greater or equal to zero and less than one.
* `CV`    The number of folds for cross-validation. It should be greater than one.
* `NCores`    The number of cores for parallel computing.
*  `NewData`  A methylation data for prediction.  


## Linear discriminant analysis model
The `LinearDiscriminantAnalysisModel` is a function to train a linear discriminant analysis model to classify medulloblastoma subgroups using the DNA methylation dataset.

### Usage

```{r}
set.seed(1234)
fac <- ncol(Data1)
NewData <- sample(data.frame(t(Data1[,-fac])),10)
NewData <- cbind(rownames(NewData), NewData)
colnames(NewData)[1] <- "ID"

lda <- LinearDiscriminantAnalysisModel(SplitRatio = 0.8, 
                                       CV = 10, 
                                       NCores = 1, 
                                       NewData = NewData)
ModelMetrics(Model = lda)
NewDataPredictionResult(Model = lda)
```

This function has the following arguments:

* `SplitRatio`    Train and test split ratio. A value greater or equal to zero and less than one.
* `CV`    The number of folds for cross-validation. It should be greater than one.
* `NCores`    The number of cores for parallel computing.
*  `NewData`  A methylation data for prediction.  

## Naive bayes model
The `NaiveBayesModel` is a function to train a Naive Bayes model to classify medulloblastoma subgroups using the DNA methylation dataset.

### Usage

```{r}
set.seed(1234)
fac <- ncol(Data1)
NewData <- sample(data.frame(t(Data1[,-fac])),10)
NewData <- cbind(rownames(NewData), NewData)
colnames(NewData)[1] <- "ID"

nb <- NaiveBayesModel(SplitRatio = 0.8, 
                      CV = 10, 
                      Threshold = 0.8, 
                      NCores = 1, 
                      NewData = NewData)
ModelMetrics(Model = nb)
NewDataPredictionResult(Model = nb)
```

This function has the following arguments:

* `SplitRatio`    Train and test split ratio. A value greater or equal to zero and less than one.
* `CV`    The number of folds for cross-validation. It should be greater than one.
* `Threshold`   The threshold for deciding class probability. A value greater or equal to zero and less than one.
* `NCores`    The number of cores for parallel computing.
*  `NewData`  A methylation data for prediction.  

## Artificial neural network model

The `NeuralNetworkModel` is a function to train an artificial neural network model to classify medulloblastoma subgroups using the DNA methylation dataset. Please uncomment the following lines and run the function. If it is the first time you run this function, set the InstallTensorFlow parameter to TRUE. It will automatically install the Python and TensorFlow library (version 2.10-cpu) in a virtual environment. Then set the parameter to FALSE. 

### Usage

```{r}
set.seed(1234)
fac <- ncol(Data1)
NewData <- sample(data.frame(t(Data1[,-fac])),10)
NewData <- cbind(rownames(NewData), NewData)
colnames(NewData)[1] <- "ID"
ann <- NeuralNetworkModel(Epochs = 100,
                          NewData = NewData,
                          InstallTensorFlow = TRUE)
ModelMetrics(Model = ann)
NewDataPredictionResult(Model = ann)
```
This function has the following arguments:

* `Epochs`    The number of epochs.
* `NewData`   A methylation data from ReadMethylFile function.
* `InstallTensorFlow`   Logical. Running this function for the first time, you need to install TensorFlow library (V 2.10-cpu). Default is TRUE.
