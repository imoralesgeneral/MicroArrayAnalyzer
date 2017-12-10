MicroArrayAnalyzer
====
Web-tool developed to analyze microarrays
----
***Thesis. Master degree in Bioinformatics.***

  Microarray analysis to calculate the differential expression of genes is one of the techniques in which bioinformatics plays an important role. However, these analyzes requires that user has some programming expertize. 
  
  This work aims to create a web application that allows users who do not have enough programming skills to perform the differential expression analysis of genes from microarrays. To do this, a pipeline is developed in the R programming language that will mark the path to be followed for analysis and this pipeline will be implemented in a web tool created using Shiny. 
  
  To verify the correct functionality of both the pipeline and the web application, three datasets are downloaded to allow the analysis to be performed using different parameters. The first dataset will need a Mus musculus annotation file, while the other two datasets will use Homo sapiens. In addition, dataset number three will contain three groups to make comparisons between them.
  
  After the implementation of the tool, it shows a correct operation for the programmed functions. However, it would be advisable to extend the tool by granting a larger number of predetermined annotation files, different types of data normalization, and RNA-seq and NGS analysis.

##### Datasets used:
1. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27174]. Platform: Affymetrix Mouse Gene 1.0 ST Array.

2. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62947]. Platform: Affymetrix Human Genome U133 Plus 2.0 Array.

3. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42771]. Platform: Affymetrix Human Genome U133 Plus 2.0 Array.
