#
# This is a Shiny web application. You can run the application by clicking.
# 
# This application allow you analyse microarray data.
# 
# Author: Isaac Morales General
# Version: 1.0
#

# Required packages

library(shiny)
library(shinydashboard)
library(oligo)
library(annotate)
library(ggplot2)
library(limma)
library(GO.db)

# Max size Upload
options(shiny.maxRequestSize=50*1024^2) 

ui <- dashboardPage(
		
		dashboardHeader(title = "MicroArray AnalyZer", titleWidth = 220), 
		skin = "red",
		
		dashboardSidebar(
				width = 220,
				sidebarMenu(
						br(),
						menuItem("Collecting Data", tabName = "collecting", icon = icon("file")),
						menuItem("Quality Control", tabName = "quality", icon = icon("thumbs-up")),
						menuItem("Selection Settings", tabName = "adj", icon = icon("gears"),
								selectInput("annotation", "Default annotation file", choices = list("I want insert my .annot file"="other", "mogene10sttranscriptcluster"="mogene10sttranscriptcluster.db", "hgu133a"="hgu133a.db",
												"hugene21sttranscriptcluster"="hugene21sttranscriptcluster.db", "hugene20sttranscriptcluster" = "hugene20sttranscriptcluster.db",
												"clariomdhumantranscriptcluster" = "clariomdhumantranscriptcluster.db", "clariomshumanhttranscriptcluster"="clariomshumanhttranscriptcluster.db",
												"clariomshumantranscriptcluster"="clariomshumantranscriptcluster.db", "clariomsmousehttranscriptcluster"="clariomsmousehttranscriptcluster.db",
												"clariomsmousetranscriptcluster"="clariomsmousetranscriptcluster.db"), 
										selected = NULL, multiple = FALSE,selectize = FALSE, width = NULL, size = NULL),
								fileInput("ann", "Annotation file", multiple = FALSE),
								radioButtons("groups", label="Groups to compare", inline = TRUE, 
										choices = list("Group 1 vs Group 2" = "1", "Group 1 vs Group 3" = "2", "Group 1 vs Group 4" = "4", 
												"Group 2 vs Group 3" = "3", "Group 2 vs Group 4" = "5", "Group 3 vs Group 4" = "6"), 
										selected = 1),
								numericInput("lfc", "Select LFC: ", value = 1, min = 0, max = 10, step = 1, width = NULL),
								selectInput("pval", "Adj. p-val", choices = list("0.001", "0.005", "0.01","0.05","0.1","0.5","1"), selected = "0.1", multiple = FALSE,
										selectize = TRUE, width = NULL, size = NULL),
								sliderInput("volcano", "Number of genes (Volcano plot):", min = 1, max = 50, value = 3),
								sliderInput("max_genes", "Number of genes (Heatmap):", min = 0, max = 500, value = 100)
						),
						menuItem("Results", tabName = "selection", icon = icon("search")),
						menuItem("Selection Settings GO", tabName = "adjGO", icon = icon("gears"),
								selectInput("pvalGO", "P. value", choices = list(0.001, 0.005, 0.01,0.05,0.1,0.5,1), selected = 0.01, multiple = FALSE,
										selectize = TRUE, width = NULL, size = NULL),
								radioButtons("ontol", label="Ontology", inline = TRUE, 
										choices = list("CC (Cellular component)" = "CC", "MF (Molecular function)" = "MF", "BP (Biological process)" = "BP"), 
										selected = "BP")
						),
						menuItem("GO Analysis", tabName = "ontology2", icon = icon("sitemap")),
						menuItem("Gene Annotation", tabName = "ontology", icon = icon("sitemap")),
						menuItem("Help", tabName = "help", icon = icon("question-circle"))
				)
		),
		
		dashboardBody(
				# tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "style.css")),
				tabItems(
						tabItem(tabName = "collecting",
								h1("Collecting Data"),
								h3("Please, select the path of the '.CEL' file and the target file"), 
								br(),
								fluidPage(splitLayout(
												textInput("CEL", "Path .CEL files", placeholder = "Path"),
												fileInput("target", "Target file", multiple = FALSE, accept = c('text/csv','text/comma-separated-values,text/plain','.csv'))
										)
								),
								fluidPage(splitLayout(
												box(title = "Files .CEL", status = "warning", solidHeader = TRUE, 
														tableOutput("celfiles"), height = 450, width = 1500), 
												box(title = "Target File", status = "warning", solidHeader = TRUE, 
														tableOutput("targetfile"), height = 450, width = 1500)
										)
								),
								fluidPage(splitLayout(
												h4("Number of groups: "),
												textOutput("gro"),
												h4("Groups: "),
												textOutput("grou")
										)
								)
						),
						tabItem(tabName = "quality",
								h1("Quality control"), 
								br(),
								fluidPage(splitLayout(
												box(title = "Boxplot before normalization", status = "warning", solidHeader = TRUE,
														plotOutput("box_bef"), height = 475, width = 600),
												box(title = "Boxplot after normalization", status = "warning", solidHeader = TRUE,
														plotOutput("box_aft"), height = 475, width = 600)
										)
								),
								fluidPage(splitLayout(
												box(title = "Dendogram before normalization", status = "warning", solidHeader = TRUE,
														plotOutput("den_bef"), height = 475, width = 600),
												box(title = "Dendogram after normalization", status = "warning", solidHeader = TRUE,
														plotOutput("den_aft"), height = 475, width = 600)
										)
								)
						),
						tabItem(tabName = "selection",
								h1("Results"), 
								br(),
								h4("Only the first items will be displayed. You can download the complete table."),
								fluidPage(splitLayout(
												box(title = "Top Table",status = "warning", solidHeader = TRUE, 
														tableOutput("selected1"), downloadButton("dwtable", "Download Genes"),  height = 450, width = 500),
												box(title = "Selected", status = "warning", solidHeader = TRUE,
														tableOutput("selected2"), downloadButton("dwtable2", "Download Genes"), height = 450, width = 500) 
										)
								),
								fluidPage(splitLayout(
												box(title = "Up-regulated Genes", status = "warning", solidHeader = TRUE, 
														tableOutput("upgenes"), downloadButton("dwup", "Download Up Genes"), height = 450, width = 500), 
												box(title = "Down-regulated Genes", status = "warning", solidHeader = TRUE, 
														tableOutput("downgenes"), downloadButton("dwdw", "Download Down Genes"), height = 450, width = 500)
										)
								),
								fluidPage(splitLayout(
												box(title = "Volcano Plot", status = "warning", solidHeader = TRUE, 
														plotOutput("volcan"), height = 475, width = 500), 
												box(title = "Heat Map", status = "warning", solidHeader = TRUE, 
														plotOutput("heatm"), height = 475, width = 500)
										
										
										)
								),
								textOutput("prueba")
						),
						tabItem(tabName = "ontology2",
								h1("Gene Ontology Analysis"),
								br(),
								h4("Only the first items will be displayed. You can download the complete table (File .html)"),
								box(title = "Ontology", status = "warning", solidHeader = TRUE, 
										#          tableOutput("go2"), downloadButton("dwgo2", "Download table"), height = 750, width = 1000)
										tableOutput("go2") , downloadButton("dwgo2", "Download table"), height = 750, width = 1000)
						),
						tabItem(tabName = "ontology",
								h1("Gene Annotation"),
								br(),
								h4("Only the first items will be displayed. You can download the complete table."),
								box(title = "Annotation", status = "warning", solidHeader = TRUE, 
										tableOutput("go"), downloadButton("dwgo", "Download table"), height = 750, width = 1000)
						),
						tabItem(tabName = "help",
								h2("What is this?"),
								p("The purpose of this application is to help people with little knowledge of bioinformatics to perform simple microarray analysis. Developed by Isaac Morales General. Microarray Analyzer v1.0."),
								h2("How can you use it?"),
								p("You need to perform this steps to analyse your data:"),
								tags$ol(
										tags$li(p(strong("Collecting data: ")," You need to indicate: "),
												tags$ul(tags$li(p(strong("Path .cel files: "), "Indicate the path in which you have storaged your 
																				.cel files"))), 
												tags$ul(tags$li(p(strong("Target file: "), "The targets is a file where you defines 
																				the levels of each variable under study. The application will use this file to
																				colour the samples in the plot. Upload a .csv file  with your variables
																				and their levels in columns. Please use the first column to put the name of the 
																				.cel file, the second column to put the name of the group, in the third column the name 
																				of the sample and in the last column the colour to 
																				identify them in the plot. Every column must be separated by commas. See" 
																		,tags$a(href="target.csv", target="_blank", "target.csv"),"to see 
																				and example of a target file.")))),
										tags$li(p(strong("Quality control: ")," In this section you can check the normalization and the grouping
																of data: "),
												tags$ul(tags$li(p(strong("Boxplots: "), "The boxplots will give an idea of the distribution of the intensities"))), 
												tags$ul(tags$li(p(strong("Dendograms: "), "Carry out a grouping of the samples by degree of similarity")))),
										tags$li(p(strong("Selection settings: "), "We need to configure some parameters before run the 
																application"),
												tags$ul(tags$li(p(strong("Select the annotation file:"), "You have two choices: Select one of the annotation files
																				listed in the dropdown menu or upload your own .annot file"))),
												tags$ul(tags$ul(tags$li(p("Predefined annotation files: "))),
														tags$ul(tags$ul(tags$li(tags$dd(p(strong("mogene10sttranscriptioncluster: "), "MOUSE GENE ARRAY 1.0 ST"))))),
														tags$ul(tags$ul(tags$li(tags$dd(p(strong("hgu133a: "), "HUMAN GENOME ARRAY 133A"))))),
														tags$ul(tags$ul(tags$li(tags$dd(p(strong("hgu21sttranscriptioncluster: "), "HUMAN GENE 2.1 ST ARRAY PLATE"))))),
														tags$ul(tags$ul(tags$li(tags$dd(p(strong("hgu20sttranscriptioncluster: "), "HUMAN GENE ARRAY 2.0 ST"))))),
														tags$ul(tags$ul(tags$li(tags$dd(p(strong("clariomdhumantranscriptcluster: "), "HUMAN CLARIOM D ARRAY"))))),
														tags$ul(tags$ul(tags$li(tags$dd(p(strong("clariomshumanhttranscriptcluster: "), "HUMAN CLARIOM S HT ARRAY"))))),
														tags$ul(tags$ul(tags$li(tags$dd(p(strong("clariomshumantranscriptcluster: "), "HUMAN CLARIOM S ARRAY"))))),
														tags$ul(tags$ul(tags$li(tags$dd(p(strong("clariomsmousehttranscriptcluster: "), "MOUSE CLARIOM S HT ARRAY"))))),
														tags$ul(tags$ul(tags$li(tags$dd(p(strong("clariomsmousetranscriptcluster: "), "MOUSE CLARIOM S ARRAY")))))),
												tags$ul(tags$ul(tags$li(p("If you select the .annot file, you must remove the head. See examples: ",tags$a(href="annot11.jpg", target="_blank", "Before remove"),
																				(" / ")    ,tags$a(href="annot21.jpg", target="_blank", "After remove"))))),
												tags$ul(tags$li(p(strong("Select the groups to compare:"), "Here you can choose
																				which groups you wish to compare"))),
												tags$ul(tags$li(p(strong("Select LFC:"), "Here you have to put the minimum 
																				absolute log2-fold-change required. (0-10)")))),
										tags$ul(tags$li(p(strong("Select Adj. p-val:"), "Select the cutoff value for adjusted p-values. Only genes with 
																		lower p-values are listed"))),
										tags$ul(tags$li(p(strong("Select the number of genes (Volcano Plot):"), "Select the number of genes showed in the 
																		Volcano Plot"))),
										tags$ul(tags$li(p(strong("Select the number of genes (Heatmap):"), "Select the number of genes showed in the 
																		Heatmap"))),
										tags$li(p(strong("Results: You can check the result of the analysis"),": "),
												tags$ul(tags$li(p(strong("Top Table: "), "Genes used for the analysis. Only the first items will be displayed. You can download the complete table."))), 
												tags$ul(tags$li(p(strong("Selected: "), "Genes that meet the specified requirements. Only the first items will be displayed. You can download the complete table."))), 
												tags$ul(tags$li(p(strong("Up-regulated Genes: "), "Genes that meet the specified requirements and are overexpressed. Only the first items will be displayed. You can download the complete table."))), 
												tags$ul(tags$li(p(strong("Down-regulated Genes: "), "Genes that meet the specified requirements and are under-expressed. Only the first items will be displayed. You can download the complete table."))), 
												tags$ul(tags$li(p(strong("Volcano Plot: "), "Will show the genes as a function of the expression changes and the p-value."))), 
												tags$ul(tags$li(p(strong("Heat Map: "), "Will display the expressions of each gene grouping them according to the sample.")))),
										tags$li(p(strong("Selection Settings GO: "),"You need to configure some parameter before run the GO analysis")),
										tags$ul(tags$li(p(strong("Select the P. value:"), "Select the cutoff value"))),
										tags$ul(tags$li(p(strong("Select the Ontolgy:"), "Select the Ontology which you wish analyse: "))),
										tags$ul(tags$ul(tags$li(p("CC: Cellular Component")))),
										tags$ul(tags$ul(tags$li(p("MF: Molecular Function")))),
										tags$ul(tags$ul(tags$li(p("BP: Biological Process")))),
										tags$li(p(strong("GO Analysis: "),"Which functions are differentially expressed with respect to the total of the analyzed genes.    ** Only works with predefined annotation files **")),
										tags$li(p(strong("Gene Annotation: "),"The functions performed by each of the genes that have been selected."))
								)
						)
				)
		)
)

# Define server logic required to draw a histogram
server <- function(input, output) {
	
	# Read CEL files
	data <- reactive({
				validate(
						need(input$CEL != "", "Please select a data set")
				)
				list <- list.files(input$CEL,full.names=TRUE)
				read.celfiles(list)
			})
	
	# Read target file
	phenodata <- reactive({
				validate(
						need(input$target != "", "Please select a target file")
				)
				file1 = input$target
				data1 = read.csv(file1$datapath)
				return(data1)
			})
	
	# Process RMA
	data.rma <- reactive({
				d <- data()
				rma(d)
			})
	
	# Data Matrix
	data.matrix <- reactive({
				d <- data.rma()
				exprs(d)
			})
	
	# Fit Data
	data.fit <- reactive({
				f <- phenodata()
				groups <- f[,2]
				vect <- unique(groups)
				fac <- factor(groups,levels=vect)
				design <- model.matrix(~ 0 + fac)
				colnames(design) <- vect
				df <- lmFit(data.matrix(),design)
				return(df)
			})
	
	# Select groups to compare
	gr_sel <- reactive({
				num <- 1
				control <- FALSE
				f <- phenodata()
				gro <- as.character(f[,2])
				vec <- unique(gro)
				cont <- length(vec)
				h <- ""
				v <- vector()
				m <- vector()
				for (i in 1:cont)
				{
					v[i] <- vec[i]
				}
				for (i in 1:cont)
				{
					for (j in cont:i)
					{
						if(v[i] != v[j])
						{
							h <- paste(v[i],v[j], sep = "-")
							m[num] <- h
							num <- num+1
						}
					}
				}
				return(m)
			})
	
	# eBayes to Fit Data
	data.fit.eb <- reactive({
				f <- phenodata()
				groups <- f[,2]
				vect <- unique(groups)
				fac <- factor(groups,levels=vect)
				design <- model.matrix(~ 0 + fac)
				colnames(design) <- vect
				gru <- gr_sel()
				df <- lmFit(data.matrix(),design)
				contrast <- makeContrasts(contrasts = gru,levels=design)
				datafitcon = contrasts.fit(df,contrast)
				dfeb = eBayes(datafitcon)
				return(dfeb)
			})
	
	# Top Table
	result <- reactive({
				options(digits=10)
				cf <- as.numeric(input$groups)
				d <- data()
				num <- nrow(d@featureData@data)
				topgenes <- topTable(data.fit.eb(),coef=cf, number = num)
				ID <- row.names(topgenes)
				topgenes <- cbind(topgenes, ID)
				annota <- input$annotation
				if(annota != "other") {
					library(annota, character.only = TRUE)
					contador <- nchar(annota)
					if(substr(annota, contador-2,contador) == ".db")
					{
						annota_lib <- substr(annota, 1, contador-3)
					}
					else {
						annota_lib <- annota
					}
					Gene.symbol <- getSYMBOL(rownames(topgenes), annota_lib)
				}
				else {
					annotacion <- input$ann
					data1 = read.table(annotacion$datapath, sep="\t", header = FALSE, fill=TRUE, stringsAsFactors = FALSE)
					data1 = data1[ ,c(1,3)]
					names(data1) <- c("ID", "Gene.symbol")
				}
				if(annota != "other") {
					results <- cbind(Gene.symbol, topgenes)
				} else {
					results <- merge(data1, topgenes, by = "ID")
				}
				results <- na.omit(results)
				row.names(results) <- results$ID
				return(results)
			})
	
	# Selected Genes
	selected <- reactive({
				tab <- result()
				lfc_n <- as.numeric(input$lfc) * -1
				lfc_p <- input$lfc 
				max_gen <- input$max_genes
				topgen <- tab[tab[, "adj.P.Val"] < input$pval, ]
				topgen1 <- topgen[topgen[, "logFC"] < lfc_n,] 
				topgen2 <- topgen[topgen[, "logFC"] > lfc_p, ]
				topgenes3 <- rbind.data.frame(topgen1,topgen2)
				topgenes3 <- topgenes3[order(topgenes3$adj.P.Val), ]
				topgen4 <- na.omit(topgenes3)
				return(topgen4)
			})
	
	
	# COLLECTING DATA
	
	# CEL
	
	
	output$gro <- renderText({
				f <- phenodata()
				grou <- as.character(length(levels(f[,2])))
				grou
			})
	
	output$grou <- renderText({
				f <- phenodata()
				gr <- as.character(levels(f[,2]))
				ob <- paste(gr, " ", sep = ";")
				ob
			})
	
	# SHOW CELFILES
	output$celfiles <- renderTable({
				validate(
						need(input$CEL != "", "Please select a data set")
				)
				l <- list.files(input$CEL,full.names=TRUE)
				h <- strsplit(l, "/", fixed=TRUE)
				Files.CEL <- vector()
				for(i in 1:length(h)){
					a <- length(h[[i]])
					Files.CEL <- c(Files.CEL, h[[i]][a])
				}
				as.data.frame(Files.CEL)
			})
	
	# PHENODATA
	output$targetfile <- renderTable({
				f <- phenodata()
				f
			})
	
	#BOXPLOT BEFORE
	output$box_bef <- renderPlot({
				d <- data()
				f <- phenodata()
				boxplot(d,which='pm',col=as.character(f[,4]),names=f[,3], target='core', main="Before", las=2)
			})
	
	#BOXPLOT AFTER
	output$box_aft <- renderPlot({
				f <- phenodata()
				dm <- data.matrix()
				boxplot(dm,names=f[,3], col=as.character(f[,4]),main="After", las=2)
			})
	
	# DENDOGRAM BEFORE
	output$den_bef <- renderPlot({
				d <- data()
				f <- phenodata()
				eset <- exprs(d)
				distance <- dist(t(eset), method = "maximum")
				clusters <- hclust(distance)
				plot(clusters, main="Before", labels = f[,3])
			})
	
	# DENDOGRAM AFTER
	output$den_aft <- renderPlot({
				f <- phenodata()
				dr <- data.rma()
				eset_rma <- exprs(dr)
				distance_rma <- dist(t(eset_rma), method = "maximum")
				clusters_rma <- hclust(distance_rma)
				plot(clusters_rma, main="After", labels = f[,3])
			})
	
	# VOLCANO PLOT
	output$volcan <- renderPlot({
				dfeb <- data.fit.eb()
				cf <- as.numeric(input$groups)
				annota <- input$annotation
				if(annota != "other") {
					a <- result()
					volcanoplot(dfeb,coef=cf, highlight=input$volcano, main = "Volcano Plot", names = a$Gene.symbol)
				} else { 
					volcanoplot(dfeb,coef=cf, highlight=input$volcano, main = "Volcano Plot", names = names(dfeb$coefficients[,1]))
				}
			})
	
	# HEAT MAP
	output$heatm <- renderPlot({
				res0 <- selected()
				max_gen <- input$max_genes
				res <- res0[1:max_gen, ]
				res <- na.omit(res)
				res1 <- res[res[, "logFC"] < -1, ]
				res2 <- res[res[, "logFC"] > 1, ]
				res <- rbind.data.frame(res1,res2)
				f <- phenodata()
				dm <- data.matrix()
				dm2 <- dm[(rownames(res)),]
				sampleNames <- vector()
				featureNames <- vector()
				heatlogs <- vector()
				num_samples <- nrow(f)
				for (i in 1:num_samples)
				{
					sampleNames <- c(sampleNames,rep(f[i,3],dim(res)[1]))
					featureNames <- c(featureNames,rownames(dm2[1:dim(res)[1],]))
					heatlogs <- c(heatlogs,dm2[1:dim(res)[1],i])
					
				}
				vector = as.character(f[,3])
				for (i in 1:num_samples)
				{
					sampleNames<-replace(sampleNames,sampleNames==i,vector[i])
					
				}
				annota <- input$annotation
				if(annota != "other") {
					heatData <- data.frame(norm_logInt=heatlogs,sampleName=sampleNames,featureName=res$Gene.symbol)
				}
				else {
					heatData <- data.frame(norm_logInt=heatlogs,sampleName=sampleNames,featureName=substr(res$Gene.symbol, 1, 10))
				}
				dataHeat <- ggplot(heatData, aes(sampleName,featureName))
				dataHeat + geom_tile(aes(fill=norm_logInt)) + scale_fill_gradient(low="green", high="red")
			})
	
	# GENE SELECTION (TopTable)
	output$selected1 <- renderTable( digits = -3,{
				res <- result()
				head(res, 10)
			})
	
	output$dwtable <- downloadHandler(filename = function() {
				"TableGenes.csv" },
			content = function(file) {
				write.csv(result(), file)
			})
	
	# GENE SELECTION (SelectedGenes)
	output$selected2 <- renderTable( digits = -3,{
				res2 <- selected()
				head(res2, 10)
			})
	
	output$dwtable2 <- downloadHandler(filename = function() {
				"TableGenes.csv" },
			content = function(file) {
				write.csv(selected(), file)
			})
	
	# UP GENES SELECTION
	output$upgenes <- renderTable(digits = -3,{
				res <- selected()
				topups = res[res[, "logFC"] > 1, ]
				dim(topups)
				head(topups)
			})
	
	output$dwup <- downloadHandler(filename = function() {
				"TableGenes.csv" },
			content = function(file) {
				res <- selected()
				topups = res[res[, "logFC"] > 1, ]
				write.csv(topups, file)
			})
	
	# DOWN GENES SELECTION
	output$downgenes <- renderTable(digits = -3,{
				res <- selected()
				topdown = res[res[, "logFC"] < -1, ]
				dim(topdown)
				head(topdown)
			})
	
	output$dwdw <- downloadHandler(filename = function() {
				"TableGenes.csv" },
			content = function(file) {
				res <- selected()
				topdown = res[res[, "logFC"] < -1, ]
				write.csv(topdown, file)
			})
	
	# Gene Ontology
	output$go2 <- renderTable({
				require(GOstats)
				annota <- input$annotation
				bd <- character()
				p <- as.numeric(input$pvalGO)
				if(annota != "other")
				{
					bd <- annota
				}
				top <- result()
				genes <- selected()
				entrezUniverse = unique(getEG(as.character(rownames(top)), bd))
				geneIds <- unique(getEG(as.character(rownames(genes)),bd))
				GOparams = new("GOHyperGParams",
						geneIds=geneIds, universeGeneIds=entrezUniverse,
						annotation=bd, ontology=input$ontol,
						pvalueCutoff=p, conditional=FALSE,
						testDirection="over")
				KEGGparams = new("KEGGHyperGParams",
						geneIds=geneIds, universeGeneIds=entrezUniverse,
						annotation=bd,
						pvalueCutoff=p, testDirection="over")
				GOhyper = hyperGTest(GOparams)
				head(summary(GOhyper), 20)
			})
	
	output$dwgo2 <- downloadHandler(filename = function() {
				"myreport.html" },
			content = function(file) {
				require(GOstats)
				annota <- input$annotation
				bd <- character()
				p <- as.numeric(input$pvalGO)
				if(annota != "other")
				{
					bd <- annota
				}
				top <- result()
				genes <- selected()
				entrezUniverse = unique(getEG(as.character(rownames(top)), bd))
				geneIds <- unique(getEG(as.character(rownames(genes)),bd))
				GOparams = new("GOHyperGParams",
						geneIds=geneIds, universeGeneIds=entrezUniverse,
						annotation=bd, ontology=input$ontol,
						pvalueCutoff=p, conditional=FALSE,
						testDirection="over")
				KEGGparams = new("KEGGHyperGParams",
						geneIds=geneIds, universeGeneIds=entrezUniverse,
						annotation=bd,
						pvalueCutoff=p, testDirection="over")
				GOhyper = hyperGTest(GOparams)
				# Creamos un informe html con los resultados
				htmlReport(GOhyper, file = file, summary.args=list("htmlLinks"=TRUE))
			}) 
	
	
	# Gene Annotation
	output$go <- renderTable({
				res <- selected()
				a <- as.vector(res$Gene.symbol)
				annota <- input$annotation
				if(annota != "other") {
					if(annota == "mogene10sttranscriptcluster.db")
					{
						table1 <- select(mogene10sttranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "hgu133a.db")
					{
						table1 <- select(hgu133a.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "hugene21sttranscriptcluster.db")
					{
						table1 <- select(hugene21sttranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "hugene20sttranscriptcluster.db")
					{
						table1 <- select(hugene20sttranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "clariomdhumantranscriptcluster.db")
					{
						table1 <- select(clariomdhumantranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "clariomshumanhttranscriptcluster.db")
					{
						table1 <- select(clariomshumanhttranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "clariomshumantranscriptcluster.db")
					{
						table1 <- select(clariomshumantranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "clariomsmousehttranscriptcluster.db")
					{
						table1 <- select(clariomsmousehttranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "clariomsmousetranscriptcluster.db")
					{
						table1 <- select(clariomsmousetranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					b <- table1$GO
				} else {
					annotacion <- input$ann
					data1 = read.table(annotacion$datapath, sep="\t", header = FALSE, fill=TRUE, stringsAsFactors = FALSE)
					data1 = data1[ ,c(1,21)]
					names(data1) = c("ID", "GO")
					res <- merge(res, data1, by = "ID")
					res <- na.omit(res)
					row.names(res) <- res$ID
					a <- res$ID
					c <- res$Gene.symbol
					d <- res$GO
					d <- substr(d, start = 1, stop = 10)
					table1 <- data.frame(a,c,d)
					names(table1) <- c("ID", "Gene.symbol", "GO")
					b <- table1$GO
					b <- sapply(b, as.character)
					b <- substr(b, start = 1, stop = 10)
				}
				table2 <- select(GO.db, keys=b, columns=c("TERM", "DEFINITION"),
						keytype="GOID")
				names(table2) <- c("GO", "TERM", "DEFINITION")
				table3 <- merge(table1,table2, by = "GO")
				table3 <- na.omit(table3)
				table3 <- unique(table3)
				head(table3, 10)
			})
	
	output$dwgo <- downloadHandler(filename = function() {
				"TableGenes.csv" },
			content = function(file) {
				res <- selected()
				a <- as.vector(res$Gene.symbol)
				annota <- input$annotation
				if(annota != "other") {
					if(annota == "mogene10sttranscriptcluster.db")
					{
						table1 <- select(mogene10sttranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "hgu133a.db")
					{
						table1 <- select(hgu133a.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "hugene21sttranscriptcluster.db")
					{
						table1 <- select(hugene21sttranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "hugene20sttranscriptcluster.db")
					{
						table1 <- select(hugene20sttranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "clariomdhumantranscriptcluster.db")
					{
						table1 <- select(clariomdhumantranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "clariomshumanhttranscriptcluster.db")
					{
						table1 <- select(clariomshumanhttranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "clariomshumantranscriptcluster.db")
					{
						table1 <- select(clariomshumantranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "clariomsmousehttranscriptcluster.db")
					{
						table1 <- select(clariomsmousehttranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					else if(annota == "clariomsmousetranscriptcluster.db")
					{
						table1 <- select(clariomsmousetranscriptcluster.db, keys=a, columns=c("SYMBOL","GENENAME", "GO"),
								keytype="SYMBOL")
					}
					b <- table1$GO
				} else {
					annotacion <- input$ann
					data1 = read.table(annotacion$datapath, sep="\t", header = FALSE, fill=TRUE, stringsAsFactors = FALSE)
					data1 = data1[ ,c(1,21)]
					names(data1) = c("ID", "GO")
					res <- merge(res, data1, by = "ID")
					res <- na.omit(res)
					row.names(res) <- res$ID
					a <- res$ID
					c <- res$Gene.symbol
					d <- res$GO
					d <- substr(d, start = 1, stop = 10)
					table1 <- data.frame(a,c,d)
					names(table1) <- c("ID", "Gene.symbol", "GO")
					b <- table1$GO
					b <- sapply(b, as.character)
					b <- substr(b, start = 1, stop = 10)
				}
				table2 <- select(GO.db, keys=b, columns=c("TERM", "DEFINITION"),
						keytype="GOID")
				names(table2) <- c("GO", "TERM", "DEFINITION")
				table3 <- merge(table1,table2, by = "GO")
				table3 <- na.omit(table3)
				table3 <- unique(table3)
				write.csv(table3, file)
			})
	
	output$prueba <- renderText({
				
			})
	
}

# Run the application 
shinyApp(ui = ui, server = server)

