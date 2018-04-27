library(oligo)
library(limma)
library(gplots)
library(pd.ragene.2.0.st)
library(calibrate)
library(oligoClasses)
#library(xps)
library(calibrate)
library("ragene10sttranscriptcluster.db")
library(AnnotationDbi)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(grid)
#
biocLite("ragene10stv1cdf")

################################# LEER DATOS ###############################################

files=c("TAC_2.CEL","TAC_3.CEL","TAT_1.CEL","TAT_2.CEL","TAT_3.CEL")
adipo=read.celfiles(files)
hist(hipo, target= "probeset",main= "datos crudos")
boxplot(hipo, target= "probeset",main= "datos crudos")

############################     NORMALIZAR LOS DATOS #######################################
hip.rma.probe=rma(hipo, target="probeset")
expresion<-exprs(hip.rma.probe)
hist(hip.rma.probe, target= "probeset",main= "datos normalizadosR")
boxplot(probe, target= "probeset",main= "datos normalizados")

############################      PCA                ###########################################
library(factoextra)
res.pca <- prcomp(t(as.matrix(hip.rma.probe)), scale = TRUE)
fviz_eig(res.pca)
         groups<-c("CONTROL","CONTROL", "CONTROL", "SM" , "SM" , "SM" )

fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE,
             title = "AT_Stem Cells: *Control vs MS*"
             )

#############################     CONTRATE ESTADISTICO     #######################################

design=matrix(0,12,2)
	design[1:6,1]=1
	design[7:12,2]=1
	colnames(design)=c("Control","SM")
contraste="SM-Control"
cont.matrix=makeContrasts(contraste,levels=design)
fit=lmFit(hip.rma.probe, design)
fitC=contrasts.fit(fit, cont.matrix)
fitCB=eBayes(fitC)
TT=topTable(fitCB, coef=1, adjust="fdr",sort.by="logFC",number=nrow(hip.rma.probe),genelist=fit$genes)
head(TT)

##########################      ANOTACIONES                ########################################

library("ragene10sttranscriptcluster.db")
library(AnnotationDbi)
annots <- select(ragene10sttranscriptcluster.db, keys=rownames(TT),
columns=c("SYMBOL","GENENAME"), keytype="PROBEID")
## ’select()’ returned 1:1 mapping between keys and columns
resultTable <- merge(TT, annots, by.x=0, by.y="PROBEID")
head(resultTable)
write.csv(resultTable, "TT_gen2.csv")

#####################  VOLVANO PLOT    para expresión diferencia      ########################################

library(calibrate)
with(resultTable, plot(logFC, -log10(P.Value), pch=20, main="Adipose tissue", repel=TRUE,xlim=c(-4,4)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
	abline(h= 1.30103,col="black")
	abline(v=-.5,col="black")
	abline(v=.5,col="black")
	#abline(h= 2,col="black", lwd=(5))
	#abline(v=-2,col="black",lwd=(5))
	#bline(v=2,col="black",lwd=(5))


	with(subset(resultTable, logFC>.5), points(logFC, -log10(P.Value), pch=20, col="darkorange"))
	with(subset(resultTable, logFC<.5), points(logFC, -log10(P.Value), pch=20, col="cyan"))
	with(subset(resultTable, abs(logFC)<.5), points(logFC, -log10(P.Value), pch=20, col="black"))
    with(subset(resultTable, P.Value< 0.05 ), points(logFC, -log10(P.Value), pch=20, col="grey"))
    with(subset(resultTable, P.Value< 0.05 & logFC< -.5), points(logFC, -log10(P.Value), pch=20, col="blue"))
	with(subset(resultTable, P.Value< 0.05 & logFC> .5), points(logFC, -log10(P.Value), pch=20, col="red"))
	with(subset(resultTable, P.Value< 0.0001 & logFC> 0.5),textxy(logFC, -log10(P.Value), labs=SYMBOL, cex=.8, repel=TRUE))
	with(subset(resultTable, P.Value< 0.0001 & logFC< 0.5),textxy(logFC, -log10(P.Value), labs=SYMBOL, cex=.8,repel=TRUE)) 
	
##################### VOLCANO PLOT 1  #####################################            
library(ggplot2)
library(ggrepel)
library(grid)

resultTable$Significant <- ifelse(resultTable$P.Value > 0.05,"logFC > 0.5", " Not sig")
ggplot(resultTable, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red","blue")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(resultTable, P.Value < 0.0001&abs(logFC)>1),
    aes(label = SYMBOL),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )


#scale_color_continuous(
#"logFC",with(resultTable,c(min(logFC),mean(logFC),max(logFC))),
labels = c("Bajo","Medio","Alto"))


###########################  SELECCIÓN DE EXPRESIÓN DIFERENCIAL   ######################################

selected<-resultTable[abs(resultTable$logFC)>.5,]
selected_symbol<-selected[abs(selected$P.Value)<.05,]

symbol<-selected_symbol$SYMBOL
dim(symbol)

selected<-TT[abs(TT$logFC)>.5,]
selected2<-selected[abs(selected$P.Value)<.05,]
dim(selected2)