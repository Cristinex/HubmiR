args = commandArgs(trailingOnly=TRUE)
library('Seurat')
assays=args[2]
output_dir=args[3]
L1000<-read.csv(file='../files/L1000gene.csv',header=FALSE)
vct <- as.vector(unlist(L1000$V1))
vct[1] <- "PSME1"
row.names(L1000) <- L1000$V1
scRNAlist <- list()
to_merge<- c()
for (i in c(1:length(assays))){
	file <- paste(assays,assays[i],sep='')
	data <- read.csv(file=file,header=TRUE,row=1)
	input <- data[-1,]
	input <-as.data.frame(input)
	input[is.na(input)] <- 0
	type<-as.matrix(data[1,])
	label <- strsplit(assays[i],split='_')
	label <- strsplit(label[[1]][2],split='.csv')[1]
	for (cell in 1:range(length(type))){
		type[cell] <- paste(type[cell],'_',label,sep = '')
		type[cell] <- paste(strsplit(type[cell],' ')[[1]],collapse = '.')

	}
	scRNAlist[[i]]<- CreateSeuratObject(input,min.cells=0,min.features=0)
	Idents(scRNAlist[[i]])<-type
	
	if (i > 1){
		to_merge <- append(to_merge,scRNAlist[[i]])
	}
	
}

if (length(scRNAlist) >  1){
	seraobj <- merge(scRNAlist[[1]], y= to_merge)
} else {
	seraobj <- scRNAlist[[1]]
}
seraobj <- NormalizeData(seraobj, normalization.method = "LogNormalize", scale.factor = 10000)
seraobj<-ScaleData(object=seraobj)
seraobj <- FindVariableFeatures(seraobj,selection.method = "vst", nfeatures = 2000)
if (length(Cells(seraobj))<50){
    seraobj <- RunPCA(seraobj,features = VariableFeatures(object = seraobj),npcs=5)
    seraobj <- RunUMAP(seraobj,dims=1:2,n.neighbors = 3)
} else{
    seraobj <- RunPCA(seraobj,features = VariableFeatures(object = seraobj))
    seraobj <- RunUMAP(seraobj,dims=1:10)
}
cell_num<- as.data.frame(table(Idents(seraobj)))
if (args[1] == "celltypeavg"){
	cluster.averages<-AverageExpression(seraobj)
	clusteravg <- cluster.averages[["RNA"]]

	pooled_sample <- as.data.frame(matrix(nrow = length(row.names(cell_num)), ncol = 1))
	row.names(pooled_sample)<-cell_num$Var1
	pooled_sample$V1 <- rep(1,length(row.names(pooled_sample)))

} else{
	cell_list <- as.data.frame(Idents(seraobj))
	cell_list$'V2' <- row.names(cell_list)
	expression_profile <- GetAssayData(seraobj,slot='data')
	pooled_sample <- as.data.frame(matrix(nrow = length(row.names(cell_num)), ncol = 1))
	row.names(pooled_sample)<-cell_num$Var1
	for (i in c(1:length(row.names(cell_num)))){
		num<-cell_num[i,2]%/%200
		pooled_sample[i,1] <- num
		cells <-cell_list[which(cell_list$'Idents(seraobj)'==cell_num[i,1]),2]
		cell_profile <- expression_profile[,cells]
		if (num == 0){
			sampled_cells <- sample(cells,200,replace=TRUE)
			sum_profile <-  cell_profile[,sampled_cells]
			expr_sum <- as.data.frame(rowSums(sum_profile))
			expr_sum <- expr_sum/200
			colnames(expr_sum) <- as.character(cell_num[i,1])
		} else if (num == 1){
			sampled_cells <- sample(cells,200,replace=FALSE)
			sum_profile <-  cell_profile[,sampled_cells]
			expr_sum <- as.data.frame(rowSums(sum_profile))
			expr_sum <- expr_sum/200
			colnames(expr_sum) <- cell_num[i,1]
		} else if (num >= 2){
			for (j in c(1:num)){
				sampled_cells <- sample(cells,200,replace=FALSE)
				if (j == 1){
					sum_profile <-  cell_profile[,sampled_cells]
					expr_sum <- as.data.frame(rowSums(sum_profile))
					expr_sum <- expr_sum/200
					colnames(expr_sum) <- paste0(cell_num[i,1],'_',j)
					remove_cells <- sample(sampled_cells, 160,replace=FALSE)
					cells <- unlist(cells[cells %in% remove_cells == FALSE])
				} else {
					sum_profile <-  cell_profile[,sampled_cells]
					expr_sum_2 <- as.data.frame(rowSums(sum_profile))
					expr_sum_2 <- expr_sum_2/200
					colnames(expr_sum_2) <- paste0(cell_num[i,1],'_',j)
					expr_sum <- cbind(expr_sum,expr_sum_2)
					remove_cells <- sample(sampled_cells, 160,replace=FALSE)
					cells <- unlist(cells[cells %in% remove_cells == FALSE])
				}
			}
		}
		
		if (i == 1){
		clusteravg = expr_sum
		} else {
		clusteravg <- cbind(clusteravg,expr_sum)
		}
	}
	row.names(pooled_sample) <- cell_num$'Var1'
}


vectorIN<-list()
vectorOUT<-list()

vectorIN <- vct[vct %in% row.names(clusteravg)]
vectorOUT <- vct[!(vct %in% row.names(clusteravg))]
data1 <- clusteravg[vectorIN,]
if (length(vectorOUT) > 0){
	candidate <- L1000[vectorOUT,2]
	vectorIN_2 <- candidate[candidate %in% row.names(clusteravg)]
	if (length(vectorIN_2) <= 0){
		data2 <- data.frame(matrix(nrow=length(vectorOUT),ncol=length(data1[2,])))
		colnames(data2) <- colnames(data1)
		row.names(data2) <- vectorOUT
		output <- rbind(data1,data2)
	} else{
		candidate_out <- candidate[!(candidate %in% vectorIN_2)]
		vectorOUT <- L1000[which(L1000$V2 %in% candidate_out),1]
		vectorIN_replace<- L1000[which(L1000$V2 %in% vectorIN_2),1]
		data3 <- clusteravg[vectorIN_2,]
		data3 <- t(as.data.frame(data3))
		row.names(data3) <- vectorIN_replace
		colnames(data3) <- colnames(data1)
		data2 <- data.frame(matrix(nrow=length(vectorOUT),ncol=length(data1[2,])))
		colnames(data2) <- colnames(data1)
		row.names(data2) <- vectorOUT
		output <- rbind(data1,data2,data3)
		
 	}	
} else{
	output <- data1
	
}

output[is.na(output)] <- 0
output_final = t(as.data.frame(output[vct,]))

output_scale <- as.data.frame(matrix(nrow=length(row.names(output_final)),ncol=length(colnames(output_final))))
for (i in c(1:length(colnames(output_final)))){
	out <- scale(as.vector(output_final[,i]))
	output_scale[,i] <- out
}
colnames(output_scale) <- colnames(output_final)
row.names(output_scale) <- row.names(output_final)
final<-as.data.frame(output_scale)
final[is.na(final)] <- 0
setwd(output_dir)
write.csv(final,file='extracted_input.csv')
write.csv(pooled_sample,file='pooling_info.csv')
saveRDS(seraobj,file='single_cell_object.rds')
png(file='UMAP.png',width=2000,heigh=2000,bg="transparent",res=200)
DimPlot(seraobj, reduction = "umap")
dev.off()





