
##########Running validated 300 gene KRAS signatures on  samples####

library(ASSIGN)
library(sva)


# get this file
gfp_kras<-read.table("~/Documents/ThesisWork/GitRepos/Lung_Cancer_Project_KRAS/36_GFP_KRAS_TPMlog2.txt", sep='\t', header=1, row.names=1)
head(gfp_kras)

dim(gfp_kras)
CCLE<-read.table("~/Documents/ThesisWork/CCLE_Data/CCLE_RNASeq_TPMlog.txt",header=1, row.names=1, sep='\t',check.names = F)
CCLE_lung<-subset(CCLE,select=c("5bea08d5-6a97-41fa-b33c-d8cb45b67051", "39c19460-909e-47df-892d-86638ddd5969", "ec0491f2-e5f8-4dd4-9ebe-197916a04d6d","b759a72e-d81b-45f1-bc66-06e3d0f2f600", "20d36af9-655c-45e7-8f32-9770dd0c53c5", "d77cb0f1-9831-4613-a29a-5796b2d3754f", "73045153-e0f8-43a6-ae22-f1ecd7ce775a", "c83c157a-080c-478b-8aeb-d3104b431af0", "a1e1283d-14e1-4ff3-9e0f-d9c468e78dc4", "e619fbdb-2cd6-46a4-8a18-7fd1e23a31cc"))
head(CCLE_lung)
colnames(CCLE_lung)
colnames(CCLE_lung) = c("H838", "H1437", "HCC4006", "H1944", "SKLU-1", "H1650", "H460", "H1395", "H1355", "H2009")
head(CCLE_lung)
dim(CCLE_lung)
write.table(CCLE_lung,"~/Documents/ThesisWork/GitRepos/Lung_Cancer_Project_KRAS/CCLE_Lung_DrugScreen_RNASeq_TPMlog.txt",sep='\t', col.names = NA,quote=F)

sub=c(9,9,9,9,10)
expr_f <-gfp_kras[apply(gfp_kras[,1:36]==0,1,mean) < 0.85,]
#merge the PE cells with the signature 
expr<-merge_drop(expr_f,CCLE_lung)
pdf("~/Documents/ThesisWork/GitRepos/Lung_Cancer_Project_KRAS/kras_CCLE_lung_pca_plots.pdf")
pcaplot(expr,sub)

bat<-as.matrix(cbind(colnames(expr),c(rep(1,ncol(gfp_kras)),rep(2,ncol(CCLE_lung)))))
# Batch adjust 
combat_expr<-ComBat(dat=expr, batch=bat[,2], mod=NULL, numCovs=NULL)
pdf("~/Documents/ThesisWork/GitRepos/Lung_Cancer_Project_KRAS/kras_CCLE_lung_pca_plots_post.pdf")
pcaplot(combat_expr,sub)
dev.off()

#subset the diff KRAS signatures and the test data
c_kras_gfp<-subset(combat_expr,select=GFP.31:GFP.39)
c_kraswt<-subset(combat_expr,select=KRASWT.1:KRASWT.9)
c_krasqh<-subset(combat_expr,select=KRASQH.1:KRASQH.9)
c_krasgv<-subset(combat_expr,select=KRASGV.1:KRASGV.9)
c_test<-combat_expr[,(ncol(gfp_kras)+1):ncol(combat_expr)]
colnames(c_test)

##getting the gene list###
load("~/Documents/ThesisWork/GitRepos/Lung_Cancer_Project_KRAS/Kras_300s/krasqv_300_gene_list/adapB_single/output.rda")
krasqh_300_genelist<-output.data$processed.data$diffGeneList
head(krasqh_300_genelist)

load("~/Documents/ThesisWork/GitRepos/Lung_Cancer_Project_KRAS/Kras_300s/krasgv_300_gene_list/adapB_single/output.rda")
krasgv_300_genelist<-output.data$processed.data$diffGeneList
head(krasgv_300_genelist)

load("~/Documents/ThesisWork/GitRepos/Lung_Cancer_Project_KRAS/Kras_300s/kraswt_300_gene_list/adapB_single/output.rda")
kraswt_300_genelist<-output.data$processed.data$diffGeneList
head(kraswt_300_genelist)

basedir="~/Documents/ThesisWork/GitRepos/Lung_Cancer_Project_KRAS/predictions"
dir.create(basedir)

trainingLabel<-list(control=list(kraswt=1:9),kraswt=(10:18))
dir.create(paste(basedir,paste("kraswt",300,"gene_list", sep="_"),sep='/'))
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_kraswt),test=c_test,trainingLabel1 = trainingLabel,geneList = kraswt_300_genelist,out_dir_base = paste(basedir,paste("kraswt",300,"gene_list", sep="_"),sep='/'),single = 1)

trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
dir.create(paste(basedir,paste("krasgv",300,"gene_list", sep="_"),sep='/'))
trainingLabel<-list(control=list(krasgv=1:9),krasgv=(10:18))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_krasgv),test=c_test,trainingLabel1 = trainingLabel,geneList = krasgv_300_genelist,out_dir_base = paste(basedir,paste("krasgv",300,"gene_list", sep="_"),sep='/'),single = 1)

trainingLabel<-list(control=list(krasqh=1:9),krasqh=(10:18))
dir.create(paste(basedir,paste("krasqh",300,"gene_list", sep="_"),sep='/'))
assign_easy_multi(trainingData = cbind(c_kras_gfp,c_krasqh),test=c_test,trainingLabel1 = trainingLabel,geneList = krasqh_300_genelist,out_dir_base = paste(basedir,paste("krasqh",300,"gene_list", sep="_"),sep='/'),single = 1)

pe_preds<-gatherFile(basedir)
colnames(pe_preds)<-gsub(colnames(pe_preds),pattern = "/pathway_activity_testset.csv/V1",replacement = "")
colnames(pe_preds)<-gsub(colnames(pe_preds),pattern = "/pathway_activity_testset.csv",replacement = "")
rownames(pe_preds)

write.table(pe_preds,"~/Documents/ThesisWork/GitRepos/Lung_Cancer_Project_KRAS/CCLE_Lung_KRAS_Predictions",sep='\t', col.names = NA,quote=F)




# Not using

for(i in 1:nrow(pe_preds)){
  pe_preds$new_names[i]<-strsplit(rownames(pe_preds)[i],"_")[[1]][2]
}
library(xlsx)
f=read.xlsx("~/Dropbox/bild_signatures/pe_preds/Copy of PE_Gray_DNA_RNA.xls",sheetName = "Sheet1")
pe_pts<- f[1:14,]
pe_pts[,1]<-toupper(pe_pts[,1])
order(rownames(pe_pts))
rownames(pe_pts)<-pe_pts$RNA
pe_pts<-pe_pts[order(rownames(pe_pts)),]
pe_preds<-pe_preds[order(pe_preds$new_names),]
pe_pred_phen<-merge(pe_preds,pe_pts,by.x=57,by.y=1)
dim(pe_pred_phen)
gsub("/SigProtein","",colnames(pe_pred_phen))
colnames(pe_pred_phen)
colnames(pe_pred_phen)[67]<-"PR.Status"
colnames(pe_pred_phen)[68]<-"HER2.Status"
cols<-gsub("/sigProtein","",c(colnames(single_pathway_best),colnames(multi_pathway_best)))
subset(pe_pred_phen,select=cols)#,"ER.Status","PR.Status","HER2.Status")]


################## for correlations example codes####
drugs<-read.delim("ICBP_drugs.txt", header=1, sep='\t',row.names=1)
icbp_drug<-merge_drop(data_icbp,drugs)
colnames(icbp_drug)
cor_mat=p_mat=matrix(0,length(filenames_icbp_multi),90)
rownames(cor_mat)=rownames(p_mat)=colnames(icbp_drug)[1:length(filenames_icbp_multi)]
colnames(cor_mat)=colnames(p_mat)=colnames(icbp_drug)[(length(filenames_icbp_multi)+11):ncol(icbp_drug)]

for(i in 1:length(filenames_icbp_multi)){
  for(j in 1:90){
    temp=cor.test(icbp_drug[,i],icbp_drug[,(j+length(filenames_icbp_multi)+10)],use="pairwise",method="spearman")
    print(j)
    print(temp)
    cor_mat[i,j]=temp$estimate
    p_mat[i,j]=temp$p.value
  }
}
#write.table(cor_mat,"~/Dropbox/bild_signatures/kras/ICBP_single_cor_drug_mat_4_21.txt",col.names = NA,quote=F,sep='\t')
colnames(p_mat)=paste(colnames(p_mat),"p_value",sep="_")
write.table(p_mat,"~/Dropbox/bild_signatures/kras/ICBP_single_p_drug_mat_4_21.txt",col.names = NA,quote=F,sep='\t')
cor_p_mat<-cbind(cor_mat,p_mat)
order(colnames(cor_p_mat))
cor_p_mat<-cor_p_mat[,order(colnames(cor_p_mat))]
write.table(cor_p_mat,"~/Dropbox/bild_signatures/kras/ICBP_single_cor_p_mat_6_18.txt",col.names = NA,quote=F,sep='\t')



