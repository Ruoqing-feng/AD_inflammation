library(Seurat)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(stringr)
library(sctransform)
library(pheatmap)
library(RColorBrewer)

all = list.files("/Volumes/T7/Analysis/",pattern = "_")
all = all[2:7]
all = paste0("/Volumes/T7/Analysis/",all,"/")

sampel = paste0("L0",1:6)
all.mtx = Read10X_h5(paste0(all[1],"filtered_feature_bc_matrix.h5"))
colnames(all.mtx) = paste(sampel[1],colnames(all.mtx),sep = "_")

for (i in 2:6){
  print(i)
  tmp = Read10X_h5(paste0(all[i],"filtered_feature_bc_matrix.h5"))
  colnames(tmp) = paste(sampel[i],colnames(tmp),sep = "_")
  all.mtx = cbind(all.mtx,tmp)
}

fad.seurat = CreateSeuratObject(all.mtx)

fad.seurat[["percent.mt"]] <- PercentageFeatureSet(fad.seurat, pattern = "^mt-")
fad.seurat@meta.data$project = "5XFAD"
Idents(fad.seurat) = "project"

fad.seurat = subset(fad.seurat,percent.mt < 10 & nCount_RNA <= 50000  & nFeature_RNA <= 7500)
pdf("/Users/ruoqing/Projects/5xFAD_Tdepletion/cutoff.pdf",width = 5,height = 5)
VlnPlot(fad.seurat,features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = 0)
dev.off()

fad.seurat <- SCTransform(fad.seurat, 
                            vars.to.regress = "percent.mt", 
                            verbose = T,
                            variable.features.rv.th = 1.4,
                            return.only.var.genes = F,
                            variable.features.n = NULL)
fad.seurat@assays$SCT@var.features = setdiff(fad.seurat@assays$SCT@var.features,str_subset(fad.seurat@assays$SCT@var.features, pattern ='^mt-'))

DefaultAssay(fad.seurat) <- "SCT"

fad.seurat <- RunPCA(object = fad.seurat, npcs = 50,verbose = T)
pdf("/Users/ruoqing/Projects/soat1_ko/elbow.pdf",width = 5,height = 5)
ElbowPlot(fad.seurat,ndims = 50)
dev.off()

fad.seurat <- FindNeighbors(fad.seurat, dims = 1:20)
fad.seurat <- FindClusters(fad.seurat, resolution = 0.5)
fad.seurat <- RunUMAP(object = fad.seurat, dims = 1:20, do.fast = TRUE,n.neighbors = 50)

fad.seurat@meta.data$group = "Control"
fad.seurat@meta.data$group[fad.seurat@meta.data$orig.ident %in% sampel[4:6]] = "T_depletion"
DimPlot(fad.seurat,label = T,split.by = "group")

VlnPlot(fad.seurat,features = c("nCount_RNA","nFeature_RNA"),pt.size = 0)
Idents(fad.seurat) = "SCT_snn_res.0.1"

fad.marker = FindAllMarkers(fad.seurat)
FeaturePlot(fad.seurat,features = "Cd8a",split.by = "group")

fad.seurat@meta.data$tmp = as.numeric(fad.seurat@meta.data$SCT_snn_res.0.1)
fad.seurat@meta.data$tmp = fad.seurat@meta.data$tmp - 1

fad.seurat.clean = subset(fad.seurat,tmp %in% c(0:10,12,13))
DimPlot(fad.seurat.clean)
save(fad.seurat,file = "/Users/ruoqing/Projects/5xFAD_Tdepletion/r_obj/fad.seurat.Rdata")
rm(fad.seurat)

fad.seurat.clean@meta.data$annotation = "Microglia"
fad.seurat.clean@meta.data$annotation[fad.seurat.clean@meta.data$tmp == 3] = "EpendySec"
fad.seurat.clean@meta.data$annotation[fad.seurat.clean@meta.data$tmp %in% c(2,6)] = "Astrocytes"
fad.seurat.clean@meta.data$annotation[fad.seurat.clean@meta.data$tmp %in% c(5,13)] = "Endothelial"
fad.seurat.clean@meta.data$annotation[fad.seurat.clean@meta.data$tmp == 7] = "T/NK"
fad.seurat.clean@meta.data$annotation[fad.seurat.clean@meta.data$tmp == 8] = "EpendyCilia"
fad.seurat.clean@meta.data$annotation[fad.seurat.clean@meta.data$tmp == 9] = "Vacular"
fad.seurat.clean@meta.data$annotation[fad.seurat.clean@meta.data$tmp == 10] = "Hippocampus neurons"
fad.seurat.clean@meta.data$annotation[fad.seurat.clean@meta.data$tmp == 12] = "Mural"

Idents(fad.seurat.clean) = "annotation"
DimPlot(fad.seurat.clean,cols = fad.color)

tmp = randomcoloR::distinctColorPalette(10)
scales::show_col(tmp)
fad.color = tmp[2:10]
names(fad.color) = unique(fad.seurat.clean@meta.data$annotation)

fad.seurat.clean@meta.data$cell = rownames(fad.seurat.clean@meta.data)
p1 = DimPlot(fad.seurat.clean)

p1 = DimPlot(subset(fad.seurat.clean,annotation == "T/NK"))
fad.seurat.clean = subset(fad.seurat.clean, cell %in% setdiff(fad.seurat.clean@meta.data$cell,CellSelector(p1)))

Idents(fad.seurat.clean) = "annotation"
p1 = DimPlot(fad.seurat.clean,cols = fad.color,label = T)+
  theme_dr()+
  theme(panel.grid = element_blank(),
        legend.position = "none")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/umap.png",p1,width = 9,height = 6,dpi = 300)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/umap.pdf",p1,width = 9,height = 6,dpi = 300)

p1 = DimPlot(fad.seurat.clean,cols = fad.color,label = F)+
  theme_dr()+
  theme(panel.grid = element_blank(),
        legend.position = "none")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/umap_NL.png",p1,width = 7,height = 6,dpi = 300)

p1 = DimPlot(fad.seurat.clean,cols = fad.color,split.by = "group")+
  theme_dr()+
  theme(panel.grid = element_blank())
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/umap_group.png",p1,width = 15,height = 6,dpi = 300)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/umap_group.pdf",p1,width = 15,height = 6,dpi = 300)

unique(fad.seurat.clean@meta.data$annotation)

fad.marker = FindAllMarkers(fad.seurat.clean)

top10 <- fad.marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

use.gene = top10$gene

tmp = AverageExpression(fad.seurat.clean,assays = "SCT",features = use.gene,return.seurat = T,slot = "data")
DoHeatmap(tmp,features = use.gene)


tmp = tmp@assays$SCT@scale.data

bk = seq(-2,2,0.1)
p1 = pheatmap(tmp,breaks = bk,
              color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)),
              cluster_rows = F,cluster_cols = F,border_color = NA)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/heat_marker.pdf",p1,width = 3,height = 7,dpi = 300)

use.gene = c("P2ry12","C1qa","Tmem119","Gpc5","Slc1a2","Cldn5","Flt1","Cd3e","Nkg7","Dcx","Pcna","Mki67","Mgp","Ptgds","Enpp2","Tmem212","Vtn")
tmp = fad.seurat.clean@assays$SCT@data[use.gene,rownames(fad.seurat.clean@meta.data)] 
tmp = t(as.matrix(tmp))
tmp = cbind(fad.seurat.clean@meta.data,tmp)

tmp = tmp[,c("annotation",use.gene)]
tmp = reshape2::melt(tmp)
tmp$annotation = factor(tmp$annotation,levels = unique(tmp$annotation),ordered = T)
p1 = ggplot(tmp,aes(x = annotation, y = value,fill = annotation))+
  geom_violin(linetype="blank",scale = "width")+
  facet_grid(variable~. ,scales = 'free')+
  scale_fill_manual(values = as.character(fad.color),breaks = names(fad.color))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/Vln_marker.pdf",p1,width = 2,height = 5.5,dpi = 300)

tmp = matrix(NA,ncol = length(unique(fad.seurat.clean@meta.data$annotation)),nrow = length(unique(fad.seurat.clean@meta.data$orig.ident)))
colnames(tmp) = unique(fad.seurat.clean@meta.data$annotation)
rownames(tmp) = unique(unique(fad.seurat.clean@meta.data$orig.ident))

for (i in colnames(tmp)){
  for (j in rownames(tmp)){
    tmp[j,i] = nrow(fad.seurat.clean@meta.data[fad.seurat.clean@meta.data$annotation == i & fad.seurat.clean@meta.data$orig.ident == j,])/nrow(fad.seurat.clean@meta.data[fad.seurat.clean@meta.data$orig.ident == j,])
  }
}
tmp = as.data.frame(tmp)
tmp$group = c(rep("Control",3),rep("T_depletion",3))

tmp = reshape2::melt(tmp)

p1 = ggplot(tmp,aes(x = variable,y = value,fill = group))+
  geom_boxplot(color = "grey66")+
  scale_fill_manual(values = c("#4f22a3","#f2d38c"))+
  theme_classic()+xlab("")+ylab("")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/clu_fraction.pdf",p1, width = 10,height = 4,dpi = 300)


t.fad.seurat = subset(fad.seurat,SCT_snn_res.0.1 == 7)

t.fad.seurat <- SCTransform(t.fad.seurat, 
                          vars.to.regress = "percent.mt", 
                          verbose = T,
                          variable.features.rv.th = 1.4,
                          return.only.var.genes = F,
                          variable.features.n = NULL)
t.fad.seurat@assays$SCT@var.features = setdiff(t.fad.seurat@assays$SCT@var.features,str_subset(t.fad.seurat@assays$SCT@var.features, pattern ='^mt-'))

DefaultAssay(t.fad.seurat) <- "SCT"

t.fad.seurat <- RunPCA(object = t.fad.seurat, npcs = 50,verbose = T)
pdf("/Users/ruoqing/Projects/soat1_ko/elbow.pdf",width = 5,height = 5)
ElbowPlot(t.fad.seurat,ndims = 50)
dev.off()

t.fad.seurat <- FindNeighbors(t.fad.seurat, dims = 1:10)
t.fad.seurat <- FindClusters(t.fad.seurat, resolution = 0.5)
t.fad.seurat <- RunUMAP(object = t.fad.seurat, dims = 1:10, do.fast = TRUE,n.neighbors = 50)

DimPlot(t.fad.seurat,label = T)
t.marker = FindAllMarkers(t.fad.seurat,min.pct = 0.15)

ggplot(t.fad.seurat@meta.data,aes(x = group, fill = SCT_snn_res.0.5))+
  geom_bar(stat = "count",position = "fill")
save(t.fad.seurat,file = "/Users/ruoqing/Projects/5xFAD_Tdepletion/r_obj/t.fad.seurat.Rdata")
rm(t.fad.seurat)

t.fad.seurat.clean = subset(t.fad.seurat,SCT_snn_res.0.5 %in% c(0,1,2,3,5,6))
t.fad.seurat.clean@meta.data$tmp = as.numeric(t.fad.seurat.clean@meta.data$SCT_snn_res.0.5)

t.fad.seurat.clean@meta.data$celltype = "Cd8_T"
t.fad.seurat.clean@meta.data$celltype[t.fad.seurat.clean@meta.data$tmp == 3] = "Cd4_T"
t.fad.seurat.clean@meta.data$celltype[t.fad.seurat.clean@meta.data$tmp == 7] = "Cxcl10"
t.fad.seurat.clean@meta.data$celltype[t.fad.seurat.clean@meta.data$tmp == 6] = "Ccr7"
t.fad.seurat.clean@meta.data$celltype[t.fad.seurat.clean@meta.data$tmp == 4] = "NK"

t.fad.seurat.clean@meta.data$cell = rownames(t.fad.seurat.clean@meta.data)
Idents(t.fad.seurat.clean) = "celltype"
t.fad.seurat.clean@meta.data$cell = rownames(t.fad.seurat.clean@meta.data$cell)
p1 = DimPlot(subset(t.fad.seurat.clean,celltype == "Cd4_T"))
p1 = DimPlot(t.fad.seurat.clean)
t.fad.seurat.clean = subset(t.fad.seurat.clean, cell %in% setdiff(t.fad.seurat.clean@meta.data$cell,CellSelector(p1)))

t.marker = FindAllMarkers(t.fad.seurat.clean,min.pct = 0.15)

tmp = randomColor(2)
t.color = tmp
names(t.color) = unique(t.fad.seurat.clean@meta.data$celltype)

tmp = cbind(t.fad.seurat.clean@meta.data,t.fad.seurat.clean@reductions$umap@cell.embeddings)

p1 = ggplot(tmp,aes(x = UMAP_1,y = UMAP_2,color = celltype))+
  geom_point(size = 5)+
  theme_dr()+
  scale_color_manual(values = as.character(t.color),breaks = names(t.color))+
  theme(panel.grid = element_blank())
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/umap_l.pdf",p1,width = 5,height = 4,dpi = 300)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/umap_l.png",p1,width = 5,height = 4,dpi = 300)


pdf("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/cd8a_exp.pdf",width = 6,height = 4)
VlnPlot(t.fad.seurat.clean,features = c("Ccr7"),split.by = "group")+
  scale_fill_manual(values = c("#4f22a3","#f2d38c"))
dev.off()

ggplot(t.fad.seurat.clean@meta.data,aes(x = group , fill = celltype))+
  geom_bar(stat = "count",position = "fill")


tmp = matrix(NA,ncol = length(unique(t.fad.seurat.clean@meta.data$celltype)),nrow = length(unique(t.fad.seurat.clean@meta.data$orig.ident)))
colnames(tmp) = unique(t.fad.seurat.clean@meta.data$celltype)
rownames(tmp) = unique(unique(t.fad.seurat.clean@meta.data$orig.ident))

for (i in colnames(tmp)){
  for (j in rownames(tmp)){
    tmp[j,i] = nrow(t.fad.seurat.clean@meta.data[t.fad.seurat.clean@meta.data$celltype == i & t.fad.seurat.clean@meta.data$orig.ident == j,])/nrow(t.fad.seurat.clean@meta.data[t.fad.seurat.clean@meta.data$orig.ident == j,])
  }
}
tmp = as.data.frame(tmp)
tmp$group = c(rep("Control",3),rep("T_depletion",3))

tmp = reshape2::melt(tmp)

p1 = ggplot(tmp,aes(x = variable,y = value,fill = group))+
  geom_boxplot(color = "grey66")+
  scale_fill_manual(values = c("#4f22a3","#f2d38c"))+
  theme_classic()+xlab("")+ylab("")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/clu_fraction.pdf",p1, width = 6,height = 4,dpi = 300)

#####
##marker plot of 
top10 <- t.marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

use.gene = top10$gene
use.gene[use.gene == "St6galnac3"] = "Cd4"
use.gene[use.gene == "Crlf3"] = "Ccr7"

tmp = AverageExpression(t.fad.seurat.clean,assays = "SCT",features = use.gene,return.seurat = T,slot = "data")
DoHeatmap(tmp,features = use.gene,slot = "data")


tmp = tmp@assays$SCT@scale.data

bk = seq(-2,2,0.1)
p1 = pheatmap(tmp,breaks = bk,
              color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)),
              cluster_rows = F,cluster_cols = F,border_color = NA)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/heat_marker.pdf",p1,width = 2,height = 5,dpi = 300)

tmp = t.fad.seurat.clean@assays$SCT@data[use.gene,rownames(t.fad.seurat.clean@meta.data)] 
tmp = t(as.matrix(tmp))
tmp = cbind(t.fad.seurat.clean@meta.data,tmp)

tmp = tmp[,c("celltype",use.gene)]
tmp = reshape2::melt(tmp)
tmp$celltype = factor(tmp$celltype,levels = c("Cd8_T","Cd4_T","Cxcl10","Ccr7","NK"),ordered = T)
p1 = ggplot(tmp,aes(x = celltype, y = value,fill = celltype))+
  geom_violin(linetype="blank",scale = "width")+
  facet_grid(variable~. ,scales = 'free')+
  scale_fill_manual(values = as.character(t.color),breaks = names(t.color))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/Vln_marker.pdf",p1,width = 2,height = 6,dpi = 300)

tmp = t.fad.seurat.clean@assays$SCT@data[use.gene,rownames(t.fad.seurat.clean@meta.data)] 
tmp = t(as.matrix(tmp))
tmp = cbind(t.fad.seurat.clean@meta.data,tmp)

tmp = tmp[,c("celltype","group",use.gene)]
tmp = reshape2::melt(tmp)
tmp$celltype = factor(tmp$celltype,levels = c("Cd8_T","Cd4_T","Cxcl10","Ccr7","NK"),ordered = T)

p1 = ggplot(tmp,aes(x = celltype, y = value,fill = group))+
  geom_split_violin(linetype="blank",scale = "width")+
  facet_grid(variable~. ,scales = 'free')+
  scale_fill_manual(values = c("#4f22a3","#f2d38c"))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/Vln_marker_group.pdf",p1,width =2.5,height = 6,dpi = 300)


tmp = CreateSeuratObject(t.fad.seurat.clean@assays$SCT@counts,meta.data = t.fad.seurat.clean@meta.data)
tmp <- NormalizeData(tmp, normalization.method =  "LogNormalize")

tmp@meta.data$pair = paste(tmp@meta.data$group,tmp@meta.data$celltype,sep = "_")

Idents(tmp) = "pair"
for (i in unique(tmp@meta.data$celltype)){
  print(i)
  tmp.m = FindMarkers(tmp,ident.1 = paste("Control",i,sep = "_"),ident.2 = paste("T_depletion",i,sep = "_"),min.pct = 0,test.use = "MAST",logfc.threshold = 0)
  tmp.m$diff = tmp.m$pct.1 - tmp.m$pct.2
  tmp.m = na.omit(tmp.m)
  tmp.m = tmp.m[order(tmp.m$avg_log2FC,decreasing = T),]
  tmp.m$gene = rownames(tmp.m)
  tmp.m = tmp.m[tmp.m$gene != "Ttr",]
  rownames(tmp.m) = NULL
  tmp.m$label = ""
  tmp.m$label[1:5] = tmp.m$gene[1:5]
  tmp.m$label[(nrow(tmp.m)-4):nrow(tmp.m)] = tmp.m$gene[(nrow(tmp.m)-4):nrow(tmp.m)]
  write.csv(tmp.m,paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/wtvsko/",i,".csv"))
  thexlim = max(abs(tmp.m$diff))
  theylim = max(abs(tmp.m$avg_log2FC))
  
  p1 = ggplot(tmp.m,aes(x = diff, y = avg_log2FC))+
    geom_point(color = "grey66",size = 1)+
    geom_text_repel(data = tmp.m[tmp.m$label != "",],aes(label = label))+
    geom_vline(xintercept = 0.0,linetype = 2)+
    geom_hline(yintercept = 0,linetype = 2)+
    theme_bw()+theme_classic()+
    xlim(-theylim,theylim)+ylim(-theylim,theylim)+
    xlab("")
  ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/wtvsko/",i,"_wtvsko.png"),p1,width = 6,height = 4,dpi = 300)
  ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/wtvsko/",i,"_wtvsko.png"),p1,width = 6,height = 4,dpi = 300)
  
  if(theylim < 0.58){
    theylim = 0.6
  }
  rownames(tmp.m) = tmp.m$gene
  tmp.m$label = ""
  tmp.m$group = "non"
  for (s in tmp.m$gene) {
    if (tmp.m[s,"p_val"] <= 0.05 & tmp.m[s,"avg_log2FC"] >= 0.58) {
      tmp.m[s,"label"] = s
      tmp.m[s,"group"] = "up"
    }
    if (tmp.m[s,"p_val"] <= 0.05 & tmp.m[s,"avg_log2FC"] <= -0.58) {
      tmp.m[s,"label"] = s
      tmp.m[s,"group"] = "down"
    }
    
  }
  maxp = max(-log10(tmp.m$p_val))
  p1 = ggplot(tmp.m,aes(x = avg_log2FC, y = -log10(p_val)))+
    geom_point(size = 1,aes(color = group))+
    scale_color_manual(breaks = c("non","up","down"),values = c("grey66","red","navy"))+
    geom_text_repel(data = tmp.m[tmp.m$label != "",],aes(label = label))+
    geom_vline(xintercept = 0.58,linetype = 2,color = "grey88")+
    geom_vline(xintercept = -0.58,linetype = 2,color = "grey88")+
    geom_hline(yintercept = -log10(0.05),linetype = 2,color = "grey88")+
    theme_bw()+theme_classic()+
    annotate(geom="text", x= theylim/2, y=maxp+1, label="Control",  fontface="bold",colour='red', size=5)+
    annotate(geom="text", x= -theylim/2, y=maxp+1, label="T_depletion",  fontface="bold",colour='navy', size=5)+
    theme(legend.position = 'none')+
    xlim(-theylim,theylim)
  ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/wtvsko/",i,"_volcano.png"),p1,width = 6,height = 4,dpi = 300)
  ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/T/wtvsko/",i,"_volcano.pdf"),p1,width = 6,height = 4,dpi = 300)
  
}




#####
##Microglia
microglia.fad.seurat = subset(fad.seurat,SCT_snn_res.0.1 %in% c(1,0,4))

SCTransform(microglia.fad.seurat, 
            vars.to.regress = "percent.mt", 
            verbose = T,
            variable.features.rv.th = 1.4,
            return.only.var.genes = F,
            variable.features.n = NULL)
microglia.fad.seurat@assays$SCT@var.features = setdiff(microglia.fad.seurat@assays$SCT@var.features,str_subset(microglia.fad.seurat@assays$SCT@var.features, pattern ='^mt-'))

DefaultAssay(microglia.fad.seurat) <- "SCT"

microglia.fad.seurat <- RunPCA(object = microglia.fad.seurat, npcs = 50,verbose = T)
pdf("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/elbow.pdf",width = 5,height = 5)
ElbowPlot(microglia.fad.seurat,ndims = 50)
dev.off()

microglia.fad.seurat <- FindNeighbors(microglia.fad.seurat, dims = 1:20)
microglia.fad.seurat <- FindClusters(microglia.fad.seurat, resolution = 0.5)
microglia.fad.seurat <- RunUMAP(object = microglia.fad.seurat, dims = 1:20, do.fast = TRUE,n.neighbors = 50)

DimPlot(microglia.fad.seurat,label = T,split.by = "group")
FeaturePlot(microglia.fad.seurat,features = c("P2ry12","Apoe"))
VlnPlot(microglia.fad.seurat,features = c("P2ry12","Apoe"))

microglia.marker = FindAllMarkers(microglia.fad.seurat,min.pct = 0.15)
microglia.marker2 = FindAllMarkers(microglia.fad.seurat,min.pct = 0.15)
Idents(microglia.fad.seurat) = "SCT_snn_res.0.1"

microglia.fad.seurat@meta.data$tmp = as.numeric(microglia.fad.seurat@meta.data$SCT_snn_res.0.5)
Idents(microglia.fad.seurat) = "celltype"
DimPlot(microglia.fad.seurat,label = T)

tmp = FindMarkers(microglia.fad.seurat,ident.1 = 13,ident.2 = 4)

microglia.fad.seurat@meta.data$celltype = "Homeostatic"
microglia.fad.seurat@meta.data$celltype[microglia.fad.seurat@meta.data$tmp == 10] = "IRMs"
microglia.fad.seurat@meta.data$celltype[microglia.fad.seurat@meta.data$tmp == 12] = "Chemokines"
microglia.fad.seurat@meta.data$celltype[microglia.fad.seurat@meta.data$tmp == 8] = "ADM"
microglia.fad.seurat@meta.data$celltype[microglia.fad.seurat@meta.data$tmp == 11] = "Phagocytosis"
microglia.fad.seurat@meta.data$celltype[microglia.fad.seurat@meta.data$tmp == 7] = "MHCII"
microglia.fad.seurat@meta.data$celltype[microglia.fad.seurat@meta.data$tmp %in% c(1,5)] = "Activated"
microglia.fad.seurat@meta.data$celltype[microglia.fad.seurat@meta.data$tmp %in% c(13,4)] = "Ribosomal"

microglia.fad.seurat.clean = subset(microglia.fad.seurat,tmp %in% c(1:13))
microglia.fad.seurat.clean@meta.data$cell = rownames(microglia.fad.seurat.clean@meta.data)

p1 = DimPlot(subset(microglia.fad.seurat.clean,celltype == "Ribosomal"))
microglia.fad.seurat.clean = subset(microglia.fad.seurat.clean,cell %in% setdiff(rownames(microglia.fad.seurat.clean@meta.data),CellSelector(p1)))

Idents(microglia.fad.seurat.clean) = "celltype"
DimPlot(microglia.fad.seurat.clean,label = T,cols = microglia.color)+
  theme_dr()+
  theme(panel.grid = element_blank())

tmp = randomcoloR::distinctColorPalette(8)
scales::show_col(tmp)

microglia.color = tmp
names(microglia.color) = unique(microglia.fad.seurat@meta.data$celltype)

tmp = cbind(microglia.fad.seurat.clean@meta.data,microglia.fad.seurat.clean@reductions$umap@cell.embeddings)

Idents(microglia.fad.seurat.clean) = "celltype"
p1 = DimPlot(microglia.fad.seurat.clean,label = T,cols = microglia.color)+
  theme_dr()+
  theme(panel.grid = element_blank())
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/umap_l.pdf",p1,width = 10,height = 6,dpi = 300)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/umap_l.png",p1,width = 10,height = 6,dpi = 300)


tmp = matrix(NA,ncol = length(unique(microglia.fad.seurat.clean@meta.data$celltype)),nrow = length(unique(microglia.fad.seurat.clean@meta.data$orig.ident)))
colnames(tmp) = unique(microglia.fad.seurat.clean@meta.data$celltype)
rownames(tmp) = unique(unique(microglia.fad.seurat.clean@meta.data$orig.ident))

for (i in colnames(tmp)){
  for (j in rownames(tmp)){
    tmp[j,i] = nrow(microglia.fad.seurat.clean@meta.data[microglia.fad.seurat.clean@meta.data$celltype == i & microglia.fad.seurat.clean@meta.data$orig.ident == j,])/nrow(microglia.fad.seurat.clean@meta.data[microglia.fad.seurat.clean@meta.data$orig.ident == j,])
  }
}
tmp = as.data.frame(tmp)
tmp$group = c(rep("Control",3),rep("T_depletion",3))

tmp = reshape2::melt(tmp)
tmp$variable = factor(tmp$variable,levels = c("Homeostatic","Activated","MHCII","Ribosomal","IRMs","ADM","Phagocytosis","Chemokines"),ordered = T)
p1 = ggplot(tmp,aes(x = variable,y = value,fill = group))+
  geom_boxplot(color = "grey66",size = 0.2)+
  scale_fill_manual(values = c("#4f22a3","#f2d38c"))+
  theme_classic()+
  #stat_compare_means(aes(x=variable,y=value,group=group),method = "t.test")+
  xlab("")+ylab("")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/fraction.pdf",p1,width =8,height = 4,dpi = 300)

microglia.fad.seurat.clean@meta.data$celltype = factor(microglia.fad.seurat.clean@meta.data$celltype,levels = c("Homeostatic","Activated","MHCII","Ribosomal","IRMs","ADM","Phagocytosis","Chemokines"),ordered = T)
Idents(microglia.fad.seurat.clean) = "celltype"
microglia.marker = FindAllMarkers(microglia.fad.seurat.clean,min.pct = 0.15)

top10 <- microglia.marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
use.gene = top10$gene

tmp = AverageExpression(microglia.fad.seurat.clean,assays = "SCT",features = use.gene,return.seurat = T)
DoHeatmap(tmp,features = use.gene)

tmp = tmp@assays$SCT@scale.data

bk = seq(-2,2,0.1)
p1 = pheatmap(tmp,breaks = bk,
              color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)),
              cluster_rows = F,cluster_cols = F,border_color = NA)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/heat_marker.pdf",p1,width = 3,height = 6,dpi = 300)

tmp = microglia.fad.seurat.clean@assays$SCT@data[unique(use.gene),rownames(microglia.fad.seurat.clean@meta.data)] 
tmp = t(as.matrix(tmp))
tmp = cbind(microglia.fad.seurat.clean@meta.data,tmp)

tmp = tmp[,c("celltype",unique(use.gene))]
tmp = reshape2::melt(tmp)
tmp$celltype = factor(tmp$celltype,levels = c("Homeostatic","Activated","MHCII","Ribosomal","IRMs","ADM","Phagocytosis","Chemokines"),ordered = T)
p1 = ggplot(tmp,aes(x = celltype, y = value,fill = celltype))+
  geom_violin(linetype="blank",scale = "width")+
  facet_grid(variable~. ,scales = 'free')+
  scale_fill_manual(values = as.character(microglia.color),breaks = names(microglia.color))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/Vln_marker.pdf",p1,width = 2.5,height = 8,dpi = 300)

tmp = microglia.fad.seurat.clean@assays$SCT@data[unique(use.gene),rownames(microglia.fad.seurat.clean@meta.data)] 
tmp = t(as.matrix(tmp))
tmp = cbind(microglia.fad.seurat.clean@meta.data,tmp)

tmp = tmp[,c("celltype","group",unique(use.gene))]
tmp = reshape2::melt(tmp)
tmp$celltype = factor(tmp$celltype,levels = c("Homeostatic","Activated","MHCII","Ribosomal","IRMs","ADM","Phagocytosis","Chemokines"),ordered = T)

p1 = ggplot(tmp,aes(x = celltype, y = value,fill = group))+
  geom_split_violin(linetype="blank",scale = "width")+
  facet_grid(variable~. ,scales = 'free')+
  scale_fill_manual(values = c("#4f22a3","#f2d38c"))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank())
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/Vln_marker_group.pdf",p1,width = 4,height = 7,dpi = 300)

p1 = VlnPlot(microglia.fad.seurat.clean,features = "Apoe",cols = microglia.color,pt.size = 0)+
  theme(legend.position = 'none')
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/apoe_exp.pdf",p1,width = 4,height = 3,dpi = 300)

tmp1 = microglia.fad.seurat.clean@meta.data[microglia.fad.seurat.clean@meta.data$celltype %in% c("ADM","Phagocytosis"),]
tmp = matrix(NA,ncol = length(unique(tmp1$celltype)),nrow = length(unique(tmp1$orig.ident)))
colnames(tmp) = unique(tmp1$celltype)
rownames(tmp) = unique(unique(tmp1$orig.ident))

for (i in colnames(tmp)){
  for (j in rownames(tmp)){
    tmp[j,i] = nrow(tmp1[tmp1$celltype == i & tmp1$orig.ident == j,])/nrow(tmp1[tmp1$orig.ident == j,])
  }
}
tmp = as.data.frame(tmp)
tmp$group = c(rep("Control",3),rep("T_depletion",3))

tmp = reshape2::melt(tmp)

p1 = ggplot(tmp,aes(x = variable,y = value,fill = group))+
  geom_boxplot(color = "grey66")+
  scale_fill_manual(values = c("#4f22a3","#f2d38c"))+
  theme_classic()+xlab("")+ylab("")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/ADM/fraction_adm.pdf",p1,width = 5,height = 3,dpi = 300)

tmp = CreateSeuratObject(microglia.fad.seurat.clean@assays$SCT@counts,meta.data = microglia.fad.seurat.clean@meta.data)
tmp = subset(tmp,celltype %in% c("ADM","Phagocytosis"))
tmp <- NormalizeData(tmp, normalization.method =  "LogNormalize")
Idents(tmp) = "celltype"

tmp.m = FindMarkers(tmp,ident.1 = "ADM",ident.2 = "Phagocytosis",test.use = "DESeq2", max.cells.per.ident = 30)
tmp.m$diff = tmp.m$pct.1 - tmp.m$pct.2
tmp.m = na.omit(tmp.m)
tmp.m = tmp.m[order(tmp.m$avg_log2FC,decreasing = T),]
tmp.m$gene = rownames(tmp.m)
rownames(tmp.m) = NULL
tmp.m$label = ""
tmp.m$label[1:5] = tmp.m$gene[1:5]
tmp.m$label[(nrow(tmp.m)-4):nrow(tmp.m)] = tmp.m$gene[(nrow(tmp.m)-4):nrow(tmp.m)]
thexlim = max(abs(tmp.m$diff))
theylim = max(abs(tmp.m$avg_log2FC))
if(theylim < 0.58){
  theylim = 0.6
}
rownames(tmp.m) = tmp.m$gene
tmp.m$label = ""
tmp.m$group = "non"
for (s in tmp.m$gene) {
  if (tmp.m[s,"p_val"] <= 0.05 & tmp.m[s,"avg_log2FC"] >= 0.58) {
    tmp.m[s,"label"] = s
    tmp.m[s,"group"] = "up"
  }
  if (tmp.m[s,"p_val"] <= 0.05 & tmp.m[s,"avg_log2FC"] <= -0.58) {
    tmp.m[s,"label"] = s
    tmp.m[s,"group"] = "down"
  }
  
}
maxp = max(-log10(tmp.m$p_val))
p1 = ggplot(tmp.m,aes(x = avg_log2FC, y = -log10(p_val)))+
  geom_point(size = 1,aes(color = group))+
  scale_color_manual(breaks = c("non","up","down"),values = c("grey66","red","navy"))+
  geom_text_repel(data = tmp.m[tmp.m$label != "",],aes(label = label))+
  geom_vline(xintercept = 0.58,linetype = 2,color = "grey88")+
  geom_vline(xintercept = -0.58,linetype = 2,color = "grey88")+
  geom_hline(yintercept = -log10(0.05),linetype = 2,color = "grey88")+
  theme_bw()+theme_classic()+
  annotate(geom="text", x=1, y=maxp, label="ADM",  fontface="bold",colour='red', size=5)+
  annotate(geom="text", x= -3, y=maxp, label="Phagocytosis",  fontface="bold",colour='navy', size=5)+
  theme(legend.position = 'none')
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/ADM/adm1vs2.pdf",p1,width = 5,height = 5,dpi = 300)
#############

astrocyte.fad.seurat = subset(fad.seurat.clean,annotation == "Astrocytes")

SCTransform(astrocyte.fad.seurat, 
            vars.to.regress = "percent.mt", 
            verbose = T,
            variable.features.rv.th = 1.4,
            return.only.var.genes = F,
            variable.features.n = NULL)
astrocyte.fad.seurat@assays$SCT@var.features = setdiff(astrocyte.fad.seurat@assays$SCT@var.features,str_subset(astrocyte.fad.seurat@assays$SCT@var.features, pattern ='^mt-'))

DefaultAssay(astrocyte.fad.seurat) <- "SCT"

astrocyte.fad.seurat <- RunPCA(object = astrocyte.fad.seurat, npcs = 50,verbose = T)
pdf("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/elbow.pdf",width = 5,height = 5)
ElbowPlot(astrocyte.fad.seurat,ndims = 50)
dev.off()

astrocyte.fad.seurat <- FindNeighbors(astrocyte.fad.seurat, dims = 1:10)
astrocyte.fad.seurat <- FindClusters(astrocyte.fad.seurat, resolution = 0.1)
astrocyte.fad.seurat <- RunUMAP(object = astrocyte.fad.seurat, dims = 1:10, do.fast = TRUE,n.neighbors = 50)

DimPlot(astrocyte.fad.seurat,label = T,split.by = "group")

##"Apoe","Cst7", "Lpl","Spp1"
##"P2ry12","Tmem119"
##"Ifit1","Iigp1","Ifitm3"
##"H2-Eb2","H2-Eb1","H2-Aa","H2-Ab1","H2-Ob","H2-DMb1","H2-DMb2","H2-DMa","H2-Oa"
for (i in c("Isg15","Ifit2")){
p1 =FeaturePlot(microglia.fad.seurat.clean,
            features = i,
            cols = brewer.pal(5, "OrRd"),raster=FALSE)+
  theme_dr()+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_blank(),
        panel.background=element_blank(),
        axis.ticks = element_blank())
ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/feature/",i,".png"),p1,width = 6,height = 5,dpi = 300)
ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/feature/",i,".pdf"),p1,width = 6,height = 5,dpi = 300)
}

microglia.fad.seurat.clean@meta.data$celltype2 = as.character(microglia.fad.seurat.clean@meta.data$celltype)
microglia.fad.seurat.clean@meta.data$celltype2[microglia.fad.seurat.clean@meta.data$celltype == "Ribosomal"] = "DAM2"
microglia.fad.seurat.clean@meta.data$celltype2[microglia.fad.seurat.clean@meta.data$celltype == "Activated"] = "DAM1"
microglia.fad.seurat.clean@meta.data$celltype2[microglia.fad.seurat.clean@meta.data$celltype == "ADM"] = "DAM1.2"
microglia.fad.seurat.clean@meta.data$celltype2[microglia.fad.seurat.clean@meta.data$celltype == "Chemokines"] = "DAM1.1"
microglia.fad.seurat.clean@meta.data$celltype2[microglia.fad.seurat.clean@meta.data$celltype == "Phagocytosis"] = "DAM3"
microglia.fad.seurat.clean@meta.data$celltype2[microglia.fad.seurat.clean@meta.data$celltype == "MHCII"] = "DAM4"

microglia.color2 = microglia.color
names(microglia.color2) = unique(microglia.fad.seurat.clean@meta.data$celltype2)
Idents(microglia.fad.seurat.clean) = "celltype2"
p1 = DimPlot(microglia.fad.seurat.clean,label = T,cols = microglia.color2)+
  theme_dr()+
  theme(panel.grid = element_blank())
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/umap2_l.pdf",p1,width = 10,height = 6,dpi = 300)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/umap2_l.png",p1,width = 10,height = 6,dpi = 300)

tmp = matrix(NA,ncol = length(unique(microglia.fad.seurat.clean@meta.data$celltype2)),nrow = length(unique(microglia.fad.seurat.clean@meta.data$orig.ident)))
colnames(tmp) = unique(microglia.fad.seurat.clean@meta.data$celltype2)
rownames(tmp) = unique(unique(microglia.fad.seurat.clean@meta.data$orig.ident))

for (i in colnames(tmp)){
  for (j in rownames(tmp)){
    tmp[j,i] = nrow(microglia.fad.seurat.clean@meta.data[microglia.fad.seurat.clean@meta.data$celltype2 == i & microglia.fad.seurat.clean@meta.data$orig.ident == j,])/nrow(microglia.fad.seurat.clean@meta.data[microglia.fad.seurat.clean@meta.data$orig.ident == j,])
  }
}
tmp = as.data.frame(tmp)
tmp$group = c(rep("Control",3),rep("T_depletion",3))

tmp = reshape2::melt(tmp)
tmp$variable = factor(tmp$variable,levels = c("Homeostatic","IRMs","DAM1","DAM1.1","DAM1.2","DAM2","DAM3","DAM4"),ordered = T)
p1 = ggplot(tmp,aes(x = variable,y = value,fill = group))+
  geom_boxplot(color = "grey66",size = 0.2)+
  scale_fill_manual(values = c("#4f22a3","#f2d38c"))+
  theme_classic()+
  stat_compare_means(aes(x=variable,y=value,group=group),method = )+
  xlab("")+ylab("")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/fraction2.pdf",p1,width =6,height = 2,dpi = 300)

microglia.fad.seurat.clean@meta.data$celltype2 = factor(microglia.fad.seurat.clean@meta.data$celltype2,levels = c("Homeostatic","IRMs","DAM1","DAM1.1","DAM1.2","DAM2","DAM3","DAM4"),ordered = T)
Idents(microglia.fad.seurat.clean) = "celltype2"
microglia.marker = FindAllMarkers(microglia.fad.seurat.clean,min.pct = 0.15)

top10 <- microglia.marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
use.gene = top10$gene

tmp = AverageExpression(microglia.fad.seurat.clean,assays = "SCT",features = use.gene,return.seurat = T)
DoHeatmap(tmp,features = use.gene)

tmp = tmp@assays$SCT@scale.data

bk = seq(-2,2,0.1)
p1 = pheatmap(tmp,breaks = bk,
              color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)),
              cluster_rows = F,cluster_cols = F,border_color = NA)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/heat_marker2.pdf",p1,width = 3,height = 6,dpi = 300)

tmp.dam = microglia.fad.seurat.clean@meta.data
tmp.dam = tmp.dam[tmp.dam$celltype2 %in% c("DAM1","DAM1.1","DAM1.2","DAM2","DAM3","DAM4"),]
tmp = matrix(NA,ncol = length(unique(tmp.dam$celltype2)),nrow = length(unique(tmp.dam$orig.ident)))
colnames(tmp) = unique(tmp.dam$celltype2)
rownames(tmp) = unique(unique(tmp.dam$orig.ident))

for (i in colnames(tmp)){
  for (j in rownames(tmp)){
    tmp[j,i] = nrow(tmp.dam[tmp.dam$celltype2 == i & tmp.dam$orig.ident == j,])/nrow(tmp.dam[tmp.dam$orig.ident == j,])
  }
}
tmp = as.data.frame(tmp)
tmp$group = c(rep("Control",3),rep("T_depletion",3))

tmp = reshape2::melt(tmp)
tmp$variable = factor(tmp$variable,levels = c("DAM1","DAM1.1","DAM1.2","DAM2","DAM3","DAM4"),ordered = T)
p1 = ggplot(tmp,aes(x = variable,y = value,fill = group))+
  geom_boxplot(color = "grey66")+
  scale_fill_manual(values = c("#4f22a3","#f2d38c"))+
  theme_classic()+
  stat_compare_means(aes(x=variable,y=value,group=group),method = "t.test")+
  xlab("")+ylab("")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/fraction_dam.pdf",p1,width =6,height = 4,dpi = 300)


tmp = CreateSeuratObject(microglia.fad.seurat.clean@assays$SCT@counts,meta.data = microglia.fad.seurat.clean@meta.data)
tmp <- NormalizeData(tmp, normalization.method =  "LogNormalize")

tmp@meta.data$pair = paste(tmp@meta.data$group,tmp@meta.data$celltype2,sep = "_")

Idents(tmp) = "pair"
for (i in unique(tmp@meta.data$celltype2)){
  print(i)
    tmp.m = FindMarkers(tmp,ident.1 = paste("Control",i,sep = "_"),ident.2 = paste("T_depletion",i,sep = "_"),min.pct = 0,test.use = "MAST",logfc.threshold = 0)
  tmp.m$diff = tmp.m$pct.1 - tmp.m$pct.2
  tmp.m = na.omit(tmp.m)
  tmp.m = tmp.m[order(tmp.m$avg_log2FC,decreasing = T),]
  tmp.m$gene = rownames(tmp.m)
  tmp.m = tmp.m[tmp.m$gene != "Ttr",]
  rownames(tmp.m) = NULL
  tmp.m$label = ""
  tmp.m$label[1:5] = tmp.m$gene[1:5]
  tmp.m$label[(nrow(tmp.m)-4):nrow(tmp.m)] = tmp.m$gene[(nrow(tmp.m)-4):nrow(tmp.m)]
  write.csv(tmp.m,paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/wtvsko/",i,".csv"))
  thexlim = max(abs(tmp.m$diff))
  theylim = max(abs(tmp.m$avg_log2FC))
  
  p1 = ggplot(tmp.m,aes(x = diff, y = avg_log2FC))+
    geom_point(color = "grey66",size = 1)+
    geom_text_repel(data = tmp.m[tmp.m$label != "",],aes(label = label))+
    geom_vline(xintercept = 0.0,linetype = 2)+
    geom_hline(yintercept = 0,linetype = 2)+
    theme_bw()+theme_classic()+
    xlim(-theylim,theylim)+ylim(-theylim,theylim)+
    xlab("")
  ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/wtvsko/",i,"_wtvsko.png"),p1,width = 6,height = 4,dpi = 300)
  ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/wtvsko/",i,"_wtvsko.png"),p1,width = 6,height = 4,dpi = 300)
  
  if(theylim < 0.58){
    theylim = 0.6
  }
  rownames(tmp.m) = tmp.m$gene
  tmp.m$label = ""
  tmp.m$group = "non"
  for (s in tmp.m$gene) {
    if (tmp.m[s,"p_val"] <= 0.05 & tmp.m[s,"avg_log2FC"] >= 0.58) {
      tmp.m[s,"label"] = s
      tmp.m[s,"group"] = "up"
    }
    if (tmp.m[s,"p_val"] <= 0.05 & tmp.m[s,"avg_log2FC"] <= -0.58) {
      tmp.m[s,"label"] = s
      tmp.m[s,"group"] = "down"
    }
    
  }
  maxp = max(-log10(tmp.m$p_val))
  p1 = ggplot(tmp.m,aes(x = avg_log2FC, y = -log10(p_val)))+
    geom_point(size = 1,aes(color = group))+
    scale_color_manual(breaks = c("non","up","down"),values = c("grey66","red","navy"))+
    geom_text_repel(data = tmp.m[tmp.m$label != "",],aes(label = label))+
    geom_vline(xintercept = 0.58,linetype = 2,color = "grey88")+
    geom_vline(xintercept = -0.58,linetype = 2,color = "grey88")+
    geom_hline(yintercept = -log10(0.05),linetype = 2,color = "grey88")+
    theme_bw()+theme_classic()+
    annotate(geom="text", x=1, y=maxp, label="Control",  fontface="bold",colour='red', size=5)+
    annotate(geom="text", x= -1, y=maxp, label="T_depletion",  fontface="bold",colour='navy', size=5)+
    theme(legend.position = 'none')+
    xlim(-theylim,theylim)
  ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/wtvsko/",i,"_volcano.png"),p1,width = 6,height = 4,dpi = 300)
  ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/wtvsko/",i,"_volcano.png"),p1,width = 6,height = 4,dpi = 300)
  
}


astrocyte.fad.seurat.clean = astrocyte.fad.seurat
Idents(astrocyte.fad.seurat.clean)="SCT_snn_res.0.1"
p1 = DimPlot(subset(astrocyte.fad.seurat.clean,SCT_snn_res.0.1 == 2))
astrocyte.fad.seurat.clean = subset(astrocyte.fad.seurat.clean, cell %in% setdiff(astrocyte.fad.seurat.clean@meta.data$cell,CellSelector(p1)))

DimPlot(astrocyte.fad.seurat.clean)

astrocyte.marker = FindAllMarkers(astrocyte.fad.seurat.clean,min.pct = 0.15)
FeaturePlot(astrocyte.fad.seurat.clean,features = "Aqp4")

save(astrocyte.fad.seurat,file = "/Users/ruoqing/Projects/5xFAD_Tdepletion/r_obj/astrocyte.fad.seurat.Rdata")
rm(astrocyte.fad.seurat)

top10 <- astrocyte.marker %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC)
use.gene = top10$gene

tmp = AverageExpression(astrocyte.fad.seurat.clean,assays = "SCT",features = use.gene,return.seurat = T)
DoHeatmap(tmp,features = use.gene)

tmp = tmp@assays$SCT@scale.data

bk = seq(-2,2,0.1)
p1 = pheatmap(tmp,breaks = bk,
              color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)),
              cluster_rows = F,cluster_cols = F,border_color = NA)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/astrocyte/heat_marker.pdf",p1,width = 2,height = 6,dpi = 300)

tmp = matrix(NA,ncol = length(unique(astrocyte.fad.seurat.clean@meta.data$SCT_snn_res.0.1)),nrow = length(unique(astrocyte.fad.seurat.clean@meta.data$orig.ident)))
colnames(tmp) = unique(astrocyte.fad.seurat.clean@meta.data$SCT_snn_res.0.1)
rownames(tmp) = unique(unique(astrocyte.fad.seurat.clean@meta.data$orig.ident))

for (i in colnames(tmp)){
  for (j in rownames(tmp)){
    tmp[j,i] = nrow(astrocyte.fad.seurat.clean@meta.data[astrocyte.fad.seurat.clean@meta.data$SCT_snn_res.0.1 == i & astrocyte.fad.seurat.clean@meta.data$orig.ident == j,])/nrow(astrocyte.fad.seurat.clean@meta.data[astrocyte.fad.seurat.clean@meta.data$orig.ident == j,])
  }
}
tmp = as.data.frame(tmp)
tmp$group = c(rep("Control",3),rep("T_depletion",3))

tmp = reshape2::melt(tmp)
tmp$variable = factor(tmp$variable,levels = c(0,1,2),ordered = T)
p1 = ggplot(tmp,aes(x = variable,y = value,fill = group))+
  geom_boxplot(color = "grey66")+
  scale_fill_manual(values = c("#4f22a3","#f2d38c"))+
  theme_classic()+
  #stat_compare_means(aes(x=variable,y=value,group=group),method = "t.test")+
  xlab("")+ylab("")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/astrocyte/fraction.pdf",p1,width =5,height = 4,dpi = 300)

tmp = randomcoloR::distinctColorPalette(4)
scales::show_col(tmp)
astrocyte.color = tmp[2:4]
names(astrocyte.color) = unique(astrocyte.fad.seurat.clean@meta.data$SCT_snn_res.0.1)
p1 = DimPlot(astrocyte.fad.seurat.clean,cols = astrocyte.color,label = T)+
  theme_dr()+
  theme(panel.grid = element_blank(),
        legend.position = "none")
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/astrocyte/umap.png",p1,width = 9,height = 6,dpi = 300)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/astrocyte/umap.pdf",p1,width = 9,height = 6,dpi = 300)


write.table(microglia.marker,"/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/microglia_marker.csv",row.names = F)



saveRDS(microglia.fad.seurat.clean@assays$RNA@counts,file = "/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/microglia_count.rds")
saveRDS(microglia.fad.seurat.clean@meta.data,file = "/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/microglia_meta.rds")
