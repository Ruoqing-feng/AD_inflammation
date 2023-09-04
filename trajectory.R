library(monocle3)

saveRDS(intergrate2.seurat,file = "/Users/ruoqing/Projects/5xFAD_Tdepletion/r_obj/integrate_seurat.rds")
saveRDS(microglia.color3,file = "/Users/ruoqing/Projects/5xFAD_Tdepletion/r_obj/integrate_seurat_color.rds")

intergrate2.seurat = read_rds("/Users/ruoqing/Projects/5xFAD_Tdepletion/r_obj/integrate_seurat.rds")
microglia.color3 = read_rds("/Users/ruoqing/Projects/5xFAD_Tdepletion/r_obj/integrate_seurat_color.rds")
intergrate2.traj.seurat = subset(intergrate2.seurat,project != "wt" & celltype_new != "Hemo2")

#####
tmp = intergrate2.traj.seurat@assays$RNA@counts
tmp = tmp[setdiff(rownames(tmp),str_subset(rownames(tmp), pattern ='^mt-|^Gm')),]

cell_metadata = intergrate2.traj.seurat@meta.data[,c("cell","celltype_new","orig.ident","group")]
gene_annotation = as.data.frame(rownames(tmp))
colnames(gene_annotation) = "gene_short_name"
rownames(gene_annotation) = gene_annotation$gene_short_name

cds <- new_cell_data_set(tmp,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 30,method = "PCA") 

plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds, preprocess_method = "PCA",reduction_method = c("UMAP"))
plot_cells(cds,color_cells_by = "celltype_new")

cds = cluster_cells(cds, resolution=1)

cds@int_colData$reducedDims$UMAP = Embeddings(intergrate2.traj.seurat, reduction = "umap")[rownames(cds@int_colData$reducedDims$UMAP),]
plot_cells(cds,color_cells_by = "celltype_new",cell_size = 1)

tmp = colData(cds)$celltype_new
names(tmp) = colData(cds)$cell
cds@clusters$UMAP$clusters = tmp
cds <- learn_graph(cds,verbose = T)

cds = order_cells(cds)

p1 = plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_roots = T,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=5,
                trajectory_graph_color = "grey55",
                group_label_size=8,
                cell_size=1)+
  theme_dr()+
  theme(panel.grid = element_blank())
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/figures/monocle3_pseudotimo_nl.png",p1,width = 6,height = 6,dpi = 300)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/figures/monocle3_pseudotimo_l.pdf",p1,width = 6,height = 6,dpi = 300)


ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/trajectory/monocle3_pseudotimo.png",p1,width = 8,height = 6,dpi = 300)

p1 = plot_cells(cds,
                color_cells_by = "celltype_new",
                label_cell_groups=T,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=5,
                trajectory_graph_color = "grey33",
                group_label_size=5,
                cell_size=1)+
  theme_dr()+
  scale_color_manual(breaks = names(microglia.color3),values = microglia.color3)+
  theme(panel.grid = element_blank(),legend.position = 'none')
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/trajectory/monocle3_pseudotimo_cl.png",p1,width = 8,height = 6,dpi = 300)

p1 = plot_cells(cds,
                color_cells_by = "celltype_new",
                label_cell_groups=F,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                show_trajectory_graph = F,
                label_roots = F,
                graph_label_size=5,
                trajectory_graph_color = "grey33",
                group_label_size=5,
                cell_size=0.8)+
  scale_color_manual(breaks = names(microglia.color3),values = microglia.color3)+
  theme(axis.text = element_blank(),
        axis.title= element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        legend.position = 'none')
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/figures/2d_cluster.png",p1,width = 6,height = 6,dpi = 300)



tmp = choose_graph_segments(cds,return_list = T)
traj.cells = tmp
tmp = cds[,tmp$cells]


modulated_genes <- graph_test(tmp, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.3 & status == "OK"))

library(ClusterGVis)
mat <- pre_pseudotime_matrix(cds_obj = tmp,
                             gene_list = setdiff(c(genes,"Cxcl10"),c("")))

ck <- clusterData(exp = mat,
                  cluster.method = "kmean",
                  cluster.num = 8)
library(org.Mm.eg.db)


s = ck$long.res
s$tmp = as.character(s$cluster_name)
s$tmp[s$tmp %in% c("cluster 7 (4)","cluster 8 (14)")] = "cluster 7 (18)" 
levels(s$tmp) = c("cluster 1 (29)","cluster 2 (45)","cluster 3 (7)" ,"cluster 4 (4)","cluster 5 (8)","cluster 6 (37)","cluster 7 (18)")
s$tmp = factor(s$tmp,levels = c("cluster 1 (29)","cluster 2 (45)","cluster 3 (7)" ,"cluster 4 (4)","cluster 5 (8)","cluster 6 (37)","cluster 7 (18)"),ordered = T)
s = s[,c(1,2,3,4,6)]
colnames(s)[5] = "cluster_name"
s$cluster[s$cluster %in% c(7,8)] = 7
ck$long.res = s
ck$wide.res$cluster[ck$wide.res$cluster == 8] = 7

enrich <- enrichCluster(object = ck,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        organism = "mmu",
                        pvalueCutoff = 0.5,
                        topn = 30)
enrich1 = enrich[c("GO:0016358","GO:0106027","GO:0050808...151",
                  "GO:0050808...122","GO:0021549","GO:0022037",
                  "GO:0051607","GO:0019221","GO:0034341...205","GO:0045089",
                  "GO:0007596","GO:0034109","GO:0032640",
                  "GO:0006898...1","GO:1900221","GO:0014009","GO:1990000","GO:0034381","GO:0045834",
                  "GO:0034113","GO:0030643","GO:2000833",
                  "GO:0019886","GO:0002478","GO:0002399"),]
add_pseudo = as.data.frame(as.numeric(as.numeric(pseudotime(tmp))))
colnames(add_pseudo) = "pseudo"
add_pseudo = add_pseudo[order(add_pseudo$pseudo),]

ck$pseudotime = add_pseudo
ck$type = "scRNAdata"
ck$geneMode = "all"
ck$geneType = "non-branched"
pdf("/Users/ruoqing/Projects/5xFAD_Tdepletion/figures/heat_monocle3.pdf",width = 12,height = 10)
p = visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           show_row_dend = F,
           markGenes.side = "left",
           line.side = "left",
           cluster.order = rev(c(4,3,1,2,7,5,6)),
           markGenes = c("Apoe","P2ry12","Cd74","Isg15","Stat1","Tmem119","Spp1","Gpnmb","Cxcl10","Rpl38","Csmd3","Cd68","Ccl3","Lpl",
                         "Ctnna3","Ccl12","Ccl4","H2-Ab1","H2-Aa","Trem2","Ifit2","Ifit3"),
           annoTerm.data = enrich1,go.size = 15,
           go.col = c(rep(c("#FA7F6F","#82B0D2"),12),"#FA7F6F")
           )
dev.off()

col_fun = colorRamp2(c(0, 32), c("purple", "yellow"))
ha = HeatmapAnnotation(foo = add_pseudo, col = list(foo = col_fun))

#####
library(slingshot)
library(grDevices)
dimred <- intergrate2.traj.seurat@reductions$umap@cell.embeddings
clustering <- intergrate2.traj.seurat$celltype_new
counts <- as.matrix(intergrate2.traj.seurat@assays$RNA@counts[intergrate2.traj.seurat@assays$SCT@var.features, ])

sce <- SingleCellExperiment(assays = List(counts = counts))

assays(sce)$norm <- FQnorm(assays(sce)$counts)
reducedDims(sce) <- SimpleList(PCA = Embeddings(intergrate2.traj.seurat,reduction = "pca"), UMAP = Embeddings(intergrate2.traj.seurat,reduction = "umap"))

cl = intergrate2.traj.seurat@meta.data[rownames(colData(sce)),]$celltype_new
colData(sce)$mm = "test"
colData(sce)$cluster = cl

colData(sce) = colData(sce)[,c("mm","cluster")]
sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce)$UMAP, col = plotcol,  pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

library(tradeSeq)
sim <- fitGAM(sce)


tmp = as.data.frame(colData(sce)[,c(1,2,4)])
tmp$cell = rownames(tmp)
umap = as.data.frame(reducedDims(sce)$UMAP)
umap$cell = rownames(umap)

curves = slingCurves(sce, as.df = TRUE)

tmp = merge(tmp,umap,by = "cell")
p1 = ggplot(tmp, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = cluster),size = 1) +
  scale_color_manual(breaks = names(microglia.color3),values = microglia.color3)+
  geom_path(data = curves[curves$Lineage == 1,] %>% arrange(Order),aes(group = Lineage), arrow = arrow(length=unit(0.2, "cm")),size = 1) +
  theme_classic()+
  theme_dr()+
  theme(panel.grid = element_blank(),legend.position = 'none')
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/microglia/trajectory/slingshot_cl.png",p1,width = 8,height = 6,dpi = 300)

head(tmp)

p1 = ggplot(tmp, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "grey66",size = 1)+
  geom_point(data = na.omit(tmp),aes(x = UMAP_1, y = UMAP_2,color = slingPseudotime_1))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('yellow',"orange",'firebrick1',"black"))+
  geom_path(data = curves[curves$Lineage == 1,] %>% arrange(Order),aes(group = Lineage), arrow = arrow(length=unit(0.2, "cm")),size = 1,color = "cyan") +
  theme_dr()+
  theme(panel.grid = element_blank())
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/figures/slingshot_pseudotime_l.pdf",p1,width = 6,height = 6,dpi = 300)
ggsave("/Users/ruoqing/Projects/5xFAD_Tdepletion/figures/slingshot_pseudotime_nl.png",p1,width = 6,height = 6,dpi = 300)


tmp = modulated_genes[modulated_genes$status == "OK",]
tmp = tmp[order(tmp$morans_I,decreasing = T),]
i = "H2-Aa"
j = 8
p1 =FeaturePlot(intergrate2.traj.seurat,features = i,
                cols = brewer.pal(j, "OrRd"),raster=FALSE)+
    theme_dr()+
    theme(axis.line = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text = element_blank(),
          panel.background=element_blank(),
          axis.ticks = element_blank())
ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/figures/feature/",i,".pdf"),p1,width = 5,height = 4.5,dpi = 300)


p1 =FeaturePlot(intergrate2.traj.seurat,features = i,
                cols = brewer.pal(j, "OrRd"),raster=FALSE)+
  theme_bw()+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_blank(),
        panel.background=element_blank(),title = element_blank(),panel.border = element_blank(),
        axis.ticks = element_blank(),legend.position = 'none')
ggsave(paste0("/Users/ruoqing/Projects/5xFAD_Tdepletion/figures/feature/",i,".png"),p1,width = 4,height = 4,dpi = 300)


