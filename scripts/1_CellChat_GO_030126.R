#### library
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(harmony)
library(ArchR)
library(hexbin)
library(biovizBase)
library(ggplot2)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(statmod)
library(data.table)
library(foreach)




        # load data
        ct_intmkr = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_grandmarker_genes_010424.txt")
        gene_expr = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/DE_gene/DEGs_list_allcount_013125.txt")
        
        ECDEGs = gene_expr[cell_type == "EC"][padj < 0.05][log2FoldChange > 0.58]
        e1e3 = ct_intmkr[cluster %in% c("EC_E1","EC_E3")][p_val_adj < 0.05]#[avg_log2FC > 0.5]
        
        library(clusterProfiler)
        library(org.Hs.eg.db)
        
        ece1_go <- enrichGO(gene    = e1e3[cluster == "EC_E1"]$gene, # 275
                            OrgDb    = org.Hs.eg.db,
                            keyType  = 'SYMBOL',
                            ont      = "CC",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable = TRUE)
        
        ece3_go <- enrichGO(gene    = e1e3[cluster == "EC_E3"]$gene, # 486
                            OrgDb    = org.Hs.eg.db,
                            keyType  = 'SYMBOL',
                            ont      = "CC",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable = TRUE)
        
        ec_go <- enrichGO(gene    = ECDEGs$gene_name, # 245
                          OrgDb    = org.Hs.eg.db,
                          keyType  = 'SYMBOL',
                          ont      = "CC",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable = TRUE)
        
        ece1_go=as.data.table(as.data.frame(ece1_go))
        ece1_go$geneID<-NULL
        ece3_go=as.data.table(as.data.frame(ece3_go))
        ece3_go$geneID<-NULL
        ec_go=as.data.table(as.data.frame(ec_go))
        ec_go$geneID<-NULL
        
        go = rbind(ece1_go[grepl("collagen|signaling|proliferation|actin",Description)][, group := "EC_E1"], ece3_go[grepl("collagen|signaling|proliferation|actin",Description)][, group := "EC_E3"], ec_go[grepl("collagen|signaling|proliferation|actin",Description)][, group := "EC"])
        
        go[, `-log10(P-value)` := -log10(pvalue)]
        q.temp = do.call("rbind", strsplit(go$GeneRatio,"\\/"))
        
        go[, `Odds Ratio` := round(as.numeric(q.temp[,1])/as.numeric(q.temp[,2]),digits = 2)]
        go[, Term := paste0(Description, " (",ID,")")]
        
        ggplot(go[`-log10(P-value)` > 3], aes(x = `-log10(P-value)`, y = reorder(Term, `-log10(P-value)`))) + 
          # Lollipop sticks
          geom_segment(aes(x = 0, xend = `-log10(P-value)`,y = reorder(Term, `-log10(P-value)`), yend = reorder(Term, `-log10(P-value)`)),
                       color = "gray80",linewidth = 0.8) +
          # Left-aligned text labels
          geom_text(aes(x = 0, y = reorder(Term, `-log10(P-value)`),label = Term),
                    color = "black",size = 4,hjust = 0,vjust = -1.5,nudge_x = -0.1) +
          # Points
          geom_point(aes(size = `Odds Ratio`, color = group),alpha = 0.8,
          ) +
          # Scales
          scale_size_continuous(name = "Gene ratio",range = c(3, 6),limits = c(0, 0.2),breaks = c(0, 0.1, 0.2)) +
          scale_color_manual(name = "Cell type",values = c("EC" = "#FC8D62","EC_E1" = "#00BFFF", "EC_E3" = "#DDA0DD")) +
          labs(x = "-log10 (P-value)",y = NULL) +
          theme_classic() +
          theme(panel.grid.major.y = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(size = 12),
                axis.title.x = element_text(size = 12),axis.text.y = element_blank(),axis.ticks.y = element_blank(),
                legend.position = c(0.9, 0.2),legend.title = element_text(size = 12),legend.text = element_text(size = 10),legend.key.size = unit(1, "lines"))
        
        
        
        
        # load data
        dgene = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/DE_gene/DEGs_list_allcount_013125.txt")
        dgene[padj<0.05][cell_type == "EC"][log2FoldChange>0.58] # log2FC > 0.5
        
        EC_E_mrkr_gene = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_grandmarker_genes_010424.txt")
        EC_E_mrkr_gene[grepl("EC_E3",cluster)][p_val_adj < 0.05]
        EC_E_mrkr_gene[grepl("EC_E1",cluster)][p_val_adj < 0.05]
        
        #write.table(EC_E_mrkr_gene[grepl("EC_E3",cluster)][p_val_adj < 0.05][,c("gene"),with=F], file = "/Users/yangyu/Desktop/EC_E3.txt", quote = F, row.names = F, col.names = F, sep = "\t")
        #write.table(EC_E_mrkr_gene[grepl("EC_E1",cluster)][p_val_adj < 0.05][,c("gene"),with=F], file = "/Users/yangyu/Desktop/EC_E1.txt", quote = F, row.names = F, col.names = F, sep = "\t")
        #write.table(dgene[padj<0.05][cell_type == "EC"][log2FoldChange>0.58][,c("gene_name"),with=F], file = "/Users/yangyu/Desktop/EC_DE.txt", quote = F, row.names = F, col.names = F, sep = "\t")
        
        
        
        # LISA
        EC_E1_lisa = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/slides/Paper/Circulation/revision/data/RevisionMinors_022326/data_results/data_IKMPresults/EC_E1_lisa_022926.csv")
        EC_E3_lisa = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/slides/Paper/Circulation/revision/data/RevisionMinors_022326/data_results/data_IKMPresults/EC_E3_lisa_022926.csv")
        EC_lisa = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/slides/Paper/Circulation/revision/data/RevisionMinors_022326/data_results/data_IKMPresults/EC_DE_lisa_022926.csv")
        
        setnames(EC_E1_lisa, c("TF","rank","P","p2","p3","p4","p5"))
        setnames(EC_E3_lisa, c("TF","rank","P","p2","p3","p4","p5"))
        setnames(EC_lisa, c("TF","rank","P","p2","p3","p4","p5"))
        
        lisa = rbind(EC_E1_lisa[1:30,][,group:="EC_E1_marker"],EC_E3_lisa[1:30,][,group:="EC_E3_marker"],EC_lisa[1:30,][,group:="EC_DE"])
        
        common_TFs <- lisa[group %in% c("EC_E1_marker", "EC_E3_marker", "EC_DE"), .N, by = TF][N == 3, TF]
        lisa[TF %in% common_TFs]
        
        # Split TFs by Group
        tf_lists <- split(lisa$TF, lisa$group)
        
        library(VennDiagram)
        library(ggvenn)
        
        ggvenn(tf_lists,fill_color = c("#becdbe", "#6f91c8", "#6f91c8"),stroke_size = 0.5,set_name_size = 4)
        
        
        
        
        # CellChat results
        
        # load data # https://github.com/JiahaoWongg/CellChat_fix/blob/main/visualization.R
        combined = readRDS("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_combined_grandannt_010424.rds")
        
        seurat.EFEs = subset(combined, subset = SAMPLE %in% c("EFEs"))
        seurat.CTRLs = subset(combined, subset = SAMPLE %in% c("CTRLs"))
        
        
        # functions
        create_chatobjt = function(seurat_object = seurat.EFEs, pop_size = TRUE){
          
          data.input <- seurat_object[["RNA"]]@data
          labels <- Idents(seurat_object)
          meta <- data.frame(labels = labels, row.names = names(labels))
          
          cellChat <- createCellChat(object = seurat_object, group.by = "ident", assay = "RNA")
          
          # set the ligand receptor interaction database
          CellChatDB <- CellChatDB.human
          
          # use a subset of CellChatDB for cell-cell communication analysis
          CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
          
          # set the used database in the object
          cellChat@DB <- CellChatDB.use
          
          # Preprocessing the expression data for cell-cell communication analysis
          cellChat <- subsetData(cellChat)
          cellChat <- identifyOverExpressedGenes(cellChat)
          cellChat <- identifyOverExpressedInteractions(cellChat)
          
          # Compute the communication probability and infer cellular communication network
          cellChat <- computeCommunProb(cellChat, type = "truncatedMean", trim = 0.25, population.size = pop_size)
          cellChat <- filterCommunication(cellChat, min.cells = 10)
          
          # Infer the cell-cell communication at a signaling pathway level
          cellChat <- computeCommunProbPathway(cellChat)
          
          
          #Calculate the aggregated cell-cell communication network
          cellChat <- aggregateNet(cellChat)
          
          # Compute the network centrality scores
          cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
          
          # Idenfity signaling groups based on functional similarity
          cellChat<-computeNetSimilarity(cellChat,type="functional")
          cellChat<-netEmbedding(cellChat,type="functional")
          cellChat<-netClustering(cellChat,type="functional",do.parallel = FALSE)
          
          
          return(cellChat)
        }
        
        # create objt
        cellchat.EFEs = create_chatobjt(seurat_object = seurat.EFEs, pop_size = TRUE)
        cellchat.CTRLs = create_chatobjt(seurat_object = seurat.CTRLs, pop_size = TRUE)
        
        # merge objt
        object.list <- list(EFEs = cellchat.EFEs, CTRLs = cellchat.CTRLs)
        cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
        
        cellchat
        
        netVisual_heatmap(cellchat)
        
        # Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets
        library(RColorBrewer)
        
        cell_type_cols <- c(
          "EC" = "#FFA500",       
          "FB" = "#6A5ACD", 
          "Macrophage" = "#8FBC8F", 
          "Neuro" = "#A0522D",    
          "CM" = "#DC143C", 
          "EC_E1" = "#00BFFF",     
          "EC_E3" = "#DDA0DD", 
          "PeriC" = "#DB7093",  
          "SMC" = "#9932CC",       
          "EndoC" = "#FF69B4",    
          "T-cell" = "#87CEFA"   
        )
        
        # Get cell type order from the first object
        cell_types <- rownames(object.list[[1]]@net$count)
        
        # Get corresponding colors in correct order
        plot_colors <- cell_type_cols[cell_types]
        
        # Get max weight for scaling
        weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
        
        # Set up plot layout
        par(mfrow = c(1, length(object.list)), xpd = TRUE, mar = c(2, 2, 4, 2))
        
        # Generate plots
        for (i in 1:length(object.list)) {
          netVisual_circle(object.list[[i]]@net$count,weight.scale = TRUE,label.edge = FALSE,edge.weight.max = weight.max[2],
                           edge.width.max = 4,color.use = plot_colors, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                           vertex.size = 20,vertex.label.cex = 0.8)
        }
        
        # Chord diagram
        pathways.show <- c("SEMA3") 
        par(mfrow = c(1,2), xpd=TRUE)
        for (i in 1:length(object.list)) {
          netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", color.use = plot_colors, 
                              signaling.name = paste(pathways.show, names(object.list)[i]))
        }
        
        bubble_plot <- netVisual_bubble(cellchat, sources.use = 6, targets.use = c(3,1,9), sort.by.target = T, comparison = c(1, 2), angle.x = 45, return.data = TRUE) #KEY
        
        # The plot data is stored in the plot object
        plot_data <- bubble_plot$communication
        
        # Get min and max values
        min(plot_data$prob, na.rm = TRUE)
        max(plot_data$prob, na.rm = TRUE)
        
        
        # Circle differential functions
        netVisual_diffInteraction <- function(object, comparison = c(1,2), 
                                              measure = c("count", "weight", "count.merged", "weight.merged"), 
                                              color.use = NULL, 
                                              color.edge = c('#b2182b','#2166ac'), 
                                              title.name = NULL, 
                                              sources.use = NULL, 
                                              targets.use = NULL, 
                                              remove.isolate = FALSE, 
                                              top = 1,
                                              weight.scale = FALSE, 
                                              vertex.weight = 20, 
                                              vertex.weight.max = NULL, 
                                              vertex.size.max = 15, 
                                              vertex.label.cex=1,
                                              vertex.label.color= "black",
                                              edge.weight.max = NULL, 
                                              edge.width.max=8, 
                                              alpha.edge = 0.6, 
                                              label.edge = FALSE,
                                              edge.label.color='black',
                                              edge.label.cex=0.8,
                                              edge.curved=0.2,
                                              shape='circle',
                                              layout=in_circle(), 
                                              margin=0.2,
                                              arrow.width=1,
                                              arrow.size = 0.2
        ){
          options(warn = -1)
          measure <- match.arg(measure)
          
          obj1_raw <- object@net[[comparison[1]]][[measure]]
          obj2_raw <- object@net[[comparison[2]]][[measure]]
          
          shared_label = intersect(rownames(obj1_raw), rownames(obj2_raw))
          obj1 = obj1_raw[shared_label, shared_label]
          obj2 = obj2_raw[shared_label, shared_label]
          
          net.diff <- obj2 - obj1
          if (measure %in% c("count", "count.merged")) {
            if (is.null(title.name)) {
              title.name = "Differential number of interactions"
            }
          } else if (measure %in% c("weight", "weight.merged")) {
            if (is.null(title.name)) {
              title.name = "Differential interaction strength"
            }
          }
          net <- net.diff
          if ((!is.null(sources.use)) | (!is.null(targets.use))) {
            df.net <- reshape2::melt(net, value.name = "value")
            colnames(df.net)[1:2] <- c("source","target")
            # keep the interactions associated with sources and targets of interest
            if (!is.null(sources.use)){
              if (is.numeric(sources.use)) {
                sources.use <- rownames(net.diff)[sources.use]
              }
              df.net <- subset(df.net, source %in% sources.use)
            }
            if (!is.null(targets.use)){
              if (is.numeric(targets.use)) {
                targets.use <- rownames(net.diff)[targets.use]
              }
              df.net <- subset(df.net, target %in% targets.use)
            }
            cells.level <- rownames(net.diff)
            df.net$source <- factor(df.net$source, levels = cells.level)
            df.net$target <- factor(df.net$target, levels = cells.level)
            df.net$value[is.na(df.net$value)] <- 0
            net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
            net[is.na(net)] <- 0
          }
          
          if (remove.isolate) {
            idx1 <- which(Matrix::rowSums(net) == 0)
            idx2 <- which(Matrix::colSums(net) == 0)
            idx <- intersect(idx1, idx2)
            net <- net[-idx, ]
            net <- net[, -idx]
          }
          
          net[abs(net) < stats::quantile(abs(net), probs = 1-top, na.rm= T)] <- 0
          
          g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
          edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
          coords<-layout_(g,layout)
          if(nrow(coords)!=1){
            coords_scale=scale(coords)
          }else{
            coords_scale<-coords
          }
          if (is.null(color.use)) {
            color.use = scPalette(length(igraph::V(g)))
          }
          if (is.null(vertex.weight.max)) {
            vertex.weight.max <- max(vertex.weight)
          }
          vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
          
          loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
          
          # ERROR: Length of new attribute value must be 1 or 45
          # igraph::V(g)$size<-vertex.weight
          
          # Fix below
          igraph::V(g)$size<-vertex.weight
          
          igraph::V(g)$color<-color.use[igraph::V(g)]
          igraph::V(g)$frame.color <- color.use[igraph::V(g)]
          igraph::V(g)$label.color <- vertex.label.color
          igraph::V(g)$label.cex<-vertex.label.cex
          if(label.edge){
            igraph::E(g)$label<-igraph::E(g)$weight
            igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
          }
          igraph::E(g)$arrow.width<-arrow.width
          igraph::E(g)$arrow.size<-arrow.size
          igraph::E(g)$label.color<-edge.label.color
          igraph::E(g)$label.cex<-edge.label.cex
          #igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
          igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1],color.edge[2])
          igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, alpha.edge)
          
          igraph::E(g)$weight <- abs(igraph::E(g)$weight)
          
          if (is.null(edge.weight.max)) {
            edge.weight.max <- max(igraph::E(g)$weight)
          }
          if (weight.scale == TRUE) {
            #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
            igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
          }else{
            igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
          }
          
          
          if(sum(edge.start[,2]==edge.start[,1])!=0){
            # add this line
            igraph::E(g)$loop.angle = rep(0, nrow(edge.start))
            igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
          }
          
          radian.rescale <- function(x, start=0, direction=1) {
            c.rotate <- function(x) (x + start) %% (2 * pi) * direction
            c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
          }
          label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
          label.dist <- vertex.weight/max(vertex.weight)+2
          plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
               vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
          if (!is.null(title.name)) {
            text(0,1.5,title.name, cex = 1.1)
          }
          # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
          # grid.echo()
          # gg <-  grid.grab()
          gg <- recordPlot()
          return(gg)
        }
        
        
        
        
        
        # Circle plots
        netVisual_diffInteraction(cellchat, 
                                  measure = c("count"),
                                  weight.scale = TRUE,
                                  edge.width.max = 4,
                                  color.use = plot_colors,
                                  vertex.size = 20,
                                  vertex.label.cex = 0.8,
                                  title.name = "Differential Network",
                                  comparison = c(1, 2))  # Add this
        
        
        
        
        
        # save
        saveRDS(cellchat, file = "/Users/yangyu/Desktop/Harvard_ChenLab/EP/slides/Paper/Circulation/revision/data/RevisionMinors_022326/data_results/cellchat_comparisonAnalysis_EFECTRL_030126.rds")
        
        
        
        
        
        
        
        
        
        
        
        
        
