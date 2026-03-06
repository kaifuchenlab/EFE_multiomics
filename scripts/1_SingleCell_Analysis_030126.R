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



            ###################
            # This is a script of the core single cell analysis
            # Summarized by YY 030126
            # In the 'Cell_Type' column, which labels annotated cell types,'EC' == 'EC_1', 'EC_E1' == 'EC_2' and 'EC_E3' == 'EC_3'
            ###################
            
            
            
            
            # Fig.B 
            
            # load data
            combined = readRDS("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_combined_grandannt_010424.rds")
            
            # colors
            #library(scRNAtoolVis)
            library(RColorBrewer)
            cell_type_cols=c("#FFA500","#00BFFF","#DDA0DD","#FF69B4","#6A5ACD","#9932CC",
                             "#DC143C","#8FBC8F","#87CEFA","#DB7093","#A0522D","#C71585",
                             "#32CD32","#00CED1","#8B0000","#BC8F8F","#F08080","#DAA520",
                             "#90EE90","#3CB371","#191970")
            
            # assign groups
            combined <- SetIdent(combined, value = "SAMPLE")
            
            # color
            combined$SAMPLE=factor(combined$SAMPLE,levels=c("EFEs","CTRLs","Fetals","Youngs"))
            condition_colors <- c("CTRLs" = "darkgreen","EFEs" = "#E41A1C","Fetals" = "#87CEFA","Youngs" = "#6A5ACD")    # Red for EFEs
            
            # plot
            DimPlot(combined, cols = condition_colors, pt.size = 0.001, raster = FALSE, shuffle = TRUE) #+ NoLegend()
            
            
            
            # Fig.C
            combined <- SetIdent(combined, value = "Cell_Type")
            cell_order = c("EC","EC_E1","EC_E3","EndoC","FB","SMC","CM","Macrophage", "T-cell","PeriC","Neuro","Mast-cell")
            combined@active.ident=factor(combined@active.ident,levels=cell_order)
            
            combined[["SP"]] <- ""
            combined@meta.data$SP[combined@meta.data$SAMPLE %in% c("EFEs")] <- "EFEs"
            combined@meta.data$SP[combined@meta.data$SAMPLE %in% c("CTRLs","Fetals","Youngs")] <- "CTRLs"
            
            combined$SP <- factor(combined$SP, levels = c("EFEs","CTRLs"))
            DimPlot(combined, cols = cell_type_cols, split.by = "SP", pt.size = 0.001) + NoLegend()
            
            
            
            # Fig.D
            T_cell=c("CD247")
            
            EC=c("PECAM1","VWF","BTNL9")
            EndoC=c("SMOC1","NPR3")
            
            Fibroblast=c("DCN","COL1A1") 
            
            cardiomyocyte=c("TNNT2", "TTN") 
            
            M2_macrophage=c("CD163","MRC1") 
            
            Pericyte=c("RGS5","PDGFRB") 
            smoothe_muscle=c("ACTA2","MYH11") 
            
            Neuro=c("XKR4") 
            
            Mast_cell=c("TPSAB1")
            
            # assign another
            young_se1=combined
            
            young_se1 <- SetIdent(young_se1, value = "Cell_Type")
            young_se1$Cell_Type=factor(young_se1$Cell_Type,levels=cell_order)
            
            # cluster Marker gene plot
            marker_gene=c(EC, EndoC,Fibroblast,smoothe_muscle,cardiomyocyte,M2_macrophage,T_cell,Pericyte,Neuro,Mast_cell)
            dat=list(EC,EndoC,Fibroblast,smoothe_muscle,cardiomyocyte,M2_macrophage,T_cell,Pericyte,Neuro,Mast_cell)
            cluster_gene=NULL
            for (i in 1:length(dat)){a=rep(i,length(dat[[i]]))
            cluster_gene=c(cluster_gene,a)}
            
            dat1=data.frame(marker_gene,cluster_gene)
            dat1=dat1[!duplicated(dat1$marker_gene),]
            dat1=dat1[order(dat1$cluster_gene,decreasing = T),]
            
            # Plot
            library(scRNAtoolVis)
            
            # plot
            p1=jjDotPlot(object = young_se1,gene = dat1$marker_gene,id = 'Cell_Type',rescale = T,ytree = F, y.sort = cell_order)
            p1
            
            
            
            
            # Fig.E
            
            # proportions
            library(speckle)
            
            # load data               
            combined = readRDS("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_combined_grandannt_010424.rds")
            
            condition_colors <- c("CTRLs" = "darkgreen","EFEs" = "#E41A1C","Fetals" = "#87CEFA","Youngs" = "#6A5ACD")    # Red for EFEs
            pred_label = as.data.table(combined@meta.data, keep.rownames = T)
            
            # Calculate proportions
            prop = as.data.table(propeller(clusters = pred_label$coarse_type, sample = pred_label$STIM, group = pred_label$SAMPLE),keep.rownames = T)
            
            # separate comparisions
            as.data.table(propeller(clusters = pred_label[SAMPLE %in% c("CTRLs","EFEs")]$coarse_type, sample = pred_label[SAMPLE %in% c("CTRLs","EFEs")]$STIM, group = pred_label[SAMPLE %in% c("CTRLs","EFEs")]$SAMPLE),keep.rownames = T)[FDR<0.05]
            as.data.table(propeller(clusters = pred_label[SAMPLE %in% c("Fetals","EFEs")]$coarse_type, sample = pred_label[SAMPLE %in% c("Fetals","EFEs")]$STIM, group = pred_label[SAMPLE %in% c("Fetals","EFEs")]$SAMPLE),keep.rownames = T)[FDR<0.05]
            as.data.table(propeller(clusters = pred_label[SAMPLE %in% c("Youngs","EFEs")]$coarse_type, sample = pred_label[SAMPLE %in% c("Youngs","EFEs")]$STIM, group = pred_label[SAMPLE %in% c("Youngs","EFEs")]$SAMPLE),keep.rownames = T)[FDR<0.05]
            
            # anova
            as.data.table(propeller(clusters = pred_label$coarse_type, sample = pred_label$STIM, group = pred_label$SAMPLE),keep.rownames = T)
            
            # for plot
            prop_plot <- pred_label[, .(count = .N), by = .(SAMPLE, STIM, coarse_type)][, proportion := count/sum(count), by = .(SAMPLE, STIM)][,mean_prop := round(mean(proportion),digits = 4), list(SAMPLE,coarse_type)]
            
            
            # Export datatables
            Pcoa_ctrl = as.data.table(propeller(clusters = pred_label[SAMPLE %in% c("CTRLs","EFEs")]$coarse_type, sample = pred_label[SAMPLE %in% c("CTRLs","EFEs")]$STIM, group = pred_label[SAMPLE %in% c("CTRLs","EFEs")]$SAMPLE),keep.rownames = T)
            Pcoa_fetal = as.data.table(propeller(clusters = pred_label[SAMPLE %in% c("Fetals","EFEs")]$coarse_type, sample = pred_label[SAMPLE %in% c("Fetals","EFEs")]$STIM, group = pred_label[SAMPLE %in% c("Fetals","EFEs")]$SAMPLE),keep.rownames = T)
            Pcoa_young = as.data.table(propeller(clusters = pred_label[SAMPLE %in% c("Youngs","EFEs")]$coarse_type, sample = pred_label[SAMPLE %in% c("Youngs","EFEs")]$STIM, group = pred_label[SAMPLE %in% c("Youngs","EFEs")]$SAMPLE),keep.rownames = T)
            
            
            
            # StackedBarplot 022526
            prop_stats <- pred_label[, .(count = .N), by = .(SAMPLE, STIM, coarse_type)][
              , proportion := count/sum(count), by = .(SAMPLE, STIM)][
                , .(mean_prop = mean(proportion),
                    se = sd(proportion)/sqrt(.N),  # Standard error
                    individual_prop = proportion),  # Keep original points
                by = .(SAMPLE, coarse_type)]
            
            # reorder the cell label
            prop_stats$SAMPLE <- factor(prop_stats$SAMPLE, levels = rev(c("EFEs","CTRLs","Fetals","Youngs")))
            prop_stats$coarse_type <- factor(prop_stats$coarse_type, levels = c("EC","EndoC","Fibroblast","SMC","Cardiomyocyte","Macrophage", "T-cell","Pericyte","Neuro","Mast-cell"))
            
            # plot
            ggplot(prop_stats[!duplicated(mean_prop)], aes(x = mean_prop, y = SAMPLE, fill = coarse_type)) +
              geom_col(position = "stack", width = 0.7) +
              labs(x = "Proportion", y = "Sample", fill = "Coarse Type") + 
              theme_classic() +
              theme(legend.position = "right",
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size = 14)) +
              scale_fill_manual(values = cell_type_cols[c(1,4:12)])
            
            # StackedBarplot 022526
            prop_stats <- pred_label[, .(count = .N), by = .(SAMPLE, STIM, Cell_Type)][
              , proportion := count/sum(count), by = .(SAMPLE, STIM)][
                , .(mean_prop = mean(proportion),
                    se = sd(proportion)/sqrt(.N),  # Standard error
                    individual_prop = proportion),  # Keep original points
                by = .(SAMPLE, Cell_Type)]
            
            # reorder the cell label
            prop_stats$SAMPLE <- factor(prop_stats$SAMPLE, levels = rev(c("EFEs","CTRLs","Fetals","Youngs")))
            prop_stats$Cell_Type <- factor(prop_stats$Cell_Type, levels = c("EC","EC_E1","EC_E3","EndoC","FB","SMC","CM","Macrophage", "T-cell","PeriC","Neuro","Mast-cell"))
            
            
            # plot
            ggplot(prop_stats[!duplicated(mean_prop)], aes(x = mean_prop, y = SAMPLE, fill = Cell_Type)) +
              geom_col(position = "stack", width = 0.7) +
              labs(x = "Proportion", y = "Sample", fill = "Coarse Type") + 
              theme_classic() +
              theme(legend.position = "right",
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size = 14)) +
              scale_fill_manual(values = cell_type_cols)
            
            
            # Fig. Sankey
            # Sankey digrams 3 x 3.8
            library(ggalluvial)
            
            comb_s = subset(combined, subset = SAMPLE %in% c("EFEs","CTRLs"))
            gc()
            
            comb_s$STIM=factor(comb_s$STIM,levels=c("EFE1","EFE2","EFE3","CTRL1","CTRL2","CTRL3"))
            only_EC = subset(comb_s, subset = Cell_Type %in% c("EC","EC_E1","EC_E3"))
            
            only_EC <- SetIdent(only_EC, value = "Cell_Type")
            pt <- table(Idents(only_EC), only_EC$STIM)
            pt <- as.data.table(as.data.frame(pt))
            pt$Var1 <- as.character(pt$Var1)
            
            # total per STIM cells
            sum_cell = as.data.table(as.data.frame(table(comb_s$STIM)))
            setnames(sum_cell, c("Var2","sum_cell"))
            
            pt = merge(setkey(pt, Var2), setkey(sum_cell, Var2))
            pt[, prop := round(Freq/sum_cell, digits = 2)]
            
            pt$Var2 <- factor(pt$Var2,levels = c("EFE1","EFE3","EFE2","CTRL2","CTRL1","CTRL3")) # change this order for showing differences from samp_order to other orders
            
            ggplot(pt,
                   aes(x = Var2,               # x-axis variable (STIM groups)
                       stratum = Var1,         # cell types
                       y = prop,               # use proportion instead of count
                       fill = Var1,            # fill by cell type
                       alluvium = Var1)) +     # grouping for flows
              geom_flow(width = 1/4, alpha = 0.6) +
              geom_stratum(width = 1/4) +
              #scale_y_continuous(labels = scales::percent) + # Show as percentages
              theme_classic(base_size = 14) +  # Base font size (affects all text)
              labs(x = "Sample",
                   y = "Cell proportion") +
              scale_fill_manual(values = cell_type_cols) + 
              theme(
                # Legend customization
                legend.title = element_blank(),
                legend.text = element_text(size = 14),  # Legend item text size
                legend.position = c(0.8, 0.8),         # Legend position (x,y)
                legend.background = element_blank(),    # Remove legend background
                
                # Axis titles
                axis.title.x = element_blank(),        # Remove x-axis title
                axis.title.y = element_text(size = 16),  # y-axis title size
                
                # Axis text (tick labels)
                axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),  # x-axis labels
                axis.text.y = element_text(size = 14))
            
            
            
            # Fig.G
            # load data
            
            # plot function
            chrom_gene_compare_plot = function(data = data_pair, pairedTest = TRUE, yPos = 20, yHeightPara = 0.8, colorRange = c("#FFA500", "#00BFFF","#DDA0DD"),yTitle = "Chromatin accessibility for\nEC marker genes",dogeWidth = 0.5){
              # the stats
              my_comparisons <- combn(levels(data$variable), 2, simplify = FALSE)
              
              if(is.null(yPos)){
                y_positions <- seq(max(data$value)*yHeightPara, 
                                   by = max(data$value)*0.1*yHeightPara, 
                                   length.out = length(my_comparisons))
              }else{
                y_positions <- seq(yPos*yHeightPara, 
                                   by = yPos*0.1*yHeightPara, 
                                   length.out = length(my_comparisons))
              }
              
              # the actual plot
              p = ggplot(data, aes(x = variable, y = value, fill = variable)) + 
                geom_boxplot(width = 0.4,alpha = 0.7,              # Remove transparency for print
                             lwd = 0.5, fatten = 1.2, outlier.shape = NA,coef = 1.5, position = position_dodge(width = dogeWidth)) + 
                stat_boxplot(geom = "errorbar",width = 0.1,size = 0.5,color = "black",
                             aes(ymin = after_stat(ymax), ymax = after_stat(ymax)), position = position_dodge(width = dogeWidth)) +
                stat_boxplot(geom = "errorbar",width = 0.1,size = 0.5,color = "black",  # Or match your box colors
                             aes(ymin = after_stat(ymin),ymax = after_stat(ymin)), position = position_dodge(width = dogeWidth)) + 
                labs(title = "",x = "Cell Type",y = yTitle) +
                theme_classic() + 
                theme(legend.position = "none", panel.spacing.x = unit(0, "cm"),
                      axis.text.x = element_text(size = 12, angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
                      axis.text = element_text(size = 12),axis.title = element_text(size = 14)) + 
                stat_compare_means(
                  comparisons = my_comparisons,
                  label = "p.signif",  # label = "p.format", actual p formats
                  #hide.ns = TRUE,
                  label.y = y_positions,
                  method = "t.test",
                  paired = pairedTest, vjust = 0.5, size = 8,
                  tip.length = 0.05, # VERY short error bars (default: 0.03)
                  bracket.size = 0.2, # Thinner comparison lines
                  step.increase = 0.1
                ) + scale_fill_manual(values = colorRange)
              
              return(p)
            }
            
            # Fig.G chromatin accessibility and gene expression
            
            # load markers
            marker_genelist = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_coarsemarker_genes_011024.txt")
            fibro_gene = c("POSTN", "FN1", "TNC", "COL1A1", "COL1A2", "COL3A1", "DEC1", "RUNX1", "ADAM12", "KIF26B", "THBS4", "FAP", "COL5A1", "SERPINE1", "SLC20A1", "KALRN", "PRICKLE1", "CDH11", "FGF14", "THBS2")
            
            combined = readRDS("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_combined_grandannt_010424.rds")
            combined = subset(combined, subset = SAMPLE %in% c("EFEs","CTRLs"))
            gc()
            combined.atac = readRDS("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/combined_seurat_object_110624.rds")
            
            # select marker gene
            log_fc = 1 
            cell_type = "EC"
            cell_marker_gene_list = marker_genelist[cluster == cell_type][avg_log2FC > log_fc]$gene 
            
            # processings
            DefaultAssay(combined.atac) <- "ACTIVITY"
            
            # chromatin accessibility
            ctrl_activity <- as.data.table(as.matrix(AverageExpression(subset(combined.atac, subset = SAMPLE == "CTRLs"),assays = "ACTIVITY",group.by = "Cell_Type",slot = "data")$ACTIVITY),keep.rownames = T)     
            efes_activity <- as.data.table(as.matrix(AverageExpression(subset(combined.atac, subset = SAMPLE == "EFEs"),assays = "ACTIVITY",group.by = "Cell_Type",slot = "data")$ACTIVITY),keep.rownames = T)
            
            # cmpare
            pair_EC = merge(setkey(ctrl_activity[rn %in% cell_marker_gene_list], rn), setkey(efes_activity[rn %in% cell_marker_gene_list], rn))
            melt_pair = as.data.table(melt(pair_EC, id.vars = "rn"))
            melt_pair[,SAMPLE:="EFEs"][grepl(".x",variable),SAMPLE:="CTRLs"]
            melt_pair$variable <- gsub(".x|.y","",melt_pair$variable)
            
            # cell colors
            library(ggpubr)
            
            cell_type_cols=c("#FFA500","#00BFFF","#DDA0DD","#FF69B4","#6A5ACD","#9932CC",
                             "#DC143C","#8FBC8F","#87CEFA","#DB7093","#A0522D","#C71585",
                             "#32CD32","#00CED1","#8B0000","#BC8F8F","#F08080","#DAA520",
                             "#90EE90","#3CB371","#191970")
            
            
            # plots take all EC and EFE subtypes
            data_pair = melt_pair[(variable == "EC" & SAMPLE == "CTRLs") | (variable == "EC" & SAMPLE == "EFEs") | (variable == "EC_E1" & SAMPLE == "EFEs") | (variable == "EC_E3" & SAMPLE == "EFEs")]
            data_pair[variable == "EC" & SAMPLE == "CTRLs", variable := "EC_CTRLs"][variable == "EC" & SAMPLE == "EFEs", variable := "EC_EFEs"][variable == "EC_E1" & SAMPLE == "EFEs", variable := "EC_E1"][variable == "EC_E3" & SAMPLE == "EFEs", variable := "EC_E3"] 
            data_pair$variable <- factor(data_pair$variable, levels = c("EC_CTRLs","EC_EFEs","EC_E1","EC_E3"))
            
            # function plot 3 x 3
            p_chr_acc_EC = chrom_gene_compare_plot(data = data_pair, pairedTest = TRUE, yPos = 7.2, yHeightPara = 0.8, colorRange = c("#FFA500","#FFA500", "#00BFFF","#DDA0DD"),yTitle = "Chromatin accessibility for\nEC marker genes",dogeWidth = 0.5)
            p_chr_acc_EC + ylim(c(0,9.5))
            
            p_chr_acc_EC + ylim(c(0,6))
            data_pair1 = data_pair
            
            # Gene expression for EC and fibrosis genes here
            DefaultAssay(combined) <- "RNA"
            
            # expression
            ctrl_expression <- as.data.table(as.matrix(AverageExpression(subset(combined, subset = SAMPLE == "CTRLs"),assays = "RNA",group.by = "Cell_Type",slot = "data")$RNA),keep.rownames = T)     
            efes_expression <- as.data.table(as.matrix(AverageExpression(subset(combined, subset = SAMPLE == "EFEs"),assays = "RNA",group.by = "Cell_Type",slot = "data")$RNA),keep.rownames = T)     
            
            # paired boxplots
            pair_EC = merge(setkey(ctrl_expression[rn %in% cell_marker_gene_list], rn), setkey(efes_expression[rn %in% cell_marker_gene_list], rn))
            
            pair_EC = pair_EC[rn %in% ctrl_activity[rn %in% cell_marker_gene_list]$rn]
            melt_pair = as.data.table(melt(pair_EC, id.vars = "rn"))
            melt_pair[,SAMPLE:="EFEs"][grepl(".x",variable),SAMPLE:="CTRLs"]
            melt_pair$variable <- gsub(".x|.y","",melt_pair$variable)
            
            # plots
            data_pair = melt_pair[(variable == "EC" & SAMPLE == "CTRLs") | (variable == "EC" & SAMPLE == "EFEs") | (variable == "EC_E1" & SAMPLE == "EFEs") | (variable == "EC_E3" & SAMPLE == "EFEs")]
            data_pair[variable == "EC" & SAMPLE == "CTRLs", variable := "EC_CTRLs"][variable == "EC" & SAMPLE == "EFEs", variable := "EC_EFEs"][variable == "EC_E1", variable := "EC_E1"][variable == "EC_E3", variable := "EC_E3"]
            data_pair$variable <- factor(data_pair$variable, levels = c("EC_CTRLs","EC_EFEs","EC_E1","EC_E3"))
            
            
            # plot 3 x 3
            p_chr_RNA_ecall = chrom_gene_compare_plot(data = data_pair, pairedTest = TRUE, yPos = 22, yHeightPara = 0.8, colorRange = c("#FFA500","#FFA500", "#00BFFF","#DDA0DD"),yTitle = "Gene expression for\nEC marker genes",dogeWidth = 0.5)
            p_chr_RNA_ecall + ylim(c(0,37))
            
            p_chr_RNA_ecall + ylim(c(0,20))
            data_pair2 = data_pair
            
            # paired boxplots for fibrosis
            pair_EC = merge(setkey(ctrl_activity[rn %in% fibro_gene], rn), setkey(efes_activity[rn %in% fibro_gene], rn))
            melt_pair = as.data.table(melt(pair_EC, id.vars = "rn"))
            melt_pair[,SAMPLE:="EFEs"][grepl(".x",variable),SAMPLE:="CTRLs"]
            melt_pair$variable <- gsub(".x|.y","",melt_pair$variable)
            
            # plots
            data_pair = melt_pair[(variable == "EC" & SAMPLE == "CTRLs") | (variable == "EC" & SAMPLE == "EFEs") | (variable == "EC_E1" & SAMPLE == "EFEs") | (variable == "EC_E3" & SAMPLE == "EFEs")]
            data_pair[variable == "EC" & SAMPLE == "CTRLs", variable := "EC_CTRLs"][variable == "EC" & SAMPLE == "EFEs", variable := "EC_EFEs"][variable == "EC_E1", variable := "EC_E1"][variable == "EC_E3", variable := "EC_E3"]
            data_pair$variable <- factor(data_pair$variable, levels = c("EC_CTRLs","EC_EFEs","EC_E1","EC_E3"))
            
            # plot
            p_chr_acc_fibrosis = chrom_gene_compare_plot(data = data_pair, pairedTest = TRUE, yPos = 3, yHeightPara = 0.8, colorRange = c("#FFA500","#FFA500",  "#00BFFF","#DDA0DD"),yTitle = "Chromatin accessibility for\nfibrosis genes",dogeWidth = 0.5)
            p_chr_acc_fibrosis
            
            p_chr_acc_fibrosis + ylim(c(0,3))
            data_pair3 = data_pair
            
            # RNA for fibrosis genes
            pair_EC = merge(setkey(ctrl_expression[rn %in% fibro_gene], rn), setkey(efes_expression[rn %in% fibro_gene], rn))
            
            pair_EC = pair_EC[rn %in% ctrl_activity[rn %in% fibro_gene]$rn]
            
            melt_pair = as.data.table(melt(pair_EC, id.vars = "rn"))
            melt_pair[,SAMPLE:="EFEs"][grepl(".x",variable),SAMPLE:="CTRLs"]
            melt_pair$variable <- gsub(".x|.y","",melt_pair$variable)
            
            # plots
            data_pair = melt_pair[(variable == "EC" & SAMPLE == "CTRLs") | (variable == "EC" & SAMPLE == "EFEs") | (variable == "EC_E1" & SAMPLE == "EFEs") | (variable == "EC_E3" & SAMPLE == "EFEs")]
            data_pair[variable == "EC" & SAMPLE == "CTRLs", variable := "EC_CTRLs"][variable == "EC" & SAMPLE == "EFEs", variable := "EC_EFEs"][variable == "EC_E1", variable := "EC_E1"][variable == "EC_E3", variable := "EC_E3"]
            data_pair$variable <- factor(data_pair$variable, levels = c("EC_CTRLs","EC_EFEs","EC_E1","EC_E3"))
            
            # plot
            p_chr_RNA_fibrosis = chrom_gene_compare_plot(data = data_pair, pairedTest = TRUE, yPos = 7, yHeightPara = 0.8, colorRange = c("#FFA500","#FFA500","#00BFFF","#DDA0DD"),yTitle = "Gene expression for\nfibrosis genes",dogeWidth = 0.5)
            p_chr_RNA_fibrosis + ylim(c(0,9))
            
            p_chr_RNA_fibrosis + ylim(c(0,8))
            data_pair4 = data_pair
            
            
            # Actual plots
            p_chr_acc_EC + geom_point(data = data_pair1[rn %in% c("PECAM1")], aes(x=variable,y=value,group=rn), color = "red") + 
              geom_line(data = data_pair1[rn %in% c("PECAM1")], aes(x=variable,y=value,group=rn), 
                        color = "red", linetype = "dashed") + ylim(c(0,6))
            
            p_chr_RNA_ecall + geom_point(data = data_pair2[rn %in% c("PECAM1")], aes(x=variable,y=value,group=rn), color = "red") + 
              geom_line(data = data_pair2[rn %in% c("PECAM1")], aes(x=variable,y=value,group=rn), 
                        color = "red", linetype = "dashed") + ylim(c(0,20))
            
            p_chr_acc_fibrosis + geom_point(data = data_pair3[rn %in% c("COL1A1")], aes(x=variable,y=value,group=rn), color = "black") + 
              geom_line(data = data_pair3[rn %in% c("COL1A1")], aes(x=variable,y=value,group=rn), 
                        color = "black", linetype = "dashed") + ylim(c(0,3))
            
            p_chr_RNA_fibrosis + geom_point(data = data_pair4[rn %in% c("COL1A1")], aes(x=variable,y=value,group=rn), color = "black") + 
              geom_line(data = data_pair4[rn %in% c("COL1A1")], aes(x=variable,y=value,group=rn), 
                        color = "black", linetype = "dashed") + ylim(c(0,8))
            
            
            # Save datatables
            #write.table(data_pair1, file = "/Users/yangyu/Desktop/FigG_ECmarker_ATAC_022226.txt", quote = F, row.names = F, col.names = T, sep = "\t")
            #write.table(data_pair2, file = "/Users/yangyu/Desktop/FigG_ECmarker_RNA_022226.txt", quote = F, row.names = F, col.names = T, sep = "\t")
            #write.table(data_pair3, file = "/Users/yangyu/Desktop/FigG_Fibrosis_ATAC_022226.txt", quote = F, row.names = F, col.names = T, sep = "\t")
            #write.table(data_pair4, file = "/Users/yangyu/Desktop/FigG_Fibrosis_RNA_022226.txt", quote = F, row.names = F, col.names = T, sep = "\t")
            
            
            
            
            
            
            
            # Fig.H endoMT score
            library(UCell)
            library(ggpubr)
            
            TGFB = c("TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3")
            SMAD = c("SMAD2","SMAD3","SMAD5") # "SMAD5"
            BMPD = c("BMP2", "BMP5", "BMP7") # "BMP5", "BMP7"
            SNAIL = c("SNAI1","SNAI2","SNAI3")
            TWIST = c("TWIST1","TWIST2")
            
            # assign markers
            endo_gene = c(TGFB, SMAD, BMPD, SNAIL, TWIST)
            
            markers <- list()
            markers$endoMT <- endo_gene
            
            # Add module
            combined <- AddModuleScore_UCell(combined, slot = "data",features = markers)
            signature.names <- paste0(names(markers), "_UCell")
            
            # fibrosis score
            Efe_ctrl = subset(combined, subset = SAMPLE %in% c("EFEs","CTRLs"))
            Efe_ctrl$SAMPLE <- factor(Efe_ctrl$SAMPLE, levels = c("EFEs","CTRLs"))
            
            # First create the base plot
            sample_levels <- unique(Efe_ctrl@meta.data$SAMPLE)
            condition_colors <- c("EFEs" = "#E41A1C","CTRLs" = "darkgreen","Fetals" = "#87CEFA","Youngs" = "#6A5ACD")[sample_levels]    # Red for EFEs
            
            # Extract the data used in the plot
            plot_data <- as.data.frame(Efe_ctrl[[c("Cell_Type", "SAMPLE")]])
            expr_data <- FetchData(Efe_ctrl, vars = signature.names)
            plot_data <- cbind(plot_data, expr_data)
            
            # compares
            endo_MT = as.data.table(plot_data, keep.rownames = T)
            
            # plots
            data_pair = endo_MT[(Cell_Type == "EC") | (Cell_Type == "EC_E1" & SAMPLE == "EFEs") | (Cell_Type == "EC_E3" & SAMPLE == "EFEs")]
            data_pair[Cell_Type == "EC" & SAMPLE == "CTRLs", Cell_Type := "EC_CTRLs"][Cell_Type == "EC" & SAMPLE == "EFEs", Cell_Type := "EC_EFEs"][Cell_Type == "EC_E1", Cell_Type := "EC_E1"][Cell_Type == "EC_E3", Cell_Type := "EC_E3"]
            data_pair$Cell_Type <- factor(data_pair$Cell_Type, levels = c("EC_CTRLs","EC_EFEs","EC_E1","EC_E3"))
            
            my_comparisons <- combn(levels(data_pair$Cell_Type), 2, simplify = FALSE)
            y_positions <- seq(max(data_pair$endoMT_UCell)*1.15, 
                               by = max(data_pair$endoMT_UCell)*0.1, 
                               length.out = length(my_comparisons))
            
            # Plot
            dogeWidth = 0.5
            
            my_comparisons <- combn(levels(data_pair$Cell_Type), 2, simplify = FALSE)
            y_positions <- seq(0.2*1.15, 
                               by = 0.2*0.1, 
                               length.out = length(my_comparisons))
            
            # Boxplot
            ggplot(data_pair, aes(x = Cell_Type, y = endoMT_UCell, fill = Cell_Type)) + 
              #geom_jitter(width = 0.05, height = 0, size = 0.001, alpha = 0.4, aes(color = Cell_Type)) + 
              geom_jitter(width = 0.05, height = 0, size = 0.001, alpha = 0.4, color = "black") + 
              geom_boxplot(width = 0.4,alpha = 0.7,              # Remove transparency for print
                           lwd = 0.5, fatten = 1.2, outlier.shape = NA,coef = 1.5, position = position_dodge(width = dogeWidth)) + 
              #geom_jitter(width = 0.05, height = 0, size = 0.001, alpha = 0.4, aes(color = Cell_Type)) + 
              stat_boxplot(geom = "errorbar",width = 0.1,size = 0.5,color = "black",
                           aes(ymin = after_stat(ymax), ymax = after_stat(ymax)), position = position_dodge(width = dogeWidth)) +
              stat_boxplot(geom = "errorbar",width = 0.1,size = 0.5,color = "black",  # Or match your box colors
                           aes(ymin = after_stat(ymin),ymax = after_stat(ymin)), position = position_dodge(width = dogeWidth)) + 
              labs(title = "",x = "Cell Type",y = "EndoMT score") +
              theme_classic() + 
              theme(legend.position = "none", panel.spacing.x = unit(0, "cm"),
                    axis.text.x = element_text(size = 12, angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
                    axis.text = element_text(size = 12),axis.title = element_text(size = 14)) + 
              stat_compare_means(
                comparisons = my_comparisons,
                label = "p.signif",  # label = "p.format", actual p formats
                #hide.ns = TRUE,
                label.y = y_positions,
                method = "t.test", 
                paired = FALSE, vjust = 0.5, size = 8,
                tip.length = 0.05, 
                bracket.size = 0.2, # Thinner comparison lines
                step.increase = 0.1
              ) + scale_fill_manual(values = c("#FFA500","#FFA500", "#00BFFF","#DDA0DD")) + 
              scale_color_manual(values = c("#FFA500","#FFA500", "#00BFFF","#DDA0DD")) + 
              ylim(c(0,0.5)) # ylim(c(0,0.25))
            
            
            
            
            
            
            # Fig.I
            # Fig.I GO:CC
            ct_intmkr = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_grandmarker_genes_010424.txt")
            gene_expr = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/DE_gene/DEGs_list_allcount_013125.txt")
            
            ECDEGs = gene_expr[cell_type == "EC"][padj < 0.05][log2FoldChange > 0.58]
            e1e3 = ct_intmkr[cluster %in% c("EC_E1","EC_E3")][p_val_adj < 0.05]
            
            library(clusterProfiler)
            library(org.Hs.eg.db)
            
            ece1_go <- enrichGO(gene    = e1e3[cluster == "EC_E1"]$gene,
                                OrgDb    = org.Hs.eg.db,
                                keyType  = 'SYMBOL',
                                ont      = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05,
                                readable = TRUE)
            
            ece3_go <- enrichGO(gene    = e1e3[cluster == "EC_E3"]$gene,
                                OrgDb    = org.Hs.eg.db,
                                keyType  = 'SYMBOL',
                                ont      = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05,
                                readable = TRUE)
            
            ec_go <- enrichGO(gene    = ECDEGs$gene_name,
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
            
            
            
            
            
            # Fig.J,L
            combined.atac=readRDS("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/combined_seurat_object_110624.rds")
            DefaultAssay(combined.atac) <- "ACTIVITY"
            
            
            # Subset the Seurat object for selected cell types
            selected_cell_types = c("EC","EC_E1","CM")
            gene_list = c("VWF","FLI1","PECAM1","BTNL9","PIK3R3","ABLIM3","DKK2","SLC26A3","LINGO1","LMOD1","MEF2D","SORBS2","CORIN","MYOM1","RBM24","FHL2","MYH7B","TTN","TNNT2")
            
            # RNA and ATAC careful
            selected_cells <- subset(combined, subset = (Cell_Type %in% selected_cell_types & STIM == "EFE1"))
            
            
            # Calculate average expression for each cell type combined
            avg_expr <- AverageExpression(
              subset(combined, subset = (Cell_Type %in% selected_cell_types & STIM == "EFE1")),
              features = gene_list,
              group.by = "Cell_Type",
              assays = "RNA",
              slot = "data"
            )$RNA
            
            
            
            # ATAC
            selected_cells <- subset(combined.atac, subset = (Cell_Type %in% selected_cell_types & STIM == "EFE1"))
            
            # ATAC combined.atac
            avg_expr <- AverageExpression(
              subset(combined.atac, subset = (Cell_Type %in% selected_cell_types & STIM == "EFE1")),
              features = gene_list,
              group.by = "Cell_Type",
              assays = "ACTIVITY",
              slot = "data"
            )$ACTIVITY
            
            
            
            
            
            # Subset for selected cell types
            avg_expr_selected <- avg_expr[, selected_cell_types]
            
            # Z-score normalize across cell types (by row/gene)
            avg_expr_z <- t(scale(t(avg_expr_selected)))
            
            
            # Create heatmap with pheatmap
            library(pheatmap)
            
            avg_expr_z_ordered <- avg_expr_z[, c("EC","EC_E1","CM")]
            
            pheatmap(
              avg_expr_z_ordered,
              color = colorRampPalette(c("blue", "yellow","darkred"))(100), # z
              #color = colorRampPalette(c("navy","lightyellow","pink","tomato","red","darkred"))(100), # scaled
              breaks = seq(-2, 2, length.out = 101), # z
              #breaks = seq(0, 1, length.out = 101), # scaled
              cluster_rows = FALSE, # TRUE
              cluster_cols = FALSE,
              show_rownames = TRUE,
              show_colnames = TRUE,
              annotation_names_row = TRUE,
              annotation_names_col = TRUE,
              fontsize_row = 8,
              fontsize_col = 10,
              angle_col = "45",
              main = "",
              display_numbers = FALSE,
              cellwidth = 20,
              cellheight = 12,
              border_color = NA
            )
            
            
            
            # EC_E3
            
            # Subset the Seurat object for selected cell types
            selected_cell_types = c("EC","EC_E3","FB")
            gene_list = c("VWF","FLI1","EGFL7","BTNL9","PIK3R3","ABLIM3","ITGA6","PTPRB","NOTCH4","ELN","COL14A1","COL1A1","ABCA8","C7","SCN7A","CDH19","DCN")
            
            selected_cells <- subset(combined, subset = (Cell_Type %in% selected_cell_types & STIM == "EFE3"))
            
            # Calculate average expression for each cell type
            avg_expr <- AverageExpression(
              subset(combined, subset = (Cell_Type %in% selected_cell_types & STIM == "EFE3")),
              features = gene_list,
              group.by = "Cell_Type",
              assays = "RNA",
              slot = "data"
            )$RNA
            
            
            # ATAC combined.atac
            selected_cells <- subset(combined.atac, subset = (Cell_Type %in% selected_cell_types & STIM == "EFE3"))
            
            # ATAC combined.atac
            avg_expr <- AverageExpression(
              subset(combined.atac, subset = (Cell_Type %in% selected_cell_types & STIM == "EFE3")),
              features = gene_list,
              group.by = "Cell_Type",
              assays = "ACTIVITY",
              slot = "data"
            )$ACTIVITY
            
            
            
            # Subset for selected cell types
            avg_expr_selected <- avg_expr[, selected_cell_types]
            
            # Z-score normalize across cell types (by row/gene)
            avg_expr_z <- t(scale(t(avg_expr_selected)))
            
            # Create heatmap with pheatmap
            library(pheatmap)
            
            avg_expr_z_ordered <- avg_expr_z[, c("EC","EC_E3","FB")]
            
            pheatmap(
              avg_expr_z_ordered,
              color = colorRampPalette(c("blue", "yellow","darkred"))(100), # z
              #color = colorRampPalette(c("navy","lightyellow","pink","tomato","red","darkred"))(100), # scaled
              breaks = seq(-2, 2, length.out = 101), # z
              #breaks = seq(0, 1, length.out = 101), # scaled
              cluster_rows = FALSE, # TRUE
              cluster_cols = FALSE,
              show_rownames = TRUE,
              show_colnames = TRUE,
              annotation_names_row = TRUE,
              annotation_names_col = TRUE,
              fontsize_row = 8,
              fontsize_col = 10,
              angle_col = "45",
              main = "",
              display_numbers = FALSE,
              cellwidth = 20,
              cellheight = 12,
              border_color = NA
            )
            
            
            
            
            # Fig.K
            # load processed
            GO_1 = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/DE/Enrichr_ec13markers/EC_E1_GO_Biological_Process_padj.txt")
            GO_3 = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/DE/Enrichr_ec13markers/EC_E3_GO_Biological_Process_padj.txt")
            
            # Let's calculate gene ratio for GO_1 and GO_3 022126
            q.temp = do.call("rbind", strsplit(GO_1$Overlap,"\\/"))
            GO_1[, GeneRatio := round(as.numeric(q.temp[,1])/as.numeric(q.temp[,2]),digits = 2)] # GO_1
            p.temp = do.call("rbind", strsplit(GO_3$Overlap,"\\/"))
            GO_3[, GeneRatio := round(as.numeric(p.temp[,1])/as.numeric(p.temp[,2]),digits = 2)] # GO_3
            
            # top5 terms for EC_E1
            #GO = GO_1[, group := "EC_E1"][1:5,] # GO_1[grepl("Contr|Muscle|Imp|Stress|Cardi|Heart",Term)][`P-value` < 0.05][c(1,2,4,6),]
            
            GO = rbind(GO_1[, group := "EC_E1"][1:5,][,class:="Top enriched"],GO_1[grepl("Contr|Muscle|Imp|Stress|Cardi|Heart",Term)][`P-value` < 0.05][c(1,2,4,6),][,class:="CM related func."])
            
            GO[Term == "Regulation of Intracellular Signal Transduction (GO:1902531)", 
               Term := "Reg. of Intracellular Sig. Transduction (GO:1902531)"][
                 Term == "Cell Migration Involved in Sprouting Angiogenesis (GO:0002042)", Term := "Cell Migr. Involved in Sprout. Angiogenesis (GO:0002042)"][
                   Term == "Calcium Import Into the Mitochondrion (GO:0036444)", Term := "Calcium Import Into the Mt. (GO:0036444)"]
            
            GO[, `-log10(P-value)` := -log10(`P-value`)]
            
            # lollipop plot
            ggplot(GO, aes(x = `-log10(P-value)`, y = reorder(Term, `-log10(P-value)`))) +
              # Lollipop sticks
              geom_segment(aes(x = 0, xend = `-log10(P-value)`,y = reorder(Term, `-log10(P-value)`), yend = reorder(Term, `-log10(P-value)`)),
                           color = "gray80",linewidth = 1) +
              # Left-aligned text labels
              geom_text(aes(x = 0, y = reorder(Term, `-log10(P-value)`),label = Term),
                        color = "black",size = 5,hjust = 0,vjust = -0.6,nudge_x = -0.1) +
              # Points
              geom_point(aes(size = GeneRatio, color = group),alpha = 0.8,                # "Odds ratio"
                         show.legend = c(color = FALSE, size = TRUE)) +                
              # Scales
              scale_size_continuous(name = "Gene ratio",range = c(3, 6),limits = c(0, 0.3),breaks = c(0, 0.15, 0.3)) + 
              scale_color_manual(name = "Cell type",values = c("EC_E1" = "#00BFFF", "EC_E3" = "#DDA0DD")) +
              labs(x = "-log10 (P-value)",y = NULL) +
              theme_classic() + facet_grid(class~., scales = "free") + 
              theme(panel.grid.major.y = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(size = 14),
                    axis.title.x = element_text(size = 14),axis.text.y = element_blank(),axis.ticks.y = element_blank(),
                    legend.position = c(0.9,0.8),legend.title = element_text(size = 16),legend.text = element_text(size = 16),
                    strip.background = element_blank(),  # Removes facet label background
                    strip.text.y = element_blank(),      # Removes facet text completely #legend.key.size = unit(1, "lines"),
                    plot.margin = margin(l = max(nchar(GO$Term)) * 0.5 + 10))
            
            
            # Fig.M
            
            # top5 terms for EC_E3
            GO = rbind(GO_3[, group := "EC_E3"][1:5,][,class:="Top enriched"],GO_3[grepl("Contr|Muscle|Imp|Stress|Cardi|Heart",Term)][`P-value` < 0.05][c(1,3,9,14),][,class:="Mesc. related func."])
            
            
            GO[Term == "Cell Surface Receptor Protein Tyrosine Kinase Signaling Pathway (GO:0007169)", 
               Term := "Cell Surface Recept. Prot. Tyrosine Kinase Sig. Pathway (GO:0007169)"][
                 Term == "Positive Regulation of Cell Migration (GO:0030335)", Term := "Pos. Reg. of Cell Migration (GO:0030335)"][
                   Term == "Regulation of Cell Migration (GO:0030334)", Term := "Reg. of Cell Migration (GO:0030334)"][
                     Term == "Positive Regulation of Cardiac Epithelial to Mesenchymal Transition (GO:0062043)", Term := "Pos. Reg. of Cardiac Epithelial to Mesenchymal Trans. (GO:0062043)"][
                       Term == "Positive Regulation of Smooth Muscle Cell Migration (GO:0014911)", Term := "Pos. Reg. of Smooth Muscle Cell Migration (GO:0014911)"][
                         Term == "Positive Regulation of Smooth Muscle Cell Proliferation (GO:0048661)", Term := "Pos. Reg. of Smooth Muscle Cell Proliferation (GO:0048661)"]
            
            GO[, `-log10(P-value)` := -log10(`P-value`)]
            
            # lollipop plot
            ggplot(GO, aes(x = `-log10(P-value)`, y = reorder(Term, `-log10(P-value)`))) +
              # Lollipop sticks
              geom_segment(aes(x = 0, xend = `-log10(P-value)`,y = reorder(Term, `-log10(P-value)`), yend = reorder(Term, `-log10(P-value)`)),
                           color = "gray80",linewidth = 0.8) +
              # Left-aligned text labels
              geom_text(aes(x = 0, y = reorder(Term, `-log10(P-value)`),label = Term),
                        color = "black",size = 4,hjust = 0,vjust = -1.5,nudge_x = -0.1) +
              # Points
              geom_point(aes(size = GeneRatio, color = group),alpha = 0.8,             # `Odds Ratio`
                         show.legend = c(color = FALSE, size = TRUE)) +
              # Scales
              scale_size_continuous(name = "Gene ratio",range = c(3, 6),limits = c(0, 0.4),breaks = c(0, 0.2, 0.4)) +
              scale_color_manual(name = "Cell type",values = c("EC_E1" = "#00BFFF", "EC_E3" = "#DDA0DD")) +
              labs(x = "-log10 (P-value)",y = NULL) +
              theme_classic() + facet_grid(class~., scales = "free") + 
              theme(panel.grid.major.y = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(size = 12),
                    axis.title.x = element_text(size = 12),axis.text.y = element_blank(),axis.ticks.y = element_blank(),
                    legend.position = c(0.9, 0.8),legend.title = element_text(size = 12),legend.text = element_text(size = 10),legend.key.size = unit(1, "lines"),
                    strip.background = element_blank(),  # Removes facet label background
                    strip.text.y = element_blank(),      # Removes facet text completely #legend.key.size = unit(1, "lines"),
                    plot.margin = margin(l = max(nchar(GO$Term)) * 0.5 + 10))
            
            
            
            
            
            # Fig.P
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
            
            
            
            
            
            
            # Fig.Q
            
            # Custom smoothing function
            library(motifmatchr)
            library(JASPAR2020)
            library(TFBSTools)
            library(BSgenome.Hsapiens.UCSC.hg38)
            library(Biostrings)
            
            # extract position frequency matrices for the motifs
            pwm <- getMatrixSet(x = JASPAR2020,opts = list(species = 9606, all_versions = FALSE))
            
            # let us do alternatively
            paper.atac = subset(combined.atac, Cell_Type %in% c("EC","EC_E1","EC_E3"))
            
            DefaultAssay(paper.atac) <- "peaks"
            paper.atac <- AddMotifs(paper.atac, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)
            
            paper.atac <- Footprint(
              object = paper.atac,
              motif.name = c("FLI1"),
              genome = BSgenome.Hsapiens.UCSC.hg38
            )
            
            p_papr = PlotFootprint(paper.atac, features = c("FLI1"), split.by = "Cell_Type", group.by = "SAMPLE",normalization = "subtract", show.expected = TRUE)
            
            
            # 1. First inspect the raw plot data
            library(ggplot2)
            library(dplyr)
            library(zoo)
            
            # 1. Extract and prepare plot data
            plot_build <- ggplot2::ggplot_build(p_papr)
            plot_data <- as.data.frame(plot_build$data[[1]])
            
            # 2. Add proper group names (replace with your actual group names)
            sample_names <- levels(as.factor(paper.atac$SAMPLE))  # From group.by
            cell_types <- levels(as.factor(paper.atac$Cell_Type)) # From split.by
            
            # 3. Apply smoothing
            smoothed_data <- plot_data %>%
              mutate(
                SAMPLE = factor(group, labels = sample_names),
                Cell_Type = factor(PANEL, labels = cell_types)
              ) %>%
              group_by(group, PANEL) %>%
              mutate(
                y_smooth = zoo::rollmean(y, k = 3, fill = "extend", align = "center")
              ) %>%
              ungroup()
            
            smoothed_data <- plot_data %>%
              mutate(
                SAMPLE = factor(group, labels = sample_names),
                Cell_Type = factor(PANEL, labels = cell_types)
              ) %>%
              group_by(group, PANEL) %>%
              mutate(
                y_smooth = if(n() >= 13) {  # Only smooth if enough data points
                  zoo::rollmean(y, k = 13, fill = "extend", align = "center")
                } else {
                  y  # Return original if not enough points
                }
              ) %>%
              ungroup()
            
            
            # 4. Create the plot
            sample_colors = condition_colors
            ggplot(smoothed_data[smoothed_data$Cell_Type %in% c("EC") | (smoothed_data$SAMPLE == "EFEs" & smoothed_data$Cell_Type %in% c("EC_E1", "EC_E3")), ], aes(x = x, y = y_smooth, color = SAMPLE)) +
              geom_line(linewidth = 0.8) +
              facet_wrap(~Cell_Type, nrow = 1) +  # One column for each cell type
              scale_color_manual(values = sample_colors, name = "Condition") +
              labs(x = "Position from motif (bp)", 
                   y = "Normalized footprint signal") +
              theme_classic() +
              #xlim(-150, 150) +  # Match PlotFootprint's default range
              theme(legend.position = "top")
            
            
            
            
            # Fig.S
            library(data.table)
            library(ComplexHeatmap)
            library(circlize)
            library(foreach)
            library(edgeR)
            library(Seurat)
            
            # Method 1: Using ComplexHeatmap (recommended)
            create_gene_expression_heatmap <- function(dt.cpmed, 
                                                       gene_column = "gene",
                                                       output_file = "gene_expression_heatmap.pdf",
                                                       show_gene_names = TRUE,
                                                       top_n_genes = NULL,
                                                       width = 12,
                                                       height = 16,
                                                       title = "Gene Expression Heatmap (Z-scores)",
                                                       cluster_columns = FALSE,
                                                       cluster_rows = TRUE,
                                                       row_font_size = 6,
                                                       color_palette = "blue-white-red",
                                                       z_cap = 3,
                                                       custom_colors = NULL,
                                                       color_breaks = NULL,       # New: specify exact break points
                                                       color_values = NULL,       # New: specify exact color values
                                                       n_color_steps = 100,       # New: number of color interpolation steps
                                                       symmetric_colors = TRUE) { # New: force symmetric color distribution
              
              # Load required libraries
              library(ComplexHeatmap)
              library(circlize)
              
              # Store original gene names from the specified column
              original_genes <- dt.cpmed[[gene_column]]
              
              # Remove gene column to get expression matrix
              expr_cols <- setdiff(names(dt.cpmed), gene_column)
              
              # Get the original column order from dt.cpmed (excluding gene column)
              original_column_order <- expr_cols
              
              # Create expression matrix
              expr_matrix <- as.matrix(dt.cpmed[, ..expr_cols])
              
              # Keep track of original gene names as row names
              rownames(expr_matrix) <- original_genes
              
              # Calculate z-scores for each gene (row-wise)
              cat("Calculating z-scores...\n")
              z_matrix <- t(scale(t(expr_matrix)))
              
              # Check for and handle NA values (genes with zero variance)
              na_rows <- apply(z_matrix, 1, function(x) any(is.na(x)))
              if (any(na_rows)) {
                cat("Removing", sum(na_rows), "genes with zero variance\n")
                z_matrix <- z_matrix[!na_rows, ]
                original_genes <- original_genes[!na_rows]
              }
              
              # Always set row names to gene names from the specified column
              rownames(z_matrix) <- original_genes
              
              # Cap extreme z-scores for better visualization
              z_matrix[z_matrix > z_cap] <- z_cap
              z_matrix[z_matrix < -z_cap] <- -z_cap
              
              # METHOD 1: Custom color breaks and values (most flexible)
              if (!is.null(color_breaks) && !is.null(color_values)) {
                if (length(color_breaks) != length(color_values)) {
                  stop("color_breaks and color_values must have the same length")
                }
                cat("Using custom color breaks and values\n")
                col_fun <- colorRamp2(color_breaks, color_values)
                
                # METHOD 2: Custom colors with automatic symmetric distribution
              } else if (!is.null(custom_colors)) {
                cat("Using custom colors with", length(custom_colors), "colors\n")
                
                if (symmetric_colors && length(custom_colors) %% 2 == 1) {
                  # For odd number of colors, middle color is at 0
                  n_colors <- length(custom_colors)
                  mid_index <- ceiling(n_colors / 2)
                  
                  # Create symmetric breaks
                  if (n_colors == 3) {
                    breaks <- c(-z_cap, 0, z_cap)
                  } else {
                    # For 5 colors: -z_cap, -z_cap/2, 0, z_cap/2, z_cap
                    # For 7 colors: -z_cap, -2*z_cap/3, -z_cap/3, 0, z_cap/3, 2*z_cap/3, z_cap
                    breaks <- seq(-z_cap, z_cap, length.out = n_colors)
                  }
                  col_fun <- colorRamp2(breaks, custom_colors)
                  
                } else {
                  # Asymmetric or even number of colors
                  if (length(custom_colors) == 2) {
                    # Two colors: low to high
                    col_fun <- colorRamp2(c(-z_cap, z_cap), custom_colors)
                  } else if (length(custom_colors) == 3) {
                    # Three colors: low, middle, high
                    col_fun <- colorRamp2(c(-z_cap, 0, z_cap), custom_colors)
                  } else {
                    # Multiple colors with equal spacing
                    breaks <- seq(-z_cap, z_cap, length.out = length(custom_colors))
                    col_fun <- colorRamp2(breaks, custom_colors)
                  }
                }
                
                # METHOD 3: Predefined palettes with multiple color points
              } else {
                cat("Using color palette:", color_palette, "\n")
                
                col_fun <- switch(color_palette,
                                  
                                  # Simple 3-color palettes
                                  "blue-white-red" = colorRamp2(c(-z_cap, 0, z_cap), c("navy", "white", "red")),
                                  "blue-yellow-red" = colorRamp2(c(-z_cap, 0, z_cap), c("darkblue", "yellow", "darkred")),
                                  "lightblue-white-red" = colorRamp2(c(-z_cap, 0, z_cap), c("lightblue", "white", "red")),
                                  "purple-white-green" = colorRamp2(c(-z_cap, 0, z_cap), c("purple4", "white", "green4")),
                                  "magenta-white-green" = colorRamp2(c(-z_cap, 0, z_cap), c("magenta4", "white", "darkgreen")),
                                  
                                  # Fallback to default
                                  colorRamp2(c(-z_cap, 0, z_cap), c("navy", "white", "red"))
                )
              }
              
              # Determine font size based on number of genes if not specified
              if (missing(row_font_size)) {
                row_font_size <- if (nrow(z_matrix) <= 50) {
                  8
                } else if (nrow(z_matrix) <= 100) {
                  6
                } else if (nrow(z_matrix) <= 200) {
                  5
                } else {
                  4  # Very small font for many genes
                }
              }
              
              # Create the heatmap
              cat("Creating heatmap with gene labels...\n")
              ht <- Heatmap(z_matrix,
                            name = "Z-score",
                            col = col_fun,
                            
                            # Row settings (genes) - ALWAYS show gene names
                            row_title = paste0("Genes (", nrow(z_matrix), ")"),
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = row_font_size),
                            cluster_rows = cluster_rows,
                            row_dend_width = unit(3, "cm"),
                            
                            # Column settings (samples) - preserve original order
                            column_title = title,
                            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                            show_column_names = TRUE,
                            column_names_gp = gpar(fontsize = 10),
                            cluster_columns = cluster_columns,
                            
                            # Preserve column order when not clustering
                            column_order = if (!cluster_columns) seq_len(ncol(z_matrix)) else NULL,
                            
                            # Heatmap legend - dynamic based on color breaks
                            heatmap_legend_param = list(
                              title = "Z-score",
                              title_gp = gpar(fontsize = 12),
                              labels_gp = gpar(fontsize = 10),
                              legend_height = unit(4, "cm"),
                              at = if (!is.null(color_breaks)) {
                                color_breaks
                              } else {
                                c(-z_cap, -z_cap/2, 0, z_cap/2, z_cap)
                              },
                              labels = if (!is.null(color_breaks)) {
                                as.character(round(color_breaks, 2))
                              } else {
                                c(paste0("≤", -z_cap), -z_cap/2, "0", z_cap/2, paste0("≥", z_cap))
                              }
                            ),
                            
                            # Performance optimization
                            use_raster = nrow(z_matrix) > 1000,
                            raster_quality = 2)
              
              # Save to PDF
              draw(ht)
              dev.off()
              
              # Return information
              return(list(heatmap = ht, 
                          z_matrix = z_matrix, 
                          genes = rownames(z_matrix), 
                          samples = colnames(z_matrix),
                          column_order_preserved = !cluster_columns,
                          color_palette = color_palette,
                          z_cap = z_cap,
                          color_function = col_fun))
            }
            
            # Usage examples:
            dt.cpmed = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/slides/Paper/Circulation/revision/data/Bulk_allsample_count_123126.txt")
            
            # Example 1: Basic heatmap
            gene_interest = c("FLI1","LINGO1","LMOD1","TAGLN",
                              "TCF4","ELN","PLXND1",
                              "FGF12","SORBS2","MEF2C","MEF2D",
                              "NEGR1","VIM",
                              "MYH11")
            
            # order labels
            dt.cpmed = dt.cpmed[,c("gene","MOCK1","MOCK2","VEGF1","VEGF2","FGF1","FGF2","SEMA1","SEMA2"), with = F]
            
            # Hmap
            dt.plot = dt.cpmed[gene %in% gene_interest]
            dt.plot$gene <- factor(dt.plot$gene, levels = gene_interest)
            dt.plot <- dt.plot[order(dt.plot$gene), ]
            
            create_gene_expression_heatmap(dt.plot,                                       # dt.cpmed[gene %in% gene_interest]
                                           color_breaks = c(-2,0,2),
                                           gene_column = "gene",
                                           show_gene_names = TRUE,  # Changed default to TRUE
                                           output_file = "/Users/yangyu/Desktop/my_heatmap.pdf",
                                           cluster_columns = FALSE,
                                           cluster_rows = FALSE,
                                           color_values = c("navy","white","red"))
            
            
            
            
            
            # Fig.X
            pred_label = fread("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/Merfish/processed/files/merged_predlabel_5samps_070825.csv")
            
            # select cols
            pred_label = pred_label[,c(36,41:55),with=F]
            
            pred_label[,SAMPLE:="CTRLs"][grepl("efe",library_id),SAMPLE:="EFEs"]
            target_cols = colnames(pred_label)[4:15]
            
            # Add a column for the name of the column with the largest value
            pred_label[, `:=`(
              largest_col = names(.SD)[max.col(.SD)],  # Column name with the largest value
              largest_value = do.call(pmax, .SD),      # Largest value
              second_largest_col = apply(.SD, 1, function(row) {
                sorted_indices <- order(row, decreasing = TRUE)
                names(.SD)[sorted_indices[2]]  # Column name with the second largest value
              }),
              second_largest_value = apply(.SD, 1, function(row) {
                sorted_values <- sort(row, decreasing = TRUE)
                sorted_values[2]  # Second largest value
              }),
              diff_largest_second = apply(.SD, 1, function(row) {
                sorted_values <- sort(row, decreasing = TRUE)
                sorted_values[1] - sorted_values[2]  # Difference between largest and second largest
              })
            ), .SDcols = target_cols]
            
            # select defining cell type prob.
            library(speckle)
            
            prop = as.data.table(propeller(clusters = pred_label$largest_col, sample = pred_label$library_id, group = pred_label$SAMPLE),keep.rownames = T)
            
            # convert to plottable format
            setnames(prop,c("P.Value"),c("Pvalue"))
            
            pred_label_long <- as.data.table(melt(prop, id.vars = c("rn", "Pvalue"), 
                                                  measure.vars = c("PropMean.CTRLs", "PropMean.EFEs"),
                                                  variable.name = "Condition", value.name = "Proportion"))
            
            # Define color mapping for conditions
            condition_colors <- c("CTRLs" = "darkgreen",   # Medium sea green for CTRLs
                                  "EFEs" = "#E41A1C")    # Red for EFEs
            
            # Clean condition names (if not already done)
            pred_label_long$Condition <- gsub("PropMean.", "", pred_label_long$Condition)
            pred_label_long$Condition <- factor(pred_label_long$Condition, levels = c("EFEs","CTRLs"))
            
            # v2 use this
            prop_stats <- pred_label[, .(count = .N), by = .(SAMPLE, sample_id, largest_col)][
              , proportion := count/sum(count), by = .(SAMPLE, sample_id)][
                , .(mean_prop = mean(proportion),
                    se = sd(proportion)/sqrt(.N),  # Standard error
                    individual_prop = proportion),  # Keep original points
                by = .(SAMPLE, largest_col)]
            
            # reorder the cell label
            prop_stats$SAMPLE <- factor(prop_stats$SAMPLE, levels = c("EFEs","CTRLs"))
            prop_stats = prop_stats[grepl("EC|EndoC|CM|FB|SMC", largest_col)]
            prop_stats$largest_col <- factor(prop_stats$largest_col, levels = c("EC_E1","EC_E3", "EC","CM","FB","SMC","EndoC"))
            
            
            # plot
            prop_stats$largest_col <- factor(prop_stats$largest_col, levels = c("EC_E3", "EC_E1", "EC","CM","EndoC", "FB","SMC"))
            
            ggplot() +
              geom_bar(data = prop_stats[grepl("EC|EndoC|CM|FB|SMC", largest_col)], 
                       aes(x = SAMPLE, y = mean_prop, fill = SAMPLE), 
                       stat = "identity", 
                       position = position_dodge(width = 0.1), 
                       width = 0.3,
                       alpha = 0.5) +  # Added transparency to bars
              
              geom_point(data = prop_stats[grepl("EC|EndoC|CM|FB|SMC", largest_col)], 
                         aes(x = SAMPLE, y = individual_prop, fill = SAMPLE),  # Use fill instead of color
                         position = position_jitterdodge(
                           jitter.height = 0,
                           dodge.width = 0.1,
                           jitter.width = 0.29
                         ),
                         size = 2, alpha = 1, shape = 21, color = "black", stroke = 1) +  # shape 21 allows fill + border
              labs(title = "", x = "", y = "Cell proportion") + 
              facet_wrap(~largest_col, scales = "free_y", nrow = 1) +
              scale_fill_manual(values = condition_colors) +
              scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
              theme_classic() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                    axis.text.y = element_text(size = 12),
                    axis.title.y = element_text(size = 14),
                    axis.title.x = element_blank(),
                    legend.title = element_blank(),
                    legend.position = "none",
                    strip.background = element_blank(),strip.text = element_text(size = 10, face = "bold"),
                    panel.spacing = unit(1, "lines"))
            
            
            
            
            # Fig.Y
            PATH = "/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/Merfish/processed/files/"
            files = data.table(files = list.files(PATH, all.files = F))
            
            files = files[grepl("nenrichment_070825",files)]
            key_matc = data.table(V1 = c(0:11), pred_celltype = c("CM","EC","EC_E1","EC_E3","EndoC","FB","Macrophage","Mast-cell","Neuro","PeriC","SMC","T-cell"))
            
            
            
            # rescale unused
            rescale_to_minus_one_to_one <- function(x) {
              max_abs <- max(abs(x))  # Find the maximum absolute value
              x / max_abs  # Rescale each value
            }
            
            enrichment = foreach(fl = unique(files$files),.combine = "rbind")%do%{
              tmps = fread(paste0(PATH,fl), header = T)
              setnames(tmps, c("V1","CM","EC","EC_E1","EC_E3","EndoC","FB","Macrophage","Mast-cell","Neuro","PeriC","SMC","T-cell"))
              
              tmps = merge(setkey(tmps, V1), setkey(key_matc, V1))
              tmps[, sample := paste(strsplit(fl,"_")[[1]][1],strsplit(fl,"_")[[1]][2],sep="_")]
              
              tmps$V1 <- NULL
              tmps_melt = as.data.table(melt(tmps, id.vars = c("sample","pred_celltype")))
              
              tmps_melt[is.na(value), value := 0]
              tmps_melt$value_resc <- rescale_to_minus_one_to_one(tmps_melt$value)
              
              tmps_melt
            }
            
            
            # Reshape the data.table
            enrichment$value_resc <- as.numeric(enrichment$value_resc)
            enrichment$value_resc <- round(enrichment$value_resc, digits = 4)
            
            # before reshaping the data, average across SAMPLE
            enrichment[, SAMPLE := "CTRLs"][grepl("efe",sample), SAMPLE := "EFEs"]
            
            # using stouffers method and recompute pvalue using fishers methods
            enrichment[, combi := paste(pred_celltype, variable, sep = "_")]
            
            our_enrichment = foreach(sp = unique(enrichment$SAMPLE),.combine = "rbind")%do%{
              tmps = enrichment[SAMPLE == sp]
              pred_combi = foreach(cb = unique(tmps$combi),.combine = "rbind")%do%{
                tmps3 = tmps[combi == cb]
                z_scores = tmps3$value
                combined_z <- sum(z_scores) / sqrt(length(z_scores))
                combined_p <- 2 * pnorm(-abs(combined_z))  # Two-tailed p-value
                chi_sq <- sum(z_scores^2)
                combined_p_fisher <- round(pchisq(chi_sq, df = length(z_scores), lower.tail = FALSE),digits = 4)
                tmps3[, stouffer_p := combined_p][, stouffer_z := combined_z][, fisher_p := combined_p_fisher]
                tmps3
              }
            }
            
            # use pvalue as stouffer p
            enrichment = our_enrichment[,list(value = unique(stouffer_z), p = unique(fisher_p), pvalue = unique(stouffer_p)), list(SAMPLE, pred_celltype, variable, combi)]
            enrichment$value <- round(enrichment$value, digits = 4)
            
            # start to prepare for plot
            
            # enrichment matrix change labels as EFEs or CTRLs
            dt_wide <- dcast(enrichment[grepl("EFEs",SAMPLE)], pred_celltype ~ variable, value.var = "value", fun.aggregate = mean)
            
            # pvalue matrix
            pval_matrix = dcast(enrichment[grepl("EFEs",SAMPLE)], pred_celltype ~ variable, value.var = "pvalue", fun.aggregate = mean)
            
            # Convert p matrix to matrix
            p_mat <- as.matrix(pval_matrix[, -1])
            p_mat[is.na(p_mat)] <- 0
            
            str(p_mat)
            rownames(p_mat) <- pval_matrix$pred_celltype
            
            # Define the order of rows and columns
            row_order <- c("EC","EC_E1","EC_E3","EndoC","T-cell","FB","Macrophage","Mast-cell","PeriC","Neuro","SMC","CM")  # Replace with your desired row order
            col_order <- c("EC","EC_E1","EC_E3","EndoC","T-cell","FB","Macrophage","Mast-cell","PeriC","Neuro","SMC","CM")  # Replace with your desired column order
            
            # Reorder the matrix
            pmat_ordered <- p_mat[row_order, col_order]
            
            # only show lower or upper
            pmat_lower <- pmat_ordered
            pmat_lower[upper.tri(pmat_lower)] <- NA  # Set upper triangle to NA (including diagonal)
            
            # Get row and column names
            col_names <- colnames(pmat_lower)
            row_names <- rownames(pmat_lower)
            
            # Convert value matrix to matrix
            mat <- as.matrix(dt_wide[, -1])
            mat[is.na(mat)] <- 0
            
            str(mat)
            rownames(mat) <- dt_wide$pred_celltype
            
            # Reorder the matrix
            mat_ordered <- mat[row_order, col_order]
            
            # only show lower or upper
            mat_lower <- mat_ordered
            mat_lower[upper.tri(mat_lower)] <- NA  # Set upper triangle to NA (including diagonal)
            
            # Get row and column names
            col_names <- colnames(mat_lower)
            row_names <- rownames(mat_lower)
            
            
            # Define value ranges and corresponding color segments
            value_ranges <- list(
              c(-100, -50),    # Extreme negative values
              c(-50, -3),    # Extreme negative values
              c(-3, 3),       # Main range of interest
              c(3, 50),       # Extreme positive values
              c(50, 100)       # Extreme positive values
            )
            
            # Define corresponding color palettes for each range
            color_palettes <- list(
              colorRampPalette(c("#6A5ACD", "blue"))(100),          # -100 to -2
              colorRampPalette(c("blue","skyblue", "lightblue"))(100),          # -100 to -2
              colorRampPalette(c("lightblue","white","lightyellow"))(20),  # -2 to 2
              colorRampPalette(c("orange","tomato"))(100),  # 2 to 7 (highlighted)
              colorRampPalette(c("tomato","orangered","red"))(100)  # 2 to 7 (highlighted)
              #colorRampPalette(c("orangered", "tomato", "red"))(200)               # 7 to 100
            )
            
            # Combine all color segments
            heatmap_colors <- unlist(color_palettes)
            
            # Create non-linear breaks that match the value ranges
            color_breaks <- c(
              seq(value_ranges[[1]][1], value_ranges[[1]][2], length.out = 101)[-1],
              seq(value_ranges[[2]][1], value_ranges[[2]][2], length.out = 101)[-1],
              seq(value_ranges[[3]][1], value_ranges[[3]][2], length.out = 21)[-1], # Remove duplicate point
              seq(value_ranges[[4]][1], value_ranges[[4]][2], length.out = 101)[-1],
              seq(value_ranges[[5]][1], value_ranges[[5]][2], length.out = 101)[-1]
            )
            
            # Plot the heatmap
            # First, convert p-values to asterisks
            signif_stars <- matrix("", nrow = nrow(pmat_lower), ncol = ncol(pmat_lower),dimnames = dimnames(pmat_lower))
            
            # Define significance levels (adjust thresholds as needed)
            signif_stars[pmat_lower < 0.00001] <- "*"
            
            # Now create the heatmap with asterisks
            pheatmap(mat_lower, 
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     color = heatmap_colors,
                     breaks = color_breaks,
                     main = "",
                     angle_col = "90",
                     fontsize_col = 14,
                     fontsize_row = 14,
                     na_col = "white",
                     gaps_col = 1:ncol(mat_lower),
                     gaps_row = 1:nrow(mat_lower),
                     border_color = NA,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     legend = TRUE,
                     display_numbers = signif_stars,  # Use the asterisk matrix instead of p-values
                     number_color = "black",
                     fontsize_number = 16,
                     labels_row = rownames(mat_lower),
                     labels_col = colnames(mat_lower),
                     cellwidth = 15,
                     cellheight = 15)
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
