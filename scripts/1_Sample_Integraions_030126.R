#### library  # integrates 3 EFE and 3 Pu Lab Crtl and 3 Young and 3 Fetal samples included atac data
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(harmony)
library(ArchR)
library(hexbin)
library(biovizBase)


#setseed(1234)
            addArchRThreads(threads = 32) 
            
            #### functions
            
            # removes doublet using DoubletFinder
            Remove_Doublet=function(seurat_objt=EFE001, pseu_rate = 0.04){
              
              DefaultAssay(seurat_objt) <- "RNA"
              
              # standard
              seurat_objt <- NormalizeData(seurat_objt)
              seurat_objt <- FindVariableFeatures(seurat_objt, selection.method = "vst", nfeatures = 2000)
              seurat_objt <- ScaleData(seurat_objt)
              seurat_objt <- RunPCA(seurat_objt) # ElbowPlot(seurat_objt)
              seurat_objt <- RunUMAP(seurat_objt, dims = 1:20)
              
              sweep.res.list_seut <- paramSweep_v3(seurat_objt, PCs = 1:20, sct = TRUE)
              sweep.stats_seut <- summarizeSweep(sweep.res.list_seut, GT = FALSE)
              bcmvn_seut <- find.pK(sweep.stats_seut)
              
              # script for visualizaion of pK parameter selection
              pK=as.numeric(as.character(bcmvn_seut$pK))
              BCmetric=bcmvn_seut$BCmetric
              pK_choose = pK[which(BCmetric %in% max(BCmetric))]
              
              nExp_poi <- round(pseu_rate*nrow(seurat_objt@meta.data))  ## Assuming 4% doublet formation rate - tailor for your dataset
              seurat_objt <- doubletFinder_v3(seurat_objt, pN = 0.25, pK = pK_choose, nExp = nExp_poi, PCs = 1:20) ##update newly
              
              # plot
              DF.name = colnames(seurat_objt@meta.data)[grepl("DF.classification", colnames(seurat_objt@meta.data))]
              plts = cowplot::plot_grid(ncol = 1, DimPlot(seurat_objt, group.by = DF.name) + NoAxes())
              
              
              # remove doublets
              seurat_objt=seurat_objt[, seurat_objt@meta.data[, DF.name] == "Singlet"]
              return(list(objt = seurat_objt, plts = plts))
            }
            
            # removes Mt and Rb genes from the seurat objects
            Remove_MtRb=function(seurat_objt=EFE001){
              
              mito.genes <- grep(pattern = "^MT-", x = rownames(seurat_objt), value = FALSE)
              ribo.genes <- grep(pattern = "^RP[SL]", x = rownames(seurat_objt), value = FALSE)
              retained <- rownames(seurat_objt)[-c(mito.genes, ribo.genes) ]
              seurat_objt[["RNA"]] <- subset(seurat_objt[["RNA"]], features = retained)
              
              return(seurat_objt)
            }
            
            # reads in disease RNA and ATAC data
            Read_Disease=function(counts = counts_efe1, fragpath = fragpath_efe1, projname = "EFE001", annt=annotation){
              dat_BU = CreateSeuratObject(
                counts = counts$`Gene Expression`,
                assay = "RNA",
                project=projname
              )
              
              dat_BU[["ATAC"]] <- CreateChromatinAssay(
                counts = counts$Peaks,
                sep = c(":", "-"),
                fragments = fragpath,
                annotation = annt
              )
              
              DefaultAssay(dat_BU) <- "ATAC"
              dat_BU <- NucleosomeSignal(dat_BU)
              dat_BU <- TSSEnrichment(dat_BU)
              
              DefaultAssay(dat_BU) <- "RNA"
              dat_BU[["percent.mt"]] <- PercentageFeatureSet(dat_BU, pattern = "^MT-")
              dat_BU[['percent.ribo']] <- PercentageFeatureSet(dat_BU, pattern = "^RP[SL]")
              
              return(dat_BU)
            }
            
            # reads in control RNA data
            Read_Control=function(path = y1, projname = "young1"){
              ctrl=Read10X(data.dir = path)
              ctrl_BU <- CreateSeuratObject(counts = ctrl,min.cells = 3, min.features = 200,project = projname)
              
              ctrl_BU[["percent.mt"]] <- PercentageFeatureSet(ctrl_BU, pattern = "^MT-")
              ctrl_BU[['percent.ribo']] <- PercentageFeatureSet(ctrl_BU, pattern = "^RP[SL]")
              
              return(ctrl_BU)
            }
            
            # process of macs2 called peaks for each sample stored in the same folder
            Peak_Macs2=function(macs.path = macs_path, macs.name = name, peaks.macs2 = peaks_macs2){
              
              
              # convert MACS2 narrowPeaks to Granges object
              
              
              for(i in 1:length(macs.name)){
                df <- read.table(file = paste0(macs.path,macs.name[i], "_peaks.narrowPeak"), 
                                 col.names = c("chr","start", "end", "name", "score", "strand", "fold_change", 
                                               "neg_log10pvalue_summit", "neg_log10qvalue_summit","relative_summit_position"))
                gr <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE, 
                                               starts.in.df.are.0based = TRUE)
                
                peaks.macs2[[macs.name[i]]] <- gr
              }
              
              #peaks.macs2
              
              # remove peaks on nonstandard chromosomes
              peaks.macs2 <- lapply(peaks.macs2, function(x){
                x <- keepStandardChromosomes(x, pruning.mode = "coarse")
              })
              
              # encode a genome blacklist regions for hg38 using the reference list from https://github.com/Boyle-Lab/Blacklist/
              blacklist_df <- fread(paste0(macs.path,"hg38-blacklist.v2.bed"))
              blacklist <- makeGRangesFromDataFrame(blacklist_df, ignore.strand = T, seqnames.field = "V1", start.field = "V2", end.field = "V3")
              
              seqlevelsStyle(blacklist) <- 'UCSC'
              
              # remove peaks in genomic blacklist regions
              peaks.macs2 <- lapply(peaks.macs2, function(x){
                x <- subsetByOverlaps(x, ranges = blacklist, invert = TRUE)
              })
              
              
              return(peaks.macs2)
            }
            
            # creating peaks assay for the seurat objects using a combined peak set from processed macs2
            Create_Peaks=function(seurat_objt=EFE001, fragpath = fragpath_efe1, comb.peaks = combined.peaks, annt=annotation){
              
              DefaultAssay(seurat_objt) <- "ATAC"
              
              # quantify counts in each peak
              macs2_counts <- FeatureMatrix(
                fragments = Fragments(seurat_objt),
                features = comb.peaks,
                cells = colnames(seurat_objt)
              )
              
              # create a new assay using the MACS2 peak set and add it to the Seurat object
              seurat_objt[["peaks"]] <- CreateChromatinAssay(
                counts = macs2_counts,
                fragments = fragpath,
                annotation = annt
              )
              
              
              return(seurat_objt)
            }
            
            #### This is preprocessing -----------------------------------------------------------------------------------------------------------------------
            
            #### 1. load the RNA and ATAC data
            
            # get gene annotations for hg38
            annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
            
            ###change the internet to guest, NOT TCH;  very important
            seqlevelsStyle(annotation) <- "UCSC"
            
            #### EFE samples
            counts_efe1 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE001/filtered_feature_bc_matrix.h5")
            fragpath_efe1 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE001/atac_fragments.tsv.gz"
            
            counts_efe2 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE002/filtered_feature_bc_matrix.h5")
            fragpath_efe2 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE002/atac_fragments.tsv.gz"
            
            counts_efe3 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE003/filtered_feature_bc_matrix.h5")
            fragpath_efe3 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE003/atac_fragments.tsv.gz"
            
            #### PuCtrl samples
            counts_puctrl1 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl001/filtered_feature_bc_matrix.h5")
            fragpath_puctrl1 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl001/atac_fragments.tsv.gz"
            
            counts_puctrl2 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl002/filtered_feature_bc_matrix.h5")
            fragpath_puctrl2 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl002/atac_fragments.tsv.gz"
            
            counts_puctrl3 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl003/filtered_feature_bc_matrix.h5")
            fragpath_puctrl3 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl003/atac_fragments.tsv.gz"
            
            # ctrl samples
            y1="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Ctrl_data/GSM4742860_Young1_processed/"
            y2="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Ctrl_data/GSM4742861_Young2_processed/"
            y3="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Ctrl_data/GSM4742862_Young3_processed/"
            
            i1="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Ctrl_data/GSM4742857_Fetal1_processed/"
            i2="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Ctrl_data/GSM4742858_Fetal2_processed/"
            i3="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Ctrl_data/GSM4742859_Fetal3_processed/"
            
            
            # load the RNA and ATAC data
            EFE001_BU = Read_Disease(counts = counts_efe1, fragpath = fragpath_efe1, projname = "EFE001",annt=annotation)
            EFE002_BU = Read_Disease(counts = counts_efe2, fragpath = fragpath_efe2, projname = "EFE002",annt=annotation)
            EFE003_BU = Read_Disease(counts = counts_efe3, fragpath = fragpath_efe3, projname = "EFE003",annt=annotation)
            
            Ctrl001_BU = Read_Disease(counts = counts_puctrl1, fragpath = fragpath_puctrl1, projname = "Ctrl001",annt=annotation)
            Ctrl002_BU = Read_Disease(counts = counts_puctrl2, fragpath = fragpath_puctrl2, projname = "Ctrl002",annt=annotation)
            Ctrl003_BU = Read_Disease(counts = counts_puctrl3, fragpath = fragpath_puctrl3, projname = "Ctrl003",annt=annotation)
            
            young1_BU = Read_Control(path = y1, projname = "young1")
            young2_BU = Read_Control(path = y2, projname = "young2")
            young3_BU = Read_Control(path = y3, projname = "young3")
            
            fetal1_BU = Read_Control(path = i1, projname = "fetal1")
            fetal2_BU = Read_Control(path = i2, projname = "fetal2")
            fetal3_BU = Read_Control(path = i3, projname = "fetal3")
            
            #### 2. Quality control including filter mt and rb perc and removing mt and rb genes removing doublets ---------------------------------------
            
            # first assign backup objs to active objs
            EFE001 <- EFE001_BU
            EFE002 <- EFE002_BU
            EFE003 <- EFE003_BU
            Ctrl001 <- Ctrl001_BU
            Ctrl002 <- Ctrl002_BU
            Ctrl003 <- Ctrl003_BU
            young1 <- young1_BU
            young2 <- young2_BU
            young3 <- young3_BU
            fetal1 <- fetal1_BU
            fetal2 <- fetal2_BU
            fetal3 <- fetal3_BU
            
            
            # 2.1 .........................................................................................................................................
            
            # RNA
            
            VlnPlot(
              object = EFE001,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo","nCount_ATAC","nFeature_ATAC","TSS.enrichment", "nucleosome_signal"),
              ncol = 8,
              log = TRUE,
              pt.size = 0
              #y.max = c(15000, 10000, 10000, 10000, 10000, 10000, 10, 0.8)
            ) # nCount_RNA: 157-11579; nFeature_RNA: 143-4273; pert.mt: 0.19%-25.2%; pert.rb: 0%-1.3%; Enr: 1.5-6.8; sig: 0.15-0.35
            
            VlnPlot(
              object = EFE002,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo","nCount_ATAC","nFeature_ATAC","TSS.enrichment", "nucleosome_signal"),
              ncol = 8,
              log = TRUE,
              pt.size = 0
              #y.max = c(15000, 10000, 10000, 10000, 10, 0.8)
            ) # nCount_RNA: 281-16315; nFeature_RNA: 221-4741; pert.mt: 0.22%-14.9%; pert.rb: 0.21%-3.72%; Enr: 1.7-6.4; sig: 0.02-0.98
            
            VlnPlot(
              object = EFE003,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo","nCount_ATAC","nFeature_ATAC","TSS.enrichment", "nucleosome_signal"),
              ncol = 8,
              log = TRUE,
              pt.size = 0
              #y.max = c(15000, 10000, 10000, 10000, 10, 0.8)
            ) # nCount_RNA: 281-16315; nFeature_RNA: 221-4741; pert.mt: 0.22%-14.9%; pert.rb: 0.21%-3.72%; Enr: 1.7-6.4; sig: 0.02-0.98
            
            VlnPlot(
              object = Ctrl001,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo","nCount_ATAC","nFeature_ATAC","TSS.enrichment", "nucleosome_signal"),
              ncol = 8,
              log = TRUE,
              pt.size = 0
              #y.max = c(15000, 10000, 10000, 10000, 10, 0.8)
            ) # nCount_RNA: 157-11579; nFeature_RNA: 143-4273; pert.mt: 0.19%-25.2%; pert.rb: 0%-1.3%; Enr: 1.5-6.8; sig: 0.15-0.35
            
            VlnPlot(
              object = Ctrl002,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo","nCount_ATAC","nFeature_ATAC","TSS.enrichment", "nucleosome_signal"),
              ncol = 8,
              log = TRUE,
              pt.size = 0
              #y.max = c(15000, 10000, 10000, 10000, 10, 0.8)
            ) # nCount_RNA: 281-16315; nFeature_RNA: 221-4741; pert.mt: 0.22%-14.9%; pert.rb: 0.21%-3.72%; Enr: 1.7-6.4; sig: 0.02-0.98
            
            VlnPlot(
              object = Ctrl003,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo","nCount_ATAC","nFeature_ATAC","TSS.enrichment", "nucleosome_signal"),
              ncol = 8,
              log = TRUE,
              pt.size = 0
              #y.max = c(15000, 10000, 10000, 10000, 10, 0.8)
            ) # nCount_RNA: 281-16315; nFeature_RNA: 221-4741; pert.mt: 0.22%-14.9%; pert.rb: 0.21%-3.72%; Enr: 1.7-6.4; sig: 0.02-0.98
            
            VlnPlot(
              object = fetal1,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo"),
              ncol = 4,
              log = TRUE,
              pt.size = 0
            ) # nCount_RNA: 1876-22322; nFeature_RNA: 1116-5330; pert.mt: 0.03%-5.4%; pert.rb: 0.21%-2.33%
            
            VlnPlot(
              object = fetal2,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo"),
              ncol = 4,
              log = TRUE,
              pt.size = 0
            ) # nCount_RNA: 796-36090; nFeature_RNA: 544-6183; pert.mt: 0.01%-23.1%; pert.rb: 0.11%-7.66%
            
            VlnPlot(
              object = fetal3,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo"),
              ncol = 4,
              log = TRUE,
              pt.size = 0
            ) # nCount_RNA: 559-27604; nFeature_RNA: 489-5623; pert.mt: 0.03%-6.4%; pert.rb: 0.15%-3.15%
            
            VlnPlot(
              object = young1,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo"),
              ncol = 4,
              log = TRUE,
              pt.size = 0
            ) # nCount_RNA: 1876-22322; nFeature_RNA: 1116-5330; pert.mt: 0.03%-5.4%; pert.rb: 0.21%-2.33%
            
            VlnPlot(
              object = young2,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo"),
              ncol = 4,
              log = TRUE,
              pt.size = 0
            ) # nCount_RNA: 796-36090; nFeature_RNA: 544-6183; pert.mt: 0.01%-23.1%; pert.rb: 0.11%-7.66%
            
            VlnPlot(
              object = young3,
              features = c("nCount_RNA","nFeature_RNA","percent.mt", "percent.ribo"),
              ncol = 4,
              log = TRUE,
              pt.size = 0
            ) # nCount_RNA: 559-27604; nFeature_RNA: 489-5623; pert.mt: 0.03%-6.4%; pert.rb: 0.15%-3.15%
            
            
            # 2.2.1 ......................................................................................................................................... 
            #### removing mt and rb genes 
            EFE001 = Remove_MtRb(seurat_objt = EFE001)
            EFE002 = Remove_MtRb(seurat_objt = EFE002)
            EFE003 = Remove_MtRb(seurat_objt = EFE003)
            Ctrl001 = Remove_MtRb(seurat_objt = Ctrl001)
            Ctrl002 = Remove_MtRb(seurat_objt = Ctrl002)
            Ctrl003 = Remove_MtRb(seurat_objt = Ctrl003)
            young1 = Remove_MtRb(seurat_objt = young1)
            young2 = Remove_MtRb(seurat_objt = young2)
            young3 = Remove_MtRb(seurat_objt = young3)
            fetal1 = Remove_MtRb(seurat_objt = fetal1)
            fetal2 = Remove_MtRb(seurat_objt = fetal2)
            fetal3 = Remove_MtRb(seurat_objt = fetal3)
            
            
            # 2.2.2 .........................................................................................................................................
            # convert MACS2 narrowPeaks to Granges object
            macs_path="/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/Peak_macs2/"
            name <- c("efe001", "efe002","efe003","ctrl001","ctrl002","ctrl003")
            peaks_macs2 <- list()
            
            peaks_macs2 = Peak_Macs2(macs.path = macs_path, macs.name = name, peaks.macs2 = peaks_macs2)
            
            # Create a unified set of peaks to quantify in each dataset
            combined.peaks <- GenomicRanges::reduce(x = c(peaks_macs2[["efe001"]], peaks_macs2[["efe002"]], peaks_macs2[["efe003"]],peaks_macs2[["ctrl001"]], peaks_macs2[["ctrl002"]], peaks_macs2[["ctrl003"]]))
            
            # Filter out bad peaks based on length
            peakwidths <- width(combined.peaks)
            combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
            combined.peaks
            
            # get individual samples
            
            # IMPORTANT: take EFE001 and EFE002 samples to include peaks assay and then integrate
            
            # create a new assay using the MACS2 peak set and add it to the Seurat object
            
            EFE001 = Create_Peaks(seurat_objt=EFE001, fragpath = fragpath_efe1, comb.peaks = combined.peaks, annt=annotation)
            EFE002 = Create_Peaks(seurat_objt=EFE002, fragpath = fragpath_efe2, comb.peaks = combined.peaks, annt=annotation)
            EFE003 = Create_Peaks(seurat_objt=EFE003, fragpath = fragpath_efe3, comb.peaks = combined.peaks, annt=annotation)
            
            Ctrl001 = Create_Peaks(seurat_objt=Ctrl001, fragpath = fragpath_puctrl1, comb.peaks = combined.peaks, annt=annotation)
            Ctrl002 = Create_Peaks(seurat_objt=Ctrl002, fragpath = fragpath_puctrl2, comb.peaks = combined.peaks, annt=annotation)
            Ctrl003 = Create_Peaks(seurat_objt=Ctrl003, fragpath = fragpath_puctrl3, comb.peaks = combined.peaks, annt=annotation)
            
            # gene activities
            DefaultAssay(EFE001) <- DefaultAssay(EFE002) <- DefaultAssay(EFE003) <- DefaultAssay(Ctrl001) <- DefaultAssay(Ctrl002) <- DefaultAssay(Ctrl003) <- "peaks"
            efe1_gene.activities <- GeneActivity(EFE001)
            efe2_gene.activities <- GeneActivity(EFE002)
            efe3_gene.activities <- GeneActivity(EFE003)
            ctrl1_gene.activities <- GeneActivity(Ctrl001)
            ctrl2_gene.activities <- GeneActivity(Ctrl002)
            ctrl3_gene.activities <- GeneActivity(Ctrl003)
            
            # add gene activities as a new assay - unnormalized
            EFE001[["ACTIVITY"]] <- CreateAssayObject(counts = efe1_gene.activities)
            EFE002[["ACTIVITY"]] <- CreateAssayObject(counts = efe2_gene.activities)
            EFE003[["ACTIVITY"]] <- CreateAssayObject(counts = efe3_gene.activities)
            Ctrl001[["ACTIVITY"]] <- CreateAssayObject(counts = ctrl1_gene.activities)
            Ctrl002[["ACTIVITY"]] <- CreateAssayObject(counts = ctrl2_gene.activities)
            Ctrl003[["ACTIVITY"]] <- CreateAssayObject(counts = ctrl3_gene.activities)
            
            # 2.3 .........................................................................................................................................
            quantile(EFE001@meta.data$nCount_RNA,c(0.35,0.95))
            quantile(EFE002@meta.data$nCount_RNA,c(0.05,0.95))
            quantile(EFE003@meta.data$nCount_RNA,c(0.05,0.95))
            quantile(Ctrl001@meta.data$nCount_RNA,c(0.10,0.95))
            quantile(Ctrl002@meta.data$nCount_RNA,c(0.20,0.95))
            quantile(Ctrl003@meta.data$nCount_RNA,c(0.15,0.95))
            quantile(young1@meta.data$nCount_RNA,c(0.05,0.95))
            quantile(young2@meta.data$nCount_RNA,c(0.05,0.95))
            quantile(young3@meta.data$nCount_RNA,c(0.05,0.95))
            quantile(fetal1@meta.data$nCount_RNA,c(0.05,0.95))
            quantile(fetal2@meta.data$nCount_RNA,c(0.05,0.95))
            quantile(fetal3@meta.data$nCount_RNA,c(0.05,0.95))
            
            quantile(EFE001@meta.data$nCount_ATAC,c(0.02,0.98))
            quantile(EFE002@meta.data$nCount_ATAC,c(0.02,0.98))
            quantile(EFE003@meta.data$nCount_ATAC,c(0.02,0.98))
            quantile(Ctrl001@meta.data$nCount_ATAC,c(0.02,0.98))
            quantile(Ctrl002@meta.data$nCount_ATAC,c(0.02,0.98))
            quantile(Ctrl003@meta.data$nCount_ATAC,c(0.02,0.98))
            
            #### subset
            EFE001 <- subset(
              x = EFE001,
              subset = nCount_RNA < 7853 &
                nCount_RNA > 527 &
                nFeature_RNA > 300 &
                nCount_ATAC < 6094 &
                nCount_ATAC > 89 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 2 &
                TSS.enrichment > 1
            )
            
            EFE002 <- subset(
              x = EFE002,
              subset = nCount_RNA < 11374 &
                nCount_RNA > 582 &
                nFeature_RNA > 300 &
                nCount_ATAC < 8813 &
                nCount_ATAC > 114 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 2 &
                TSS.enrichment > 1
            )
            
            EFE003 <- subset(
              x = EFE003,
              subset = nCount_RNA < 5567 &
                nCount_RNA > 624 &
                nFeature_RNA > 300 &
                nCount_ATAC < 2456 &
                nCount_ATAC > 153 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 2 &
                TSS.enrichment > 1
            )
            
            Ctrl001 <- subset(
              x = Ctrl001,
              subset = nCount_RNA < 6535 &
                nCount_RNA > 590 &
                nFeature_RNA > 300 &
                nCount_ATAC < 19215 &
                nCount_ATAC > 250 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 4 &
                TSS.enrichment > 1
            )
            
            Ctrl002 <- subset(
              x = Ctrl002,
              subset = nCount_RNA < 4691 &
                nCount_RNA > 647 &
                nFeature_RNA > 300 &
                nCount_ATAC < 10919 &
                nCount_ATAC > 104 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 4 &
                TSS.enrichment > 1
            )
            
            Ctrl003 <- subset(
              x = Ctrl003,
              subset = nCount_RNA < 4786 &
                nCount_RNA > 516 &
                nFeature_RNA > 300 &
                nCount_ATAC < 12764 &
                nCount_ATAC > 102 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 4 &
                TSS.enrichment > 1
            )
            
            young1 <- subset(
              x = young1,
              subset = nCount_RNA < 17510 &
                nCount_RNA > 3652 &
                nFeature_RNA > 300 &
                percent.mt < 20 &
                percent.ribo < 10
            )
            
            young2 <- subset(
              x = young2,
              subset = nCount_RNA < 28954 &
                nCount_RNA > 1093 &
                nFeature_RNA > 300 &
                percent.mt < 20 &
                percent.ribo < 10
            )
            
            young3 <- subset(
              x = young3,
              subset = nCount_RNA < 18868 &
                nCount_RNA > 605 &
                nFeature_RNA > 300 &
                percent.mt < 20 &
                percent.ribo < 10
            )
            
            fetal1 <- subset(
              x = fetal1,
              subset = nCount_RNA < 19233 &
                nCount_RNA > 3576 &
                nFeature_RNA > 300 &
                percent.mt < 20 &
                percent.ribo < 10
            )
            
            fetal2 <- subset(
              x = fetal2,
              subset = nCount_RNA < 18958 &
                nCount_RNA > 2628 &
                nFeature_RNA > 300 &
                percent.mt < 20 &
                percent.ribo < 10
            )
            
            fetal3 <- subset(
              x = fetal3,
              subset = nCount_RNA < 25538 &
                nCount_RNA > 583 &
                nFeature_RNA > 300 &
                percent.mt < 20 &
                percent.ribo < 10
            )
            
            #saveRDS(EFE001,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/EFE001_beforedb_raw_121923.rds")
            #saveRDS(EFE002,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/EFE002_beforedb_raw_121923.rds")
            #saveRDS(EFE003,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/EFE003_beforedb_raw_121923.rds")
            #saveRDS(Ctrl001,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Ctrl001_beforedb_raw_121923.rds")
            #saveRDS(Ctrl002,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Ctrl002_beforedb_raw_121923.rds")
            #saveRDS(Ctrl003,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Ctrl003_beforedb_raw_121923.rds")
            #saveRDS(young1,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/young1_beforedb_raw_121923.rds")
            #saveRDS(young2,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/young2_beforedb_raw_121923.rds")
            #saveRDS(young3,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/young3_beforedb_raw_121923.rds")
            #saveRDS(fetal1,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/fetal1_beforedb_raw_121923.rds")
            #saveRDS(fetal2,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/fetal2_beforedb_raw_121923.rds")
            #saveRDS(fetal3,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/fetal3_beforedb_raw_121923.rds")
            
            # go back to doublet cells
            library(DoubletFinder)
            
            EFE001.sg = Remove_Doublet(seurat_objt = EFE001,pseu_rate = 0.040)
            EFE002.sg = Remove_Doublet(seurat_objt = EFE002,pseu_rate = 0.075)
            EFE003.sg = Remove_Doublet(seurat_objt = EFE003,pseu_rate = 0.100)
            Ctrl001.sg = Remove_Doublet(seurat_objt = Ctrl001,pseu_rate = 0.040)
            Ctrl002.sg = Remove_Doublet(seurat_objt = Ctrl002,pseu_rate = 0.023)
            Ctrl003.sg = Remove_Doublet(seurat_objt = Ctrl003,pseu_rate = 0.016)
            young1.sg = Remove_Doublet(seurat_objt = young1,pseu_rate = 0.031)
            young2.sg = Remove_Doublet(seurat_objt = young2,pseu_rate = 0.040)
            young3.sg = Remove_Doublet(seurat_objt = young3,pseu_rate = 0.046)
            fetal1.sg = Remove_Doublet(seurat_objt = fetal1,pseu_rate = 0.054)
            fetal2.sg = Remove_Doublet(seurat_objt = fetal2,pseu_rate = 0.075)
            fetal3.sg = Remove_Doublet(seurat_objt = fetal3,pseu_rate = 0.062)
            
            
            
            #### 3. Integrate samples------------------------------------------------------------------------------------------------------------------------------------------------------------
            # normalize and identify variable features for each dataset independently
            EFE1<-EFE001.sg$objt
            EFE2<-EFE002.sg$objt
            EFE3<-EFE003.sg$objt
            CTRL1<-Ctrl001.sg$objt
            CTRL2<-Ctrl002.sg$objt
            CTRL3<-Ctrl003.sg$objt
            Y1<-young1.sg$objt
            Y2<-young2.sg$objt
            Y3<-young3.sg$objt
            F1<-fetal1.sg$objt
            F2<-fetal2.sg$objt
            F3<-fetal3.sg$objt
            
            #saveRDS(EFE1,file="/Users/yangyu/Desktop/EFE1_adb_raw_121923.rds")
            #saveRDS(EFE2,file="/Users/yangyu/Desktop/EFE2_adb_raw_121923.rds")
            #saveRDS(EFE3,file="/Users/yangyu/Desktop/EFE3_adb_raw_121923.rds")
            #saveRDS(CTRL1,file="/Users/yangyu/Desktop/Ctrl1_adb_raw_121923.rds")
            #saveRDS(CTRL2,file="/Users/yangyu/Desktop/Ctrl2_adb_raw_121923.rds")
            #saveRDS(CTRL3,file="/Users/yangyu/Desktop/Ctrl3_adb_raw_121923.rds")
            #saveRDS(Y1,file="/Users/yangyu/Desktop/Y1_adb_raw_121923.rds")
            #saveRDS(Y2,file="/Users/yangyu/Desktop/Y2_adb_raw_121923.rds")
            #saveRDS(Y3,file="/Users/yangyu/Desktop/Y3_adb_raw_121923.rds")
            #saveRDS(F1,file="/Users/yangyu/Desktop/F1_adb_raw_121923.rds")
            #saveRDS(F2,file="/Users/yangyu/Desktop/F2_adb_raw_121923.rds")
            #saveRDS(F3,file="/Users/yangyu/Desktop/F3_adb_raw_121923.rds")
            
            # This ordering is essential
            all.list <- list(EFE1, EFE2, EFE3, CTRL1, CTRL2, CTRL3, Y1, Y2, Y3, F1, F2, F3)
            all.list[[1]][["STIM"]] <- "EFE1"
            all.list[[2]][["STIM"]] <- "EFE2"
            all.list[[3]][["STIM"]] <- "EFE3"
            all.list[[4]][["STIM"]] <- "CTRL1"
            all.list[[5]][["STIM"]] <- "CTRL2"
            all.list[[6]][["STIM"]] <- "CTRL3"
            all.list[[7]][["STIM"]] <- "young1"
            all.list[[8]][["STIM"]] <- "young2"
            all.list[[9]][["STIM"]] <- "young3"
            all.list[[10]][["STIM"]] <- "fetal1"
            all.list[[11]][["STIM"]] <- "fetal2"
            all.list[[12]][["STIM"]] <- "fetal3"
            
            # same as lapply because we did save data in doublet finder process
            #for(i in seq_along(all.list)) {
            #DefaultAssay(all.list[[i]]) <- "RNA"
            #all.list[[i]] <- NormalizeData(all.list[[i]])
            #all.list[[i]] <- FindVariableFeatures(
            #all.list[[i]],
            #selection.method = "vst",
            #nfeatures = 2000
            #)
            #}
            #
            
            for(i in 1:length(all.list)){
              DefaultAssay(all.list[[i]]) <- "RNA"
            }
            
            all.list <- lapply(X = all.list, FUN = function(x) {
              x <- NormalizeData(x)
              x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
            })
            
            # select features that are repeatedly variable across datasets for integration
            features <- SelectIntegrationFeatures(object.list = all.list)
            
            anchors <- FindIntegrationAnchors(object.list = all.list, anchor.features = features)
            
            
            
            # Save anchors done on 032725 which is already deposited to google drive rna_only_inte_12samp/ folder storage
            #saveRDS(anchors, file = "/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_anchors_rm_scaledata032725.rds")
            
            # this command creates an 'integrated' data assay
            
            # 
            #anchors <- FindIntegrationAnchors(
            #object.list = all.list,
            #anchor.features = features,
            #normalization.method = "LogNormalize"  # Force RNA assay's data slot
            #scale = TRUE  # Internal scaling during CCA (default; unrelated to scale.data)
            #)
            #
            
            combined.rna <- IntegrateData(anchorset = anchors)
            
            # specify that we will perform downstream analysis on the corrected data note that the
            # original unmodified data still resides in the 'RNA' assay
            # Perform an intergrated analysis
            DefaultAssay(combined.rna) <- "integrated"
            
            # Run the standard workflow for visualization and clustering
            combined.rna <- ScaleData(combined.rna, verbose = FALSE)
            combined.rna <- RunPCA(combined.rna, npcs = 50, verbose = FALSE) # ElbowPlot(combined.rna,ndims=50)
            combined.rna <- RunUMAP(combined.rna, reduction = "pca", dims = 1:30)
            combined.rna <- FindNeighbors(combined.rna, reduction = "pca", dims = 1:30)
            combined.rna <- FindClusters(combined.rna, resolution = 0.5)
            DimPlot(combined.rna, reduction = "umap", split.by = "STIM", pt.size = 0.01,label = T) + NoLegend()
            # save saveRDS(combined.rna,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_combined_raw_121923.rds")
            DefaultAssay(combined.rna) <- "RNA"
            RNA_all.markers <- FindAllMarkers(combined.rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
            inte.mkrs<-data.table(RNA_all.markers)
            # write.table(inte.mkrs, file = "/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_markers_by25clusters_121923.txt",quote=F,row.names=F,col.names=T,sep="\t")
            
            # Annotation on RNA data
            # DISCO Annotate ----------------------------------------------------------------------------------------------------------------------
            young_se1<-combined.rna
            rna.data.average = AverageExpression(young_se1)
            
            # This will generate average expression for each cluster
            rna.data.average = round(rna.data.average$RNA, 2)
            #write.table(rna.data.average, file="/Users/yangyu/Downloads/CELLiD_input.txt", quote = F, col.names = F, row.names = T, sep="\t")
            
            # Then, you can upload this file to predict cell types DISCO
            
            # After you get the results, you can add the predicted cell type to seurat object as follows:
            predicted.ct = read.csv("/Users/yangyu/Downloads/spreadsheet.csv")
            young_se1$primary.predict = predicted.ct[as.numeric(young_se1$seurat_clusters),2]
            young_se1$secondary.predict = predicted.ct[as.numeric(young_se1$seurat_clusters),3]
            DimPlot(young_se1, group.by = "primary.predict", label = T)+NoLegend()
            DimPlot(young_se1, group.by = "secondary.predict", label = T)
            
            #### 4. Annotate ----------------------------------------------------------------------------------------------------------------------
            library(scCATCH)
            
            T_cell=c("CD247","THEMIS") #~ #https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-023-04224-1/figures/1
            # T_cell=c("CD247","THEMIS","CD3D","CD3E","CD8A","CD4","TRBC2","CD3G","LTB","IL7R","LEF1","GZMK","TRAC","CXCR6","CD69","CCL5")
            #FeaturePlot(combined,reduction = "umap",features = T_cell,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =T_cell,pt.size = 0,ncol = 2)
            
            #NK=c("FGFBP2","KLRD1","GNLY","DOCK2") #~ #https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-023-04224-1/figures/1
            #NK=c("FGFBP2","KLRD1","GNLY","KLRF1","NKG7","DOCK2","GZMA","GZMB","NCR1","TGFB1")
            #FeaturePlot(combined,reduction = "umap",features = NK,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =NK,pt.size = 0,ncol = 2)
            
            #NK_T=c("NCAM1","IL2RB","CD44","IL12RB2","CXCR4") #~ #https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-023-04224-1/figures/1
            #FeaturePlot(combined,reduction = "umap",features = NK_T,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =NK_T,pt.size = 0,ncol = 2)
            
            #PKHD1L1 NRG3 CDH11 LINGO1
            
            EC=c("PECAM1","VWF","CDH5","BTNL9") #~ #Cdh5, Pecam1 EC=c("PECAM1","VWF","SMOC1","CD36","MEOX2","TCF15","FABP4") # Meox2/Tcf15, Fabp4, and Cd36
            #EC_1=c("PECAM1","VWF","CDH5","CD36","CD31","CD34","CD93","EGFL7","ID3","FLT1","GNG11","MCAM","FLT4","PLVAP","ADGRF5","ABCG2","MECOM","ETV2")
            #FeaturePlot(combined,reduction = "umap",features = EC,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =EC,pt.size = 0,ncol = 2)
            
            EC_2=c("PECAM1","EGFL7","MECOM") #~ #Cdh5, Pecam1 EC=c("PECAM1","VWF","SMOC1","CD36","MEOX2","TCF15","FABP4") # Meox2/Tcf15, Fabp4, and Cd36
            #EC_2=c("PECAM1","VWF","CDH5","CD36","CD31","CD34","CD93","EGFL7","ID3","FLT1","GNG11","MCAM","FLT4","PLVAP","ADGRF5","ABCG2","MECOM","ETV2")
            #FeaturePlot(combined,reduction = "umap",features = EC_2,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =EC_2,pt.size = 0,ncol = 2)
            
            
            
            EC_Card=c("LINGO1","BCOR","MYH11","MYL9")
            #FeaturePlot(combined,reduction = "umap",features = EC_Card,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =EC_Card,pt.size = 0,ncol = 2)
            
            #Fibroblast_Card=c("LINGO1","BCOR","MYH11","MYL9")
            #FeaturePlot(combined,reduction = "umap",features = Fibroblast_Card,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =Fibroblast_Card,pt.size = 0,ncol = 2)
            
            
            EndoC=c("SMOC1","NPR3") # Npr3
            #EndoC=c("SMOC1","NPR3","VWF")
            #FeaturePlot(combined,reduction = "umap",features = EndoC,split.by = "orig.ident")
            #VlnPlot(combined,features =EndoC,pt.size = 0,ncol = 2)
            
            Fibroblast=c("DCN","COL1A1","LUM") #~ #https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-023-04224-1/figures/1 
            #Fibroblast=c("DCN","ABCA6","COL4A4","LUM","COL1A1","COL6A2", "VTN", "MFAP5", "VIM", "PDGFRB", "POSTN", "ASPN", "PRRX1", "ENG",  "COL1A2","MKI67","PTX3")
            #Fibroblast2=c("NFATC1","FOSB","FABP4","FLI1","CXCL14","INMT","CKAP4","P4HTM","FGF2","VEGFA","NGF","GSTM5","MDK","TBX20","DKK3","LAMB1","MEDAG","LAMC1","ZEB2","TCF4","STAT3",
            #"RUNX1","PBX1","NR4A2","NR4A1","NFKB1","NFAT5","KLF9","DCN","COCH","TNN","DCD","MUCL1","APOD","PTGDS","BNC2","SCGB2A2","MUCL1","MKI67")
            #FeaturePlot(combined,reduction = "umap",features = Fibroblast,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features = Fibroblast,pt.size = 0,ncol = 2)
            
            cardiomyocyte=c("ACTN2", "TTN","FHL2") # Actn2, Tnni3, Ttn #~
            #cardiomyocyte=c("ACTN2","TTN","TECRL","FHL2","TNNI3", "MYH7","FGF2", "NPPB", "TNNT2", "ANKRD1", "NPPC", "GATA4", "DES", "HAND2", "GATA6", "MYL2", "CD36", "MYL3","NOVA1","CR2","ANKRD29","MKI67","RYR2","TNNI1")
            #cardiomyocyte2=c("ACTN2","TTN","NOVA1","CR2","ANKRD29","MKI67","DES","TNNI1","TNNT2","RYR2")
            #FeaturePlot(combined,reduction = "umap",features = cardiomyocyte,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features = cardiomyocyte,pt.size = 0,ncol = 2)
            
            M2_macrophage=c("CD163","MRC1") # ~
            #M2_macrophage=c("CD163","CD163L1","MRC1","RBM47","CD68","NAAA","JAML","TYROBP","CXCL16","SCIMP","PARP14","CD247","MARCH1")
            #FeaturePlot(combined,reduction = "umap",features = M2_macrophage,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features = M2_macrophage,pt.size = 0,ncol = 2)
            
            #EpiC=c("PDGFRB","ABCC9","NOTCH3") #~ # Wt1, Tbx18
            #EpiC=c("PDGFRB","ABCC9","NOTCH3","WT1", "TBX18", "SEMA3D", "ALDH1A2", "GATA5", "TCF21")
            #FeaturePlot(combined,reduction = "umap",features = EpiC,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =EpiC,pt.size = 0,ncol = 2)
            
            Pericyte=c("RGS5","KCNJ8","PDGFRB") #~ #https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-023-04224-1/figures/1
            #Pericyte=c("EGFLAM","RGS5","GUCY1A2","KCNJ8","LHFP","ABCC9","PDGFRB","ACTA2","TBX18","MCAM","HIGD1B","ZIC1","ANGPT1","CSPG4","DES")
            #FeaturePlot(combined,reduction = "umap",features = Pericyte,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =Pericyte,pt.size = 0,ncol = 2)
            
            smoothe_muscle=c("ACTA2","MYH11") #~  ##RGS5 is a marker of smoothe_muscle cells. #https://blog.csdn.net/zengwanqin/article/details/115964144
            #smoothe_muscle=c("ACTA2","MYH11","TAGLN","MYL9","RGS5","MYLK","NOV","SH3BGR","NOX4","KCNMB1","ITGA9")
            #FeaturePlot(combined,reduction = "umap",features = smoothe_muscle,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =smoothe_muscle,pt.size = 0,ncol = 2)
            
            Neuro=c("XKR4","PTPRZ1","KCNH8") #~ ##RGS5 is a marker of smoothe_muscle cells. #https://blog.csdn.net/zengwanqin/article/details/115964144
            #FeaturePlot(combined,reduction = "umap",features = Neuro,split.by = "orig.ident",pt.size = 0.0001, min.cutoff = "q9")
            #VlnPlot(combined,features =Neuro,pt.size = 0,ncol = 2)
            
            Erythroid_like=c("ALAS2","CD55")
            #FeaturePlot(combined,reduction = "umap",features = Inte,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =Inte,pt.size = 0,ncol = 2)
            
            Mast_cell=c("TPSB2","TPSAB1","KIT")
            #FeaturePlot(combined,reduction = "umap",features = Mast_cell,split.by = "orig.ident",pt.size = 0.0001)
            #VlnPlot(combined,features =Mast_cell,pt.size = 0,ncol = 2)
            
            # 4.3 adding cluster name ..............................................
            # 4 sample integrated 
            young_se1=combined.rna
            
            C0="Cardiomyocyte" # 
            C1="Cardiomyocyte"  #         
            C2="EC"                       
            C3="Fibroblast"  #!   
            C4="Fibroblast"  #                
            C5="Pericyte"  ###                     
            C6="Fibroblast"  #                    
            C7="EC"  #                      
            C8="EC-Card"  #            
            C9="Macrophage"    #      
            C10="Cardiomyocyte" #                      
            C11="SMC"  #!             
            C12="Cardiomyocyte" #      {PVT1; MTHFD1L}      dt.mg_mkrs[cluster==12][gene%in% dt.mg_mkrs[,list(n=.N),list(gene)][n==1]$gene]          
            C13="EC"  #!  Pericyte                  
            C14="EC"  #       
            C15="Cardiomyocyte"  # Cardiomyocyte       
            C16="EndoC" # EC                    
            C17="Macrophage"   #                
            C18="Fibroblast"   #  
            C19="Cardiomyocyte"  #                 
            C20="Fibroblast" # Fibroblast
            C21="Neuro"  # 
            C22="T-cell"
            C23="Mast-cell"
            C24="Macrophage" # 
            
            
            # 4.5 Annotate .........................................................
            
            young_se1 <- SetIdent(young_se1, value = "seurat_clusters")
            new.cluster.ids=c(C0,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22,C23,C24)
            names(new.cluster.ids) <- levels(young_se1)
            young_se1 <- RenameIdents(young_se1, new.cluster.ids)
            fact=c("EC","EndoC","Fibroblast","SMC","Cardiomyocyte","Macrophage", "T-cell","Pericyte","Neuro","Mast-cell")
            fact=c("EC","EC-Card","EndoC","Fibroblast","SMC","Cardiomyocyte","Macrophage", "T-cell","Pericyte","Neuro","Mast-cell")
            
            young_se1@active.ident=factor(young_se1@active.ident,levels=fact)
            young_se1@active.ident
            
            library(scRNAtoolVis)
            library(RColorBrewer)
            cell_type_cols=c("#FFA500","#00BFFF","#DDA0DD","#FF69B4","#6A5ACD","#9932CC",
                             "#DC143C","#8FBC8F","#87CEFA","#DB7093","#A0522D","#C71585",
                             "#32CD32","#00CED1","#8B0000","#BC8F8F","#F08080","#DAA520",
                             "#90EE90","#3CB371","#191970")
            
            young_se1$cell_type=young_se1@active.ident
            
            clusterCornerAxes(object = young_se1,reduction = 'umap',clusterCol = "cell_type",pSize=0.1,cellLabel=T,cellLabelSize=3,relLength = 0.3,arrowType="closed")+scale_color_manual(values = cell_type_cols) 
            #clusterCornerAxes(object = young_se1,reduction = 'umap',clusterCol = "cell_type",pSize=0.1,cellLabel=T,cellLabelSize=3,relLength = 0.3,arrowType="closed")+scale_color_manual(values = cell_type_cols) 
            
            
            # cluster Marker gene plot
            EC=EC
            #EC_Card=EC_Card
            #Fibroblast_Card=Fibroblast_Card
            EndoC=EndoC
            Fibroblast=Fibroblast
            SMC=smoothe_muscle
            Cardiomyocyte=cardiomyocyte
            Macrophage=M2_macrophage
            T_cell=T_cell
            Pericyte=Pericyte
            Neuro=Neuro
            #EpiC=EpiC
            #Erythroid_like=Erythroid_like
            Mast_cell=Mast_cell
            
            # fact=c("EC","EndoC","Fibroblast","Fibroblast-Card","SMC","Cardiomyocyte","Macrophage", "T-cell","Pericyte","Neuro")
            marker_gene=c(EC, EndoC,Fibroblast,SMC,Cardiomyocyte,Macrophage,T_cell,Pericyte,Neuro,Mast_cell)
            dat=list(EC,EndoC,Fibroblast,SMC,Cardiomyocyte,Macrophage,T_cell,Pericyte,Neuro,Mast_cell)
            cluster_gene=NULL
            for (i in 1:length(dat)){
              a=rep(i,length(dat[[i]]))
              cluster_gene=c(cluster_gene,a)
            }
            
            dat1=data.frame(marker_gene,cluster_gene)
            dat1=dat1[!duplicated(dat1$marker_gene),]
            dat1=dat1[order(dat1$cluster_gene,decreasing = T),]
            
            # additional steps
            #coarse_order = c("EC","EndoC","Fibroblast","SMC","Cardiomyocyte","Macrophage", "T-cell","Pericyte","Neuro","Mast-cell")
            #young_se1$coarse_type <- factor(young_se1$coarse_type, levels = coarse_order)
            
            p1=jjDotPlot(object = young_se1,gene = dat1$marker_gene,id = 'coarse_type',rescale = T,ytree = F)
            p1
            
            young_se1 <- SetIdent(young_se1, value = "cell_type")
            library(dplyr)
            integrated.markers <- FindAllMarkers(young_se1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
            
            dt_intmkr<-data.table(integrated.markers)
            
            #### save 
            # write.table(dt_intmkr,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_marker_genes_121923.txt",quote=F,row.names=F,col.names=T,sep="\t")
            # saveRDS(young_se1,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_combined_annotated_121923.rds")
            
            # group samples
            young_se1[["SAMPLE"]] <- ""
            young_se1@meta.data$SAMPLE[young_se1@meta.data$STIM %in% c("EFE1","EFE2","EFE3")] <- "EFEs"
            young_se1@meta.data$SAMPLE[young_se1@meta.data$STIM %in% c("CTRL1","CTRL2","CTRL3")] <- "CTRLs"
            young_se1@meta.data$SAMPLE[young_se1@meta.data$STIM %in% c("fetal1","fetal2","fetal3")] <- "Fetals"
            young_se1@meta.data$SAMPLE[young_se1@meta.data$STIM %in% c("young1","young2","young3")] <- "Youngs"
            
            Samp_order = c("EFEs","CTRLs","Fetals","Youngs")
            young_se1@meta.data$SAMPLE <- factor(young_se1@meta.data$SAMPLE, levels = Samp_order)
            
            #clusterCornerAxes(object = young_se1,split.by = "SAMPLE", reduction = 'umap',clusterCol = "cell_type",pSize=0.1,cellLabel=T,cellLabelSize=3,relLength = 0.3,arrowType="closed")+scale_color_manual(values = cell_type_cols)
            DimPlot(young_se1, split.by = "SAMPLE",cols = cell_type_cols[c(1,4,5,6,7,8,9,10,11,12,13,14)], label = T) + NoLegend()
            table(Idents(young_se1), young_se1@meta.data$STIM)
            
            # compare cell percentage in each cluster
            pt <- table(Idents(young_se1), young_se1$STIM)
            pt <- as.data.frame(pt)
            pt$Var1 <- as.character(pt$Var1)
            
            pt$Var1 <- factor(pt$Var1,levels = fact)
            ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
              theme_bw(base_size = 15) +
              geom_col(position = "fill", width = 0.5) +
              xlab("Sample") +
              ylab("Proportion") +
              scale_fill_manual(values = cell_type_cols[c(1,4,5,6,7,8,9,10,11,12,13,14)]) +
              theme(legend.title = element_blank())
            
            
            # name coarse cell type
            young_se1[["coarse_type"]] <- ""
            young_se1@meta.data$coarse_type[young_se1@meta.data$cell_type %in% c("EC","EC-Card")] <- "EC"
            young_se1@meta.data$coarse_type[young_se1@meta.data$cell_type %in% c("Cardiomyocyte")] <- "CM" # Cardiomyocyte
            young_se1@meta.data$coarse_type[young_se1@meta.data$cell_type %in% c("Fibroblast")] <- "FB" # Fibroblast
            young_se1@meta.data$coarse_type[young_se1@meta.data$cell_type %in% c("SMC")] <- "SMC"
            young_se1@meta.data$coarse_type[young_se1@meta.data$cell_type %in% c("Pericyte")] <- "PeriC" # Pericyte
            young_se1@meta.data$coarse_type[young_se1@meta.data$cell_type %in% c("Neuro")] <- "Neuro"
            young_se1@meta.data$coarse_type[young_se1@meta.data$cell_type %in% c("Macrophage")] <- "Macrophage"
            young_se1@meta.data$coarse_type[young_se1@meta.data$cell_type %in% c("T-cell")] <- "T-cell"
            young_se1@meta.data$coarse_type[young_se1@meta.data$cell_type %in% c("Mast-cell")] <- "Mast-cell"
            young_se1@meta.data$coarse_type[young_se1@meta.data$cell_type %in% c("EndoC")] <- "EndoC"
            
            young_se1 <- SetIdent(young_se1, value = "coarse_type")
            coarse_order = c("EC","EndoC","Fibroblast","SMC","Cardiomyocyte","Macrophage", "T-cell","Pericyte","Neuro","Mast-cell")
            young_se1@active.ident=factor(young_se1@active.ident,levels=coarse_order)
            
            DimPlot(young_se1, split.by = "SAMPLE",cols = cell_type_cols[c(1,4,5,6,7,8,9,10,11,12)], label = T) + NoLegend()
            
            coarse.markers <- FindAllMarkers(young_se1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
            coarse_intmkr<-data.table(coarse.markers)
            
            #### save 
            # write.table(coarse_intmkr,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_coarsemarker_genes_011024.txt",quote=F,row.names=F,col.names=T,sep="\t")
            
            
            # compare cell percentage in each cluster
            pt <- table(Idents(young_se1), young_se1$STIM)
            pt <- as.data.frame(pt)
            pt$Var1 <- as.character(pt$Var1)
            
            pt$Var1 <- factor(pt$Var1,levels = coarse_order)
            pt$Var2 <- factor(pt$Var2,levels = samp_order)
            ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
              theme_bw(base_size = 15) +
              geom_col(position = "fill", width = 0.5) +
              xlab("Sample") +
              ylab("Proportion") +
              scale_fill_manual(values = cell_type_cols[c(1,4,5,6,7,8,9,10,11,12,13,14)]) +
              theme(legend.title = element_blank())
            
            
            
            # name specific cell type
            young_se1[["refined_type"]] <- ""
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("7")] <- "EFE3.EC"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("8")] <- "EFE1.EC"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("13")] <- "Art.EC"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("14")] <- "Ven.EC"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("2")] <- "Cap.EC"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("0","1","10","12","19")] <- "Vent.CM"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("15")] <- "CM"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("3","4","6","18","20")] <- "FB"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("11")] <- "SMC"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("5")] <- "Pericyte"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("21")] <- "Neuro"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("9","17","24")] <- "Macrophage"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("22")] <- "T-cell"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("23")] <- "Mast-cell"
            young_se1@meta.data$refined_type[young_se1@meta.data$seurat_clusters %in% c("16")] <- "EndoC"
            
            young_se1 <- SetIdent(young_se1, value = "refined_type")
            DimPlot(young_se1, split.by = "SAMPLE",cols = cell_type_cols, label = T) + NoLegend()
            
            # compare cell percentage in each cluster
            pt <- table(Idents(young_se1), young_se1$STIM)
            pt <- as.data.frame(pt)
            pt$Var1 <- as.character(pt$Var1)
            
            library(RColorBrewer)
            #pt$Var1 <- factor(pt$Var1,levels = fact)
            ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
              theme_bw(base_size = 15) +
              geom_col(position = "fill", width = 0.5) +
              xlab("Sample") +
              ylab("Proportion") +
              scale_fill_manual(values = cell_type_cols) +
              theme(legend.title = element_blank())
            
            #saveRDS(young_se1,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_combined_refinedannt_121923.rds")
            
            
            # load data
            young_se1 = readRDS("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_combined_refinedannt_121923.rds")
            
            young_se1 <- SetIdent(young_se1, value = "seurat_clusters")
            DimPlot(young_se1, split.by = "SAMPLE",label = T) + NoLegend()
            
            # name coarse cell type
            young_se1[["Cell_Type"]] <- ""
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("Cardiomyocyte")] <- "CM"
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("Fibroblast")] <- "FB"
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("SMC")] <- "SMC"
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("Pericyte")] <- "PeriC"
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("Neuro")] <- "Neuro"
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("Macrophage")] <- "Macrophage"
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("T-cell")] <- "T-cell"
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("Mast-cell")] <- "Mast-cell"
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("EndoC")] <- "EndoC"
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("EC")] <- "EC"
            
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("EC") & young_se1@meta.data$seurat_clusters %in% c("7")] <- "EC_3"
            young_se1@meta.data$Cell_Type[young_se1@meta.data$cell_type %in% c("EC-Card")] <- "EC_2"
            # young_se1@meta.data$Cell_Type[young_se1@meta.data$Cell_type %in% c("EC")] <- "EC_1"
            
            # re-order clusters
            young_se1 <- SetIdent(young_se1, value = "Cell_Type")
            cell_order = c("EC_1","EC_2","EC_3","EndoC","FB","SMC","CM","Macrophage", "T-cell","PeriC","Neuro","Mast-cell")
            samp_order = c("EFE1","EFE2","EFE3","CTRL1","CTRL2","CTRL3","fetal1","fetal2","fetal3","young1","young2","young3")
            young_se1@active.ident=factor(young_se1@active.ident,levels=cell_order)
            
            DimPlot(young_se1, split.by = "SAMPLE",cols = cell_type_cols, label = T) + NoLegend()
            
            # compare cell percentage in each cluster
            pt <- table(Idents(young_se1), young_se1$STIM)
            pt <- as.data.frame(pt)
            pt$Var1 <- as.character(pt$Var1)
            
            library(RColorBrewer)
            pt$Var1 <- factor(pt$Var1,levels = cell_order)
            pt$Var2 <- factor(pt$Var2,levels = samp_order)
            ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
              theme_bw(base_size = 15) +
              geom_col(position = "fill", width = 0.5) +
              xlab("Sample") +
              ylab("Proportion") +
              scale_fill_manual(values = cell_type_cols) +
              theme(legend.title = element_blank())
            
            
            #DimPlot(young_se1, cols = cell_type_cols, label = T) + NoLegend()
            DimPlot(young_se1, cols = cell_type_cols, label = T) + NoLegend()
            
            # Cell_Type markers
            young_se1 <- SetIdent(young_se1, value = "Cell_Type")
            
            library(dplyr)
            
            CT.markers <- FindAllMarkers(young_se1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
            ct_intmkr<-data.table(CT.markers)
            
            # write.table(ct_intmkr,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_grandmarker_genes_010424.txt",quote=F,row.names=F,col.names=T,sep="\t")
            # saveRDS(young_se1,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_combined_grandannt_010424.rds")
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
