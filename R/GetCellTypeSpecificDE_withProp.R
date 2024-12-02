#' Identify cross-regional cell type-specific differential analysis with proportion
#'
#' @param sce A single cell experiment object.
#' @param Regional1ID A vector of spot IDs for comparison region 1.
#' @param Regional2ID A vector of spot IDs for comparison region 2.
#' @param n_markers A number specifying the top DE gene number.
#' @param angle A number for angle when plotting.
#' @param hjust A number for horizontal justification when plotting.
#' @param size A number for text font size.
#' @param padj_filter A number for filtering adjusted p values.
#' @param doHeatmap Logical variable for whether drawing the heatmap.
#'
#' @importFrom SingleCellExperiment logcounts
#'
#' @return A list including the top DE genes (topDE), and all DE genes (allDE).
#' @export
#' @examples
#' data("example_sce")
#' example_sce <- mySpatialPreprocess(example_sce, platform="Visium")
#' thisID1 <- S4Vectors::metadata(example_sce)$selectCenters$selectID[1]
#' thisRadius1 <- S4Vectors::metadata(example_sce)$selectCenters$selectRadius[1]
#' OutRegRes1 <- RegionalST::FindRegionalCells(example_sce,
#'                                             centerID = thisID1,
#'                                             radius = thisRadius1,
#'                                             enhanced = FALSE,
#'                                             doPlot = FALSE,
#'                                             returnPlot = FALSE)
#' thisID2 <- S4Vectors::metadata(example_sce)$selectCenters$selectID[2]
#' thisRadius2 <- S4Vectors::metadata(example_sce)$selectCenters$selectRadius[2]
#' OutRegRes2 <- RegionalST::FindRegionalCells(example_sce,
#'                                             centerID = thisID2,
#'                                             radius = thisRadius2,
#'                                             enhanced = FALSE,
#'                                             doPlot = FALSE,
#'                                             returnPlot = FALSE)
#' Regional1ID <- OutRegRes1$closeID
#' Regional2ID <- OutRegRes2$closeID
#' CTS_DE <- GetCellTypeSpecificDE_withProp(example_sce,
#'                                          Regional1ID = Regional1ID,
#'                                          Regional2ID = Regional2ID,
#'                                          n_markers = 10,
#'                                          angle = 30,
#'                                          hjust = 0,
#'                                          size = 3,
#'                                          padj_filter = 0.05,
#'                                          doHeatmap = FALSE)
#'                                          
GetCellTypeSpecificDE_withProp <- function(sce,
                                           Regional1ID,
                                           Regional2ID,
                                           n_markers = 10,
                                           angle = 30,
                                           hjust = 0,
                                           size = 3,
                                           padj_filter = 0.05,
                                           doHeatmap = FALSE) {
    
    stopifnot(exprs = {
        is.numeric(n_markers)
        is.numeric(angle)
        is.numeric(hjust)
        is.numeric(size)
        is.null(padj_filter) | is.numeric(padj_filter)
    })
    
    tmp <- colnames(sce@metadata$Proportions)
    tmp1 <- gsub(" ", "_", tmp)
    propName <- gsub("-", "_", tmp1)
    
    newcoldata <- data.frame(group = c(rep("Reg1", length(Regional1ID)),
                                       rep("Reg2", length(Regional2ID))))
    newcoldata$group <- as.factor(newcoldata$group)
    
    Prop1 <- sce@metadata$Proportions[match(Regional1ID, sce$spot), ]
    Prop2 <- sce@metadata$Proportions[match(Regional2ID, sce$spot), ]
    colnames(Prop1) <- propName
    colnames(Prop2) <- propName
    AllProp <- rbind(Prop1, Prop2)
    message("Building design matrix...")
    Design_out <- TOAST::makeDesign(newcoldata, AllProp)
    
    allcloseID <- c(Regional1ID, Regional2ID)
    Y_raw <- as.matrix(logcounts(sce)[, match(allcloseID, colnames(SingleCellExperiment::logcounts(sce)))])
    message("Fitting model...")
    fitted_model <- TOAST::fitModel(Design_out, Y_raw)
    
    allDE <- c()
    allCT <- unique(colnames(Prop1))
    message("Test...")
    for(i in seq_len(length(allCT))) {
        message(paste0("Processing to cell type:", allCT[i]))
        if(i == 1) {
            newDEres <- TOAST::csTest(fitted_model,
                                      coef = "group",
                                      cell_type = allCT[i])
            newDEres$`Comparison` <- paste0(allCT[i], ": Region 1 vs Region 2")
            newDEres$`gene` <- rownames(newDEres)
            allDE <- newDEres
            
        } else {
            newDEres <- TOAST::csTest(fitted_model,
                                      coef = "group",
                                      cell_type = allCT[i])
            newDEres$`Comparison` <- paste0(allCT[i], ": Region 1 vs Region 2")
            newDEres$`gene` <- rownames(newDEres)
            allDE <- rbind(allDE, newDEres)
        }
    }
    
    if(is.null(padj_filter)) {
        top_markers <- allDE %>%
            dplyr::group_by(Comparison) %>%
            dplyr::top_n(n_markers, effect_size)
    } else {
        allDE <- allDE[allDE$fdr <= padj_filter,]
        top_markers <- allDE %>%
            dplyr::group_by(Comparison) %>%
            dplyr::top_n(n_markers, effect_size)
    }
    
    if (doHeatmap) {
        thiscoldata <- SummarizedExperiment::colData(sce)
        thiscoldata$TwoGroupInfo <- rep(NA, nrow(thiscoldata))
        
        thiscoldata$TwoGroupInfo[match(OneRegOut1$closeID, sce$spot)] <- twoCenter[1]
        thiscoldata$TwoGroupInfo[match(OneRegOut2$closeID, sce$spot)] <- twoCenter[2]
        thiscoldata$NewCellType <- thiscoldata$TwoGroupInfo
        
        seuratObj <- Seurat::CreateSeuratObject(counts=SingleCellExperiment::logcounts(sce)[, match(unique(allcloseID), colnames(SingleCellExperiment::logcounts(sce)))],
                                                assay='Spatial')
        tmpmetadata <- as.data.frame(thiscoldata)[match(unique(allcloseID), colnames(SingleCellExperiment::logcounts(sce))), ]
        rownames(tmpmetadata) <- colnames(seuratObj)
        seuratObj@meta.data <- tmpmetadata
        Seurat::Idents(seuratObj) <- seuratObj@meta.data$NewCellType
        
        seuratObj <- Seurat::NormalizeData(seuratObj)
        seuratObj <- Seurat::ScaleData(seuratObj)
        
        ## Scale data
        seuratObj@assays$Spatial$scale.data <-
            seuratObj@assays$Spatial$counts %>% as.matrix %>% t %>% scale %>% t
        
        palette <- RColorBrewer::brewer.pal(12, "Paired")
        Seurat::Idents(seuratObj) <- factor(Seurat::Idents(seuratObj), levels = c(twoCenter))
        
        if(length(top_markers$gene) == 0) {
            stop("There is no significant genes to plot! Try increase your filtering criteria.")
        }
        ## Plot expression of markers
        p1 <- Seurat::DoHeatmap(seuratObj, features = top_markers$gene, slot='counts',
                                group.colors=palette,
                                angle=angle, size=size, hjust = hjust, label = TRUE, raster=FALSE) +
            ggplot2::guides(col = FALSE)
        print(p1)
    }
    return(list(allDE = na.omit(allDE),
                topDE = na.omit(top_markers)))
}
