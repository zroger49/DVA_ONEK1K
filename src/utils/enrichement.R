suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(org.Hs.eg.db))

#Main analysis function
perform_enrichment <- function(genes_up, genes_down, background, FDR_enrichement = 0.05){
  res <- c()  
  
  deg_tissue_list_entrez.up <- c()
  deg_tissue_list_entrez.down <- c()
  
  out <- tryCatch({
    deg_tissue_list_entrez.up <- bitr(genes_up, fromType ="SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)["ENTREZID"]
  },error=function(cond) {
    deg_tissue_list_entrez.up <- c()
  })
  
  out <- tryCatch({
    deg_tissue_list_entrez.down <- bitr(genes_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)["ENTREZID"]
  },error=function(cond) {
    deg_tissue_list_entrez.down <- c()
  })
  

  background_entrez <- bitr(background, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)["ENTREZID"]
  
  # Enrichment analysis
  if (length(genes_up) > 5){
    go.bp.up <- enrichGO(gene = genes_up, keyType = "SYMBOL", OrgDb='org.Hs.eg.db', pvalueCutoff = FDR_enrichement, ont = "BP", pAdjustMethod = "BH", universe = background)
    go.mf.up <- enrichGO(gene = genes_up, keyType = "SYMBOL", OrgDb='org.Hs.eg.db', pvalueCutoff = FDR_enrichement, ont = "MF", pAdjustMethod = "BH", universe = background)
    go.cc.up <- enrichGO(gene = genes_up, keyType = "SYMBOL", OrgDb='org.Hs.eg.db', pvalueCutoff = FDR_enrichement, ont = "CC", pAdjustMethod = "BH", universe = background)
    res[["BP_up"]] <- go.bp.up
    res[["MF_up"]] <- go.mf.up
    res[["CC_up"]] <- go.cc.up
  }
  
  if (length(genes_down) > 5){
    go.bp.down <- enrichGO(gene = genes_down, keyType = "SYMBOL", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont = "BP", pAdjustMethod = "BH", universe = background)
    go.mf.down <- enrichGO(gene = genes_down, keyType = "SYMBOL", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont = "MF", pAdjustMethod = "BH", universe = background)
    go.cc.down <- enrichGO(gene = genes_down, keyType = "SYMBOL", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont = "CC", pAdjustMethod = "BH", universe = background)
    res[["BP_down"]] <- go.bp.down
    res[["MF_down"]] <- go.mf.down
    res[["CC_down"]] <- go.cc.down
  }
  
  if (length(deg_tissue_list_entrez.up$ENTREZID) > 5){
    kegg.up <- enrichKEGG(gene = deg_tissue_list_entrez.up$ENTREZID,  organism = "hsa", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
    do.up <- enrichDO(gene = deg_tissue_list_entrez.up$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
    dgn.up <- enrichDGN(gene = deg_tissue_list_entrez.up$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
    res[["kegg_up"]] <- kegg.up
    res[["do_up"]] <- do.up
    res[["dng_up"]] <- dgn.up
  }
  
  if (length(deg_tissue_list_entrez.down$ENTREZID) > 5){
    kegg.down <- enrichKEGG(gene = deg_tissue_list_entrez.down$ENTREZID,  organism = "hsa", pvalueCutoff = FDR_enrichement, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
    do.down <- enrichDO(gene = deg_tissue_list_entrez.down$ENTREZID, pvalueCutoff = FDR_enrichement, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
    dgn.down <- enrichDGN(gene = deg_tissue_list_entrez.down$ENTREZID, pvalueCutoff = FDR_enrichement, pAdjustMethod = "BH", universe = background_entrez$ENTREZID)
    
    res[["kegg_down"]] <- kegg.down
    res[["do_down"]] <- do.down
    res[["dng_down"]] <- dgn.down
  }
  res
}
