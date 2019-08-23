# 加载 R 包
# install.packages("WGCNA")
library(WGCNA)
#=====================================================================================
#
#  1.导入表达数据
#
#=====================================================================================
{
# 不要将文件中的字符串转换为因子
  options(stringsAsFactors = FALSE)
# 过滤掉在所有样品中表达之和小于10 的基因
  # exp <- exp0[rowSums(exp0) > 10, ]
  
  # 转置
  datExpr <- exp[order(apply(exp,1,mad), decreasing = T)[1:10000],]
  #之前转换一下SYMBOL_id,并选择筛选30d表达
  datExpr = t(exp_30d_dup)
  dim(datExpr)
  rownames(datExpr)
  
  # 保存数据
  write.table(datExpr, file = "datExpr.txt",
              sep = "\t", quote = F, row.names = T)
}

#=====================================================================================
#
#  2.寻找最佳 β 值
#
#=====================================================================================
{
  # 开启多线程模式
  enableWGCNAThreads(nThreads = 20)
  
  # 通过对 power 的多次迭代，确定最佳 power
  powers <- c(1:30)
  sft = pickSoftThreshold(datExpr, 
                          powerVector = powers, 
                          verbose = 5,
                          networkType = "signed"
                          )
  # 计算出来的最佳 β 存放在
  sft$powerEstimate
  
  # 画图
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9
  #  R2 ~ soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  
  abline(h=0.8,col="red")
  # Mean connectivity ~ soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  sft$powerEstimate = 12
}

#=====================================================================================
#
#  3. 构建网络
##  3.1. 计算相关系数
##  3.2. 计算邻接矩阵
##  3.3. 计算 TOM 矩阵
##  3.4. 聚类并划分模块
##  3.5. 合并相似模块
#=====================================================================================
{
  net = blockwiseModules(
    # 0.输入数据
    datExpr, maxBlockSize=6338,
    
    # 1. 计算相关系数
    corType = "pearson", # 相关系数算法，pearson|bicor
    
    # 2. 计算邻接矩阵
    power = sft$powerEstimate, # 前面得到的 soft power
    networkType = "unsigned", # unsigned | signed | signed hybrid
    
    # 3. 计算 TOM 矩阵
    TOMType = "unsigned", # none | unsigned | signed
    saveTOMs = TRUE, # 是否保存
    saveTOMFileBase = "blockwiseTOM", # 保存文件前缀
    
    # 4. 聚类并划分模块
    deepSplit = 2, # 0|1|2|3|4, 值越大得到的模块就越多越小
    minModuleSize = 30,
    
    # 5. 合并相似模块
    ## 5.1 计算模块特征向量（module eigengenes， MEs），即第一个主成分（1st PC）
    ## 5.2 计算 MEs 与 datTrait 之间的相关性
    ## 5.3 对距离小于 mergeCutHeight （1-cor）的模块进行合并
    mergeCutHeight = 0.15, 

    # 其他参数
    numericLabels = TRUE, # 以数字命名模块
    nThreads = 0 # 0 则使用所有可用线程
    )
  # 查看每个模块包含基因数目
  table(net$colors) 
}


#=====================================================================================
#
#  4. 可视化
#
#=====================================================================================
{
  # 模块名称修改为颜色
  moduleColors = labels2colors(net$colors)
  # 同时绘制聚类图和模块颜色
  
  plotDendroAndColors(
    net$dendrograms[[1]], 
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, 
    addGuide = TRUE)
  dev.off()
}

save(datExpr, sft, net, moduleColors, file = "wgcna-network.Rdata")
#=====================================================================================
#
#  5. 导入性状数据
#
#=====================================================================================
{
  datTraits <- read.delim(f_trait, row.names=1, 
                          quote="")
  # 确保 datTraits 与 datExpr 的 rownames 顺序一致
  rownames(datTraits)
  traitRows = match(rownames(datExpr), rownames(datTraits))
  datTraits = datTraits[traitRows, ]
}

#=====================================================================================
#
#  6. 模块与性状关系
## 6.1 计算模块特征向量（moduleEigengenes， MEs）
## 6.2 计算模块（MEs）与性状之间的相关性
## 6.3 相关性 heatmap
#=====================================================================================
{
  # 6.0 获得样本数目和基因数目
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  # 6.1 计算模块特征向量 MEs
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  
  # 6.2 计算模块（MEs）与性状之间的相关性
  moduleTraitCor = cor(MEs, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

  # 6.3 相关性 heatmap
  sizeGrWindow(10,6)
  ## 连接相关性和 pvalue
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  
  
  ## heatmap 画图
  pdf(file = "Module-trait_relationships.pdf")
  par(mar = c(8, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  }

#=====================================================================================
#
#  6. 基因与性状关系（GS）& 基因与模块关系（MM）
## 6.1 计算 module membership (MM): 基因（TMP）与模块（MEs）的相关性
## 6.2 计算 Gene Significance (GS): 基因（TMP）与性状的相关性
#   
#=====================================================================================
{
  modNames = substring(names(MEs), 3)
  traitNames = names(datTraits)
  
  # 计算 MM
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                            nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  
  # 计算 GS
  geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                            nSamples))
  names(geneTraitSignificance) = paste("GS.", traitNames, sep="");
  names(GSPvalue) = paste("p.GS.", traitNames, sep="");

  geneInfo<-cbind(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)
  
  nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#=====================================================================================
#
#  7. 基因共表达网络可视化
#   
#=====================================================================================
{
  ## 7.1. 使用所有基因绘制
  # 计算 dissTOM, dissTOM = 1 - TOM
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate);
  
  # 取 7 次方，仅为展示更显著
  plotTOM = dissTOM^7
  
  # 绘图
  diag(plotTOM) = NA;
  sizeGrWindow(9,9)
  geneTree = net$dendrograms[[1]]
  TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  
  ## 7.2 随机选择 400 条绘图
  nSelect = 400
  set.seed(10)
  select = sample(nGenes, size = nSelect)
  selectTOM = dissTOM[select, select]
  # 重新聚类
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select]
  # Open a graphical window
  sizeGrWindow(9,9)
  plotDiss = selectTOM^7
  # 绘图
  diag(plotDiss) = NA;
  
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
}

#=====================================================================================
#
#  8. 模块特征向量网络可视化
#   
#=====================================================================================
{
  # 重新计算MEs
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5)
  par(cex = 0.9)
  plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                        = 90)
  
  # Plot the dendrogram
  sizeGrWindow(6,6)
  par(cex = 1.0)
  plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                        plotHeatmaps = FALSE)
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90)
}

#=====================================================================================
#
#  9. 生成 Cytoscape 的输入文件
#   
#=====================================================================================
{
  # 重新计算 TOM
  TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate);
  # 要可视化的模块
  modules = c("brown")
  # 要可视化的基因
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule]
  # 候选基因的 TOM
  modTOM = TOM[inModule, inModule]
  
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])
  
  
  
  
  
  modules = c("green")
   # 要可视化的基因
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule]
  # 候选基因的 TOM
  modTOM = TOM[inModule, inModule]
  
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])

