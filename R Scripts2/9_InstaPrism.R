library(BayesPrism)
library(InstaPrism)
library(Seurat)
library(tidyverse)


counts <- read.csv("Dataframes/Bulk Expression data/TCGA/tcga.paad.raw.counts.csv", row.names = 1, check.names = F)


bk.dat <- as.matrix(t(baseline.riaz.counts.anno))
dim(bk.dat)

sample.atlas <- readRDS("Consensus method/Figure 3/Dataframes/reference.atlas.rds")



Idents(sample.atlas) <- "GP"

levels(sample.atlas)

sample.atlas <- subset(sample.atlas, subset = GP != c("Mast"))

sample.atlas <- subset(sample.atlas, subset = GP != c("Plasma"))

sample.atlas <- subset(sample.atlas, subset = GP != c("Schwann"))


sc.dat <- GetAssayData(object = sample.atlas,slot = "counts")


dim(sc.dat)

sc.dat <- as.matrix(t(sc.dat))

anno.df <- sample.atlas@meta.data %>% 
  filter(str_detect(GP, "_")) %>% 
  rownames_to_column(var = "Cell.ID") %>% 
  select(GP, Cell.ID) %>% 
  separate(GP, sep = "_", into = c("Cell type", "GP2"), remove = F ) %>% 
  mutate(`Cell type` = case_when(GP2 %in% c("TFF1", "KRT17") ~ "Malignant",
                                 T ~ `Cell type`))

identical(anno.df$Cell.ID, rownames(sc.dat))


cell.type.labels <- anno.df$`Cell type`
cell.state.labels <- anno.df$GP

sort(table(cell.type.labels))
sort(table(cell.state.labels))

table(cbind.data.frame(cell.state.labels, cell.type.labels))

plot.cor.phi (input=sc.dat, 
              input.labels=cell.type.labels, 
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=1, cexCol=1,
)

plot.cor.phi (input=sc.dat, 
              input.labels=cell.state.labels, 
              title="cell state correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=1, cexCol=1,
)



sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.state.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)


sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

#note this function only works for human data. For other species, you are advised to make plots by yourself.
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)

sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")


myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key = "Malignant",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

InstaPrism.res.initial = InstaPrism(input_type = 'prism',
                                    prismObj = myPrism,n.core = 8)


bp.res <- InstaPrism.res.initial

save(bp.res, file="Consensus method/Figure 6/Assets/tcga.bp.rdata")


