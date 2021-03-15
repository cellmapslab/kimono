###############################################
###############################################
###### Minimal working example - KiMONo #######
###############################################
###############################################


###############################################
### 1 libraries
###############################################
library(kimono)
library(data.table)


###############################################
### 2 Data
###############################################
# read in data
layer1 <- fread("data/expr.csv")
layer2 <- fread("data/pheno.csv")

# read in mapping (prior)
mapping11 <- fread("data/mapping_expr.csv")
mapping12 <- fread("data/mapping_expr_pheno.csv")

# make sure samples are in the same order
idorder <- layer1$sample
layer2 <- layer2[match(idorder, layer2$sample),]

# remove sample column
layer1 <- layer1[,-"sample"]
layer2 <- layer2[,-"sample"]


###############################################
# 3 Assemble into lists
###############################################
# data list
input_list <- list(
  as.data.table(layer1),
  as.data.table(layer2)
)
names(input_list) <- c('expr',
                       'pheno')
#########################
# mapping list
mapping_list <- list(
  as.data.table(mapping11),
  as.data.table(mapping12)
)
#########################
# meta info
metainfo <-data.frame('ID'   = c( 'prior_expr', 'expr_pheno'),
                      'main_to'   =  c(1,2)
)


###############################################
# 4 Run MONI
###############################################

results <- kimono(input_list, mapping_list, metainfo,  main_layer = 1, min_features = 2, stab_sel = F)
  
