###############################################
###############################################
###### Minimal working example - KiMONo #######
###############################################
###############################################


###############################################
### 1 libraries
###############################################
#library(kimono)

library(data.table)


###############################################
### 1 Input Data
###############################################
# read in data
transcriptome <- fread("../kimono/example/data/expr.csv")
phenotype <- fread("../kimono/example/data/pheno.csv")

# make sure samples are in the same order
idorder <- transcriptome$sample
phenotype <- phenotype[match(idorder, phenotype$sample),]

# remove sample column
transcriptome <- transcriptome[,-"sample"]
phenotype <- phenotype[,-"sample"]

# input data list
# IMPORTANT - list element names MUST match mapping names!!
input_data <- list(
  'gene' = as.data.table(transcriptome),
  'phenotype' = as.data.table(phenotype)
)

#########################
# 2 Prior Network
# generated from mapping files

# known gene-gene mappings
# IMPORTANT - load_mapping
gene_gene <- fread("../kimono/example/data/mapping_expr.csv")
transcriptome_map <- load_mapping(gene_gene, layers = c('gene','gene'))

# known gene-phenotype relations
gene_phenotype <- fread("../kimono/example/data/mapping_expr_pheno.csv")
phenotype_map <- load_mapping(gene_phenotype, layers = c('gene','phenotype'))

# combine mappings
prior_map <- rbind(transcriptome_map,phenotype_map)
prior_network <- create_prior_network(prior_map)

###############################################
# 4 Run MONI
###############################################

results <- kimono(input_data, prior_network, min_features = 0,core = 1 )
