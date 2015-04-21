

# Using scater for analysis -----------------------------------------------

library(scater)
library(devtools)
library(data.table)
library(devtools)

embeddr_path <- '~/oxford/embeddr/embeddr' # "/net/isi-scratch/kieran/embeddr/embeddr"
load_all(embeddr_path)

data_path <- '~/oxford/embeddr/data//data.txt' # "/net/isi-scratch/kieran//datasets/lung/data.txt"
d <- fread(data_path, data.table=FALSE)

bulk_i <- which(d$putative_cell_type == 'bulk')
d <- d[-bulk_i,]

d <- dplyr::filter(d, putative_cell_type != 'bulk')

d_numeric <- select(d, -cell_name, -time_point, -sample, -putative_cell_type)
d_pheno <- select(d, cell_name, time_point, sample, putative_cell_type)

pd <- new('AnnotatedDataFrame', data=d_pheno)
sce <- newSCESet(cellData = t(d_numeric), phenoData = pd)

