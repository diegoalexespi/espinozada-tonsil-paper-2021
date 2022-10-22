require(readxl)
require(magrittr)
require(dplyr)
require(biomaRt)

blood.data <- read_excel("zh2240-sup-tables4.xlsx", skip = 4)
blood.data %>% mutate(affy_hg_u133_plus_2 = `Probe Set ID`) -> blood.data

ensembl.db <- useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl")
ensembl.db.table <- getBM(attributes=c('affy_hg_u133_plus_2', 'hgnc_symbol', 'ensembl_gene_id', 'external_gene_name'), filters ='affy_hg_u133_plus_2', values = blood.data$`Probe Set ID`, mart = ensembl.db) %>% mutate(consensus_name = ifelse(hgnc_symbol == "", external_gene_name, hgnc_symbol))

full_join(blood.data, ensembl.db.table, by = "affy_hg_u133_plus_2") %>% mutate(final_gene_name = ifelse(is.na(consensus_name), `Gene Symbol`, consensus_name)) %>% dplyr::select(final_gene_name, `Z-score`, `Enriched in`, affy_hg_u133_plus_2) -> final.table

write.table(x = final.table, file = "blood_DZ_LZ_genes.txt", sep = '\t')
