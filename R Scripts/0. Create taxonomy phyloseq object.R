# This script will combine our QIIME2 outputs to create a taxonomy phyloseq object.

library(tidyverse)
library(phyloseq)

# Load the data
tax = read.delim('Datasets/taxonomy.tsv', header = TRUE, row.names = 1)
otu = read.delim('Datasets/feature-table.txt', skip=1, header = TRUE, row.names = 1) # skip first line
meta = read.delim('Datasets/sample-metadata.tsv', header = TRUE, row.names = 1)
meta = meta[-1, ] # Remove the first row, which does not contain data
tree = read_tree('Datasets/tree.nwk')

# Format taxonomy and count data for integration into phyloseq
# Metadata and tree are already in the correct format
tax_formatted = tax %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  select(-Confidence) %>% 
  as.matrix()

otu_formatted = otu %>% as.matrix()

# Combine into phyloseq object
ps = phyloseq(otu_table(otu_formatted, taxa_are_rows = TRUE),
              tax_table(tax_formatted),
              sample_data(meta),
              phy_tree(tree))

# Save formatt4d phyloseq object
saveRDS(ps, file = 'Datasets/phyloseq_taxonomy.rds')
