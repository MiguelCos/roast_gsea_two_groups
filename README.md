# Limma-ROAST on a two-group dataset 

## Script to run ROAST-Limma on a two-group dataset

### Requirements

- It needs to be executed in a RProject environment (this is because of the use of the `here` package to define file locations).
- It requires the presence of the "msigdb.v7.0.symbols.gmt' file, which is the newest version of the MSIG database (available here:  http://software.broadinstitute.org/gsea/downloads.jsp. This file should be in the same RStudio Project folder from which this script would be executed.
- It requires an input TXT file with the next structure: 1 protein per rown, 1 column of Uniprot Protein IDs and 1 column per experimental condition and/or replicate. 
- The helper script should be also present in the same folder where the main script would be executed. 
