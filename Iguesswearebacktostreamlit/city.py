import streamlit as st
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import numpy
from rpy2.robjects import r, pandas2ri
import pandas as pd
pandas2ri.activate()



Bioc = rpackages.importr('BiocManager')
BioM = rpackages.importr('biomaRt')



st.title("Add One to Your Number")

# Create a numeric input field
name = st.text_input('Enter SNP ID: ', "rs1121980, rs9939609, rs7903146, rs7903146")




robjects.r('''

        f <- function(name) {
    
            snp_list <- strsplit(name, " ")[[1]]

            snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
            snp_sequence <- getBM(attributes = c('refsnp_id', 'snp'),
            filters = c('snp_filter', 'upstream_flank', 'downstream_flank'),
            checkFilters = FALSE,
            values = list(snp_list,  500, 500),
            mart = snpmart,
            bmHeader = TRUE)
        }

        ''')

r_f = robjects.r['f']
res  = pd.DataFrame(list(r_f(name)[1]))

st.write('The snp sequnce ', res)
