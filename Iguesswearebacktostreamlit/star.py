
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import streamlit as st
import pandas as pd
import numpy as np


Bioc = rpackages.importr('BiocManager')
BioM = rpackages.importr('biomaRt')


st.title("Add One to Your Number")

# Create a numeric input field
input_number = st.number_input("Enter a number:")

# Add 1 to the input value and display the result
if input_number:
    output_number = input_number + 1
    st.write("Input value:", input_number)
    st.write("Output value:", output_number)
    
# robjects.r('''
#     add_nums <- function(x, y) {
#         return(x + y)
#     }


#     snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

#     snp_sequence <- getBM(attributes = c('refsnp_id', 'snp'),
#       filters = c('snp_filter', 'upstream_flank', 'downstream_flank'),
#       checkFilters = FALSE,
#       values = list("rs25",  500, 500),
#       mart = snpmart,
#       bmHeader = TRUE)

#       print(snp_sequence)

#     print(add_nums(x = 5, y = 10))
#     print(add_nums(x = 10, y = 20))
# ''')


