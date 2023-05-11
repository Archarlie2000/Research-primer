import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import streamlit as st

Bioc = rpackages.importr('BiocManager')
BioM = rpackages.importr('biomaRt')



st.title("Add One to Your Number")

# Create a numeric input field
input_number = st.number_input("Enter a number:")


output_number = robjects.r['pi']

output_number = robjects.conversion.rpy2py(output_number)

# Add 1 to the input value and display the result
if input_number:
    st.write("Input value:", input_number)
    st.write("Output value:", output_number)


