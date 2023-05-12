import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.2.3" # change as needed

import streamlit as st
import rpy2.robjects as robjects

st.write(robjects.r("pi")[0])
st.write(robjects.r("1+1")[0])