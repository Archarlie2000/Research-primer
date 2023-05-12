import streamlit as st
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

st.write(robjects.r("pi")[0])
st.write(robjects.r("1+1")[0])