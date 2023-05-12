import streamlit as st
import rpy2.robjects as robjects

st.write(robjects.r("pi")[0])
st.write(robjects.r("1+1")[0])