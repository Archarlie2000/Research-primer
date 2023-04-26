
import pandas
import numpy as np


def get_filter(df, left_TM, right_TM, 
                       left_hair_TM, 
                       right_hair_TM, 
                       diff,
                       Homodimer_left_dg, 
                       Homodimer_right_dg, 
                       Heterodimer_dg):
                         
  print("Python get fiiiiiiiiiilter activated")
  df = df[df['TM_left (°C)'] <= left_TM]
  df = df[df["TM_right (°C)"] <= right_TM]
  df = df[df["Hairpin_left (°C)"] <= left_hair_TM]
  df = df[df["Hairpin_right (°C)"] <= right_hair_TM]
  df = df[df["TM_Diff (°C)"] <= diff]
  df = df[df['Heterodimer (kcal/mol)'] <= Homodimer_left_dg]
  df = df[df['Homodimer_Left (kcal/mol)'] <= Homodimer_right_dg]
  df = df[df['Homodimer_Right (kcal/mol)'] <= Heterodimer_dg]
  
  return(df)
