import pandas
import primer3
import numpy as np


def get_data(df, left_TM, 
                       right_TM, 
                       left_hair_TM, 
                       right_hair_TM, 
                       diff,
                       Homodimer_left_dg, 
                       Homodimer_right_dg, 
                       Heterodimer_dg):

  mylist_1, mylist_2, mylist_3, mylist_4, mylist_5, mylist_6, mylist_7, mylist_25 = [], [], [], [], [], [], [], []
  
  print("Check point 4")
  for i in range(len(df["primer"])):
    mylist_1.append(round(primer3.calc_tm(df["primer"][i]),0))
    mylist_2.append(round(primer3.calc_tm(df["rightPrimers"][i]),0))
    mylist_25.append(abs(round(primer3.calc_tm(df["primer"][i])-
                              primer3.calc_tm(df["rightPrimers"][i]))))
    mylist_3.append(round(primer3.calc_hairpin(df["primer"][i]).tm,0))
    mylist_4.append(round(primer3.calc_hairpin(df["rightPrimers"][i]).tm,0))
    mylist_5.append(round(primer3.bindings.calc_heterodimer(df["primer"][i],df["rightPrimers"][i]).dg/1000,1))
    mylist_6.append(round(primer3.bindings.calc_homodimer(df["primer"][i]).dg/1000,1))
    mylist_7.append(round(primer3.bindings.calc_homodimer(df["rightPrimers"][i]).dg/1000,1))

  df["TM_left (°C)" ] = mylist_1
  df["TM_right (°C)"] = mylist_2
  df["TM_Diff (°C)"] = mylist_25
  df["Hairpin_left (°C)"] = mylist_3
  df["Hairpin_right (°C)"] = mylist_4
  df["Heterodimer (kcal/mol)"] = mylist_5
  df["Homodimer_Left (kcal/mol)"] = mylist_6
  df["Homodimer_Right (kcal/mol)"] = mylist_7

  
  df.rename(columns = {'primer':'Forward (bp)'}, inplace = True)
  df.rename(columns = {'rightPrimers':'Reversed (bp)'}, inplace = True)
  
  df = df[df['TM_left (°C)'] <= left_TM]
  df = df[df["TM_right (°C)"] <= right_TM]
  df = df[df["Hairpin_left (°C)"] <= left_hair_TM]
  df = df[df["Hairpin_right (°C)"] <= right_hair_TM]
  df = df[df["TM_Diff (°C)"] <= diff]
  df = df[df['Heterodimer (kcal/mol)'] <= Homodimer_left_dg]
  df = df[df['Homodimer_Left (kcal/mol)'] <= Homodimer_right_dg]
  df = df[df['Homodimer_Right (kcal/mol)'] <= Heterodimer_dg]


  return(df)


