import pandas
import primer3
import numpy as np


def get_data(df):

  mylist_1, mylist_2, mylist_3, mylist_4, mylist_5, mylist_6, mylist_7, mylist_25 = [], [], [], [], [], [], [], []
  
  print("Python get data activated")
  for i in range(len(df["Forward"])):
    mylist_1.append(round(primer3.calc_tm(df["Forward"][i]),0))
    mylist_2.append(round(primer3.calc_tm(df["Reversed"][i]),0))
    mylist_25.append(abs(round(primer3.calc_tm(df["Forward"][i])-
                              primer3.calc_tm(df["Reversed"][i]))))
    mylist_3.append(round(primer3.calc_hairpin(df["Forward"][i]).tm,0))
    mylist_4.append(round(primer3.calc_hairpin(df["Reversed"][i]).tm,0))
    mylist_5.append(round(primer3.bindings.calc_heterodimer(df["Forward"][i],df["Reversed"][i]).dg/1000,1))
    mylist_6.append(round(primer3.bindings.calc_homodimer(df["Forward"][i]).dg/1000,1))
    mylist_7.append(round(primer3.bindings.calc_homodimer(df["Reversed"][i]).dg/1000,1))

  df["TM_left (°C)" ] = mylist_1
  df["TM_right (°C)"] = mylist_2
  df["TM_Diff (°C)"] = mylist_25
  df["Hairpin_left (°C)"] = mylist_3
  df["Hairpin_right (°C)"] = mylist_4
  df["Heterodimer (kcal/mol)"] = mylist_5
  df["Homodimer_Left (kcal/mol)"] = mylist_6
  df["Homodimer_Right (kcal/mol)"] = mylist_7
  
  df.rename(columns = {'Forward':'Forward (bp)'}, inplace = True)
  df.rename(columns = {'Reversed':'Reversed (bp)'}, inplace = True)


  return(df)


