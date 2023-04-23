import pandas
import primer3


def get_data(df):

  mylist_1, mylist_2, mylist_3, mylist_4, mylist_5, mylist_6, mylist_7 = [], [], [], [], [], [], []
  print("Check point 4")
  for i in range(len(df["primer"])):
    mylist_1.append(round(primer3.calc_tm(df["primer"][i]),0))
    mylist_2.append(round(primer3.calc_tm(df["rightPrimers"][i]),0))
    mylist_3.append(round(primer3.calc_hairpin(df["primer"][i]).tm,0))
    mylist_4.append(round(primer3.calc_hairpin(df["rightPrimers"][i]).tm,0))
    mylist_5.append(round(primer3.bindings.calc_heterodimer(df["primer"][i],df["rightPrimers"][i]).dg/1000,1))
    mylist_6.append(round(primer3.bindings.calc_homodimer(df["primer"][i]).dg/1000,1))
    mylist_7.append(round(primer3.bindings.calc_homodimer(df["rightPrimers"][i]).dg/1000,1))

  df["TM_left (°C)" ] = mylist_1
  df["TM_right (°C)"] = mylist_2
  df["Hairpin_left (°C)"] = mylist_3
  df["Hairpin_right (°C)"] = mylist_4
  df["Heterodimer (kcal/mol)"] = mylist_5
  df["Homodimer_Left (kcal/mol)"] = mylist_6
  df["Homodimer_Right (kcal/mol)"] = mylist_7
  df.rename(columns={"primer": "leftPrimer (bp)"})
  df.rename(columns={"rightPrimers": "rightPrimer (bp)"})

  return(df)



