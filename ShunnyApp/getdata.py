

def getdata (df):
  mylist_1, mylist_2, mylist_3, mylist_4, mylist_5, mylist_6 = [], [], [], [], [], []

  for i in range(len(df["primer"])):
    mylist_1.append(primer3.calc_tm(df["primer"][i]))
    mylist_2.append(primer3.calc_tm(df["rightPrimers"][i]))
    mylist_3.append(primer3.calc_hairpin(df["primer"][i]).tm)
    mylist_4.append(primer3.calc_hairpin(df["rightPrimers"][i]).tm)
    mylist_5.append(primer3.bindings.calc_end_stability(df["rightPrimers"][i], df["primer"][i]).dg)
    mylist_6.append(primer3.bindings.calc_end_stability(df["rightPrimers"][i], 
    df["leftPrimers"][i]).structure_found)
    
  df["TM_left"] = mylist_1
  df["TM_right"] = mylist_2
  df["Hairpin_left"] = mylist_3
  df["Hairpin_right"] = mylist_4
  df["Stability"] = mylist_5
  df["Structure"] = mylist_6
  return(df)
  
