# SNPselect

This dashboard select appropriate primers set for multiplexing of SNP DNA mutation within human.  body. The primers are used to build a metabolic panel for assessing genetic predisposition of metabolic diseases such as obesity, high blood pressure, and hyperglycemia. This R shiny served as a substitute for primer3 program itself by using both R and python to run the related packages and API (MART).

## References



Python wrapper (Mainly about the data programming structure): 
https://libnano.github.io/primer3-py/api/bindings.html#primer3.bindings.calc_heterodimer


Primer3 manual (I use it to how it calculates the functions):
https://primer3.org/manual.html#PRIMER_MAX_HAIRPIN_TH

Primer 3 in python module

https://github.com/jensenlab/primer3

## Adjustable Parameters

1. Mutation Upstream Length (bp)
2. Mutation Downstream Length (bp)
3. Distance Between Primers (bp)
4. Forward Primer Length (bp)
5. Reversed Primer Length (bp)
6. Maximum Forward Primer MT (c)
7. Maximum Reverse Primer MT (c)
6. Maximum Forward Primer Hairpin MT (c)
7. Maximum Reverse Primer Hairpin MT (c)
8. Maximum Difference Between TMs (c)
9. Minimum Homodimer delta G (Cal)
10. Minimum Heterodimer delta G (Cal)



## Installation

System Requirement:
1. R 4.2.3 or above
2. The code does not use any virtual environment
3. Previous installation of python is **not** required
4. **Do not use** minoconda since it is not compatible with shiny.io

Step 1 - Install the following R packages. uncomment (crl + shift + c) to run it in the document. Remeber to comment it back after you are done.

```bash
install.packages("reticulate")
install.packages("DT")
install.packages("shiny")
install.packages("spgs")
install.packages("rsconnect")
install.packages("dplyr")
install.packages("stringi")
install.packages("tidyverse")
install.packages("shinydashboard")
```

Step 2 - Run this chunk separately (also included in the code (but I am not sure what will happen if run everything together all at once)

```bash
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
```

Step 3 - Go to the console in R studio and install the following libraries 

```bash
pip install primer3-py
pip install pandas
```
## Data Table Interpretation
1. Identification - Description about where and how this set of primer is produced
2. Forward Primer - From 5' to 3'
3. Forward Primer - From 3' to 5'
4. Hairpin TM - This is the most stable monomer structure of internal oligo calculated by thermodynamic approach. The hairpin loops, bulge loops, internal loops, internal single mismatches, dangling ends, terminal mismatches have been considered. This parameter is calculated only if PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1. The default value is 10 degrees lower than the default value of PRIMER_MIN_TM
5. Heterodimer - The entropy of  this primer binds to another primer in the same set.
6. Homodimer - The entropy of this primer binds itself.
7. GC contents - yap

Step 4 - Restart the R session

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIIIIIIIT](https://choosealicense.com/licenses/mit/)
