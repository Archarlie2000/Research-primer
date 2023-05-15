# SNPselect

This dashboard select appropriate primers set for multiplexing of SNP DNA mutation within human.  body. The primers are used to build a metabolic panel for assessing genetic predisposition of metabolic diseases such as obesity, high blood pressure, and hyperglycemia. This R shiny served as a substitute for primer3 program itself by using both R and python to run the related packages and API (MART).

Web version: https://vw8ids-archarlie-chou.shinyapps.io/Clone_move_on/

## References

Python wrapper (Mainly about the data programming structure): 
https://libnano.github.io/primer3-py/api/bindings.html#primer3.bindings.calc_heterodimer


Primer3 manual (I use it to how it calculates the functions):
https://primer3.org/manual.html#PRIMER_MAX_HAIRPIN_TH

Primer 3 in python module  https://github.com/jensenlab/primer3

## User Manual

Section: DashBoard

1. Enter SNP - It takes it multple SNP ID. Which can be found in https://www.ncbi.nlm.nih.gov/snp/docs/RefSNP_about/
2. Amplicant Length (°C) - the distance of the reversed primer will be based on
3. Shift (bp) - For each shift, the amplicant length will be shorten by 1, and a new batch of new primers will be produced
4. Forward (bp) - Number of base pair in forward primer
5. Reversed (bp) - Number of base pair in reverd primer
6. Left TM max (°C) - This TM result should be identical to primer3. We assume primer concentration = 50nM, Na+ = 50nM
7. Right TM max (°C) - Tm of reverse primer
8. Left hairpin max (°C) - This TM result should be identical to primer3
9. Right hairpin max (°C) - This TM result should be identical to primer3
10. Max difference in Tm (°C) - abs( Tm forward - Tm reversed )
11. Homodimer left (°C) - The melting temperature of forward primer to forward primer
12. Homodimer right (°C) - The melting temperature of reversed primer to reversed primer
13. Heterodimer (°C) - The melting temperature of forward primer to reversed primer

Section Analysis:

Press download to , well, dowload

## Installation

System Requirement:
1. R 4.2.3 or above

Step 1 - Install the following R packages. uncomment (crl + shift + c) to run it in the document. Remeber to comment it back after you are done.

```bash
install.packages("DT")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("stringi")
install.packages("ggplot2")
install.packages("hexbin")
install.packages("patchwork")
install.packages("plotly")
install.packages("devtools")
devtools::install_github("jensenlab/primer3")
install.packages("TmCalculator")
install.packages("BiocManager")
BiocManager::install("biomaRt")
install.packages("spgs")
install.packages("shiny")
install.packages("rsconnect")
install.packages("shinydashboard")
```





## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIIIIIIIT](https://choosealicense.com/licenses/mit/)
