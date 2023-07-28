# SNPselect

This dashboard select appropriate primers set for multiplexing of SNP DNA mutation within human.  body. The primers are used to build a metabolic panel for assessing genetic predisposition of metabolic diseases such as obesity, high blood pressure, and hyperglycemia. This R shiny served as a substitute for primer3 program itself.

Web version: https://vw8ids-archarlie-chou.shinyapps.io/Clone_move_on/
## Attention

Welcome to the documentation for the multiplexing algorithm used in my project. As you delve into this codebase, I want to highlight some important issues that you should be aware of and consider addressing:

Validation and Parameter Adjustment:
It's important to note that this program has not been validated extensively. The accuracy and reliability of the SNP location are still uncertain. However, I have set up various parameters that allow for systematic adjustment by modifying a single variable. Make sure to thoroughly validate and verify the SNP positions. Additionally, pay attention to the functions responsible for string wrangling and ensure that they are in the correct position. Due to time constraints, I couldn't thoroughly test their placement, so it's crucial to confirm their accuracy.

Tree Search Implementation:
Currently, the algorithm only utilizes the primers from the left flanking when implementing the tree search. This decision was made to avoid the complexity of growing multiple trees when incorporating both left and right flanking regions. However, it would be beneficial to explore efficient methods for streamlining this process and incorporating primers from both flanking regions. Consider the implications and potential improvements to enhance the overall efficiency and accuracy of the algorithm.

Incomplete Filter Application:
The filtering process in the current implementation is not fully comprehensive. When filtering out invalid hairpin structures and Tms, you may end up with only one remaining SNP. However, multiplexing requires multiple SNPs to be present. In the code, I have temporarily commented out some of the filtering steps to gather more data. It's essential to revisit this issue and develop a robust filtering mechanism that ensures valid SNPs are retained while maintaining an adequate number for multiplexing.



## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIIIIIIIT](https://choosealicense.com/licenses/mit/)
