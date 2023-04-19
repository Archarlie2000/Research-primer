# SNPselect

This dashboard select appropriate primers set for multiplexing of SNP DNA mutation within human.  body. The primers are used to build a metabolic panel for assessing genetic predisposition of metabolic diseases such as obesity, high blood pressure, and hyperglycemia. This R shiny served as a substitute for primer3 program itself by using both R and python to run the related packages and API (MART).

## Adjustable Parameters

1. Mutation Upstream Length (bp)
2. Mutation Downstream Length (bp)
3. Distance Between Primers (bp)
4. Forward Primer Length (bp)
5. Reversed Primer Length (bp)
6. Filter - Maximum Forward Primer MT (c)
7. Filter - Reverse Forward Primer MT (c)
8. Filter - Maximum Difference Between TMs (c)
9. Filter - Minimum Homodimer delta G (Cal)
10. Filter - Minimum Heterodimer delta G (Cal)



## Installation

System Requirement:
1. python 3.11 or above
2. R 4.2.3 or above

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install foobar.

```bash
pip install primer3-py
pip install pandas
```

## Usage

```python
import foobar

# returns 'words'
foobar.pluralize('word')

# returns 'geese'
foobar.pluralize('goose')

# returns 'phenomenon'
foobar.singularize('phenomena')
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
