# Character-Polynomial Codes

This is an implementation of character-polynomial (CP) codes using [Sage](https://doc.sagemath.org/html/en/reference/index.html). 
The following classes are implemented.

## `GRSBase`

This represents the GRS subcode corresponding to the CP code. 
The actual code is returned by its `grs_subcode()` method, which returns a 
[`LinearCode`](https://doc.sagemath.org/html/en/reference/coding/sage/coding/linear_code.html#sage.coding.linear_code.LinearCode). 

## `CharacterPolynomialCode`

This is the CP code class. Currently, it only supports CP codes over prime fields. 

See [`cpcode.sage`](https://github.com/nayel71/cpcode/blob/main/cpcode.sage) for documentations and 
[`test.ipynb`](https://github.com/nayel71/cpcode/blob/main/test.ipynb) for example usages of both classes. 
