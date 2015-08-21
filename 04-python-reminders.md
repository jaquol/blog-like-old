# 04. Python remindes



## Get variables from script arguments

```
from sys import argv

var1, var2 = argv

# Alternatively, for many variables...
var1 = argv[0]
var2 = argv[1]
```


## Multi-lines printing

```
print """
There's something going on here.
With the three double-quotes.
We'll be able to type as much as we like.
Even 4 lines if we want, or 5, or 6.
"""
```


## Formatting strings with many formats

```
formatter = "%s %i %f %r"
print formatter % ('hello', 2, 3.00, 4)
```
