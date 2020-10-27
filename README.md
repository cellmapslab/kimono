# KiMONo 

Knowledge guided multi-level network inference

## About
KiMONo is a network inference tool for multi-omics datasets. 
The multi-omics network is calculated by sparse-group-LASSO regression and can incorporate prior biological information using protein-protein Interactions.

## Installation

You can either install KiMONo locally by cloning the repository or using the devtools package.
#### Local Installation

1. In your terminal change the working directory to the location you want to install KiMONo 

2. Clone the repository:
```sh
git clone https://github.com/cellmapslab/kimono.git
```

3. Install KiMONo in R and load the package 
```R
install.packages("yourpath/kimono/", repos = NULL, type = "source")
library(kimono)
```



#### Github Installation

1. Install the devtools package and load it in R

```R
install.packages("devtools")
library(devtools)
```
2. Install KiMONo in R and load the package
```R
install_github("cellmapslab/kimono")
library(kimono)
```

### Installing dependencies CentOS
#### oem
 CentOS needs a different version of RcppArmadillo(https://www.gitmemory.com/RcppCore) -> install.packages("RcppArmadillo", repos="https://rcppcore.github.io/drat")
