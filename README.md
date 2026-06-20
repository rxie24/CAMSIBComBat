**Maintainer**: Ryan Xie, ryan.xie@pennmedicine.upenn.edu

**License**: Artistic License 2.0

## Overview

The **CAMSIBComBat** R package contains the **cam_combat** and **sib_combat** functions which implement Covariance-Aware Multivariate (CAM) ComBat and and Spatially-Informated Iterative Block (SIB) ComBat, respectively. These are methods for harmonizing neuroimaging data.

## Installation

The package can be installed from GitHub as follows:

```r
# Install devtools package if not already installed
install.packages("devtools")

# Install package
devtools::install_github("rxie24/CAMSIBComBat")
```

The package can then by loaded in R as follows:

```r
library(CAMSIBComBat)
```

Example of implementation of each function can be found in the package documentation.

## Citation (APA)

Xie, R., Srinivasan, D., Harman, G. A., Davatzikos, C., Shinohara, R. T., & Shou, H. (2026). Batch Effect Correction for Neuroimaging Data with Heterogeneous Spatial Correlations. bioRxiv, 2026-06.
