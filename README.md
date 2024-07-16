# Studying the Association Between Longitudinal Non-Dense Breast Tissue and Breast Cancer Risk: A Joint Modelling Approach

This repository contains the material for the manuscript titled _Studying the Association Between Longitudinal Non-Dense Breast Tissue and Breast Cancer Risk: A Joint Modelling Approach_.
This manuscript is currently in press at the [American Journal of Epidemiology](https://academic.oup.com/aje) and can be accessed at the following DOI: <https://doi.org/10.1093/aje/kwae196>

The content of this repository is organised as follows.

* The `programs/` folder includes R code that could be used to fit a univariate and a bivariate joint model using the {VAJointSurv} package.
  The two joint models include delayed entry, markers with a non-linear effect of time, and assume the expected value association structure; additional covariates are included in the survival submodel as well. 

* The `data/` folder includes synthetic data that could be used together with the R code to fit these joint models; data is included in both `.csv` and R serialised (`.RDS`) formats. 
  Note that these synthetic datasets should not be used for any scientific purpose, as they do not represent any real data and are provided for illustrative purposes only.

* The `data-raw/` folder includes R code that was used to generate the synthetic datasets, if you'd like to double-check how the data was generated.
  Note that this code is not reproducible unless you have access to the study data, as we inform the simulation parameters on real data features to obtain a realistic synthetic dataset.

If you have any questions regarding the code or data, please [file an issue](https://github.com/ellessenne/nond-bc-jm/issues).

## License

The code from this repository is licensed under the [MIT license](https://github.com/ellessenne/nond-bc-jm/blob/main/LICENSE).
All other included content is released under a [Creative Commons Zero v1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) license.
