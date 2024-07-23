# R CMD check results ──────────────────────────────────────── samplr 1.0.0 ────
```
❯ checking CRAN incoming feasibility ... [20s] NOTE
  Maintainer: 'Lucas Castillo <lucas.castillo-marti@warwick.ac.uk>'
  
  New submission
```

```
❯ checking C++ specification ... NOTE
    Specified C++11: please drop specification unless essential
```
We decide to keep the C++11 specification as it is essential. The package makes heavy use of lambda expressions, introduced in C++11, in the `getPDF()` function in the `src/pdf_manage.cpp` file. Without this C++11 feature, the package code would be much harder to read and maintain. 