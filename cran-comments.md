## R CMD check results
0 errors ✔ | 0 warnings ✔ | 2 notes ✖

── R CMD check results ──────────────────────────────────────────── samplr 0.0.0.9000 ────
Duration: 2m 34s

❯ checking examples ... ERROR
  Running examples in 'samplr-Ex.R' failed
  The error most likely occurred in:
  
  > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
  > ### Name: Mean_Variance
  > ### Title: Mean Variance Estimates
  > ### Aliases: Mean_Variance
  > 
  > ### ** Examples
  > 
  > library(dplyr)
  
  Attaching package: 'dplyr'
  
  The following objects are masked from 'package:stats':
  
      filter, lag
  
  The following objects are masked from 'package:base':
  
      intersect, setdiff, setequal, union
  
  > library(tidyr)
  > library(magrittr)
  
  Attaching package: 'magrittr'
  
  The following object is masked from 'package:tidyr':
  
      extract
  
  > library(samplrData)
  > data <- sundh2023e3 %>%
  +   group_by(ID, querydetail) %>% 
  +   mutate(iteration = LETTERS[1:n()]) %>% 
  +   pivot_wider(id_cols = c(ID, querydetail), 
  +       values_from = estimate, names_from = iteration) %>% 
  +   mutate(across(where(is.numeric), \(x){x/100})) %>% 
  +   ungroup %>% 
  +   select(-querydetail)
  Error: object 'sundh2023e3' not found
  Execution halted

❯ checking C++ specification ... NOTE
    Specified C++11: please drop specification unless essential

