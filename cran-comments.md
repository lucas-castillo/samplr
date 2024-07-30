# R CMD check results ──────────────────────────────────────── samplr 1.0.0 ────
```
❯ checking CRAN incoming feasibility ... [20s] NOTE
  Maintainer: 'Lucas Castillo <lucas.castillo-marti@warwick.ac.uk>'
  
  New submission
```
Thank you for reviewing this submission. 

## 29/07/2024 Comments 
Please omit the redundant "Tools for"/"A set of tools" at the estart of
your title and description.

> Done.

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year, ISBN:...)
or if those are not available: <[https:...]https:...>
with no space after 'doi:', 'https:' and angle brackets for
auto-linking. (If you want to add a title as well please put it in
quotes: "Title")

Please provide a link to the used webservices to the description field
of your DESCRIPTION file in the form
<http:...> or <[https:...]https:...>
with angle brackets for auto-linking and no space after 'http:' and
'https:'.

> Done. 

It seems like you have too many spaces in your description field.
Probably because linebreaks count as spaces too.
Please remove unecassary ones.

> Done. 

\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user. Does not seem necessary.
Please replace \dontrun with \donttest.

Please unwrap the examples if they are executable in < 5 sec, or replace
dontrun{} with \donttest{}.

> Thank you. These examples run in less than 5s so we've removed the \dontrun{} wrapper completely. 

Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited.
e.g.:
...
oldpar <- par(no.readonly = TRUE) # code line i
on.exit(par(oldpar)) # code line i + 1
...
par(mfrow=c(2,2)) # somewhere after
...
If you're not familiar with the function, please check ?on.exit. This
function makes it possible to restore options before exiting a function
even if the function breaks. Therefore it needs to be called immediately
after the option change within a function.
-> R/calc_functions.R