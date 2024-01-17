Z_equalities <- function(
    a=NULL, 
    b=NULL, 
    a_and_b=NULL, 
    a_or_b=NULL, 
    a_given_b=NULL,
    b_given_a=NULL,
    a_given_not_b=NULL,
    b_given_not_a=NULL,
    a_and_not_b=NULL,
    b_and_not_a=NULL,
    not_a=NULL,
    not_b=NULL
    ){
  # Check inputs
  is_prob <- function(x){
    x <- x[!is.na(x)]
    if (any(x<0)) stop("Probabilites cannot be negative")
    if (any(x>1)) stop("Probabilites cannot be larger than 1")
  }
  for (p in c(a,
             b,
             not_a,
             not_b,
             a_and_b,
             a_or_b,
             a_given_b,
             b_given_a,
             a_given_not_b,
             b_given_not_a,
             a_and_not_b,
             b_and_not_a)){
    if (!is.null(p)){is_prob(p)}
  }
  
  # add complements if not given
  if (is.null(not_a)){
    not_a = 1-a
  }
  if (is.null(not_b)){
    not_b = 1-b
  }
  
  safe_eval <- function(exp){
    if (length(eval(exp)) == 0) return(NA) else return(eval(exp))
  }
  
  data.frame(
    z1 = safe_eval(expression(a + b - a_and_b - a_or_b)),
    z2 = safe_eval(expression(a + b_and_not_a - b - a_and_not_b)),
    z3 = safe_eval(expression(a + b_and_not_a - a_or_b)),
    z4 = safe_eval(expression(b + a_and_not_b - a_or_b)),
    z5 = safe_eval(expression(a_and_not_b + a_and_b - a)),
    z6 = safe_eval(expression(b_and_not_a + a_and_b - b)),
    z7 = safe_eval(expression(a_and_not_b + b_and_not_a + a_and_b - a_or_b)),
    z8 = safe_eval(expression(a_and_not_b + b_and_not_a + 2*a_and_b - a - b)),
    z9 = safe_eval(expression(a_given_b*b - b_given_a*a)),
    z10= safe_eval(expression(a_given_b*b + a_given_not_b*not_b - a)),
    z11= safe_eval(expression(b_given_a*a + b_given_not_a*not_a - b)),
    z12= safe_eval(expression(b_given_a*a + a_given_not_b*not_b - a)),
    z13= safe_eval(expression(a_given_b*b + b_given_not_a*not_a - b)),
    z14= safe_eval(expression(a_given_not_b*not_b + b - b_given_not_a*not_a - a)),
    z15= safe_eval(expression(a_and_b - a_given_b*b)),
    z16= safe_eval(expression(a_and_b - b_given_a*a)),
    z17= safe_eval(expression(a_and_b - a + a_given_not_b*not_b)),
    z18= safe_eval(expression(a_and_b - b + b_given_not_a*not_a))
  )
}

