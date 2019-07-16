library(Ryacas);
x <- Sym("x");
s <- expression(a*x^2);
deriv(s,x)
der <- deriv(s,x);
x <- 1:2
ev <- eval(der)
ev
g(1)
#
#
a1 <- expression(1*x)
b1 <- expression(2*x)
as.expression(a1+b1)
f <- function(x){eval(a1[[1]])+eval(b1[[1]])}
f(1)
#
deriv3(expr = s,namevec = x);
d <- D(s,x)
typeof(d)

f <- function(x) (a*x^2)
g <- function(x) {}
body(g) <- D(body(f), 'x')
h <- function(x){}
body(h) <- D(body(g),'x')
g
h
i <- function(x){
  f(x)+g(x)
}

a <- expression(a*x);
b <- expression(b*x);
c <- expression(a+b)
c
