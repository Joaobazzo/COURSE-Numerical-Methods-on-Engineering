library(Ryacas)
comx <- Sym("x");
g_1 <- expression(1+
                    a*(2*(cos(x)-1))+
                    c*((4*cos(x)-cos(2*x)-3)/8)-
                    ci*(sin(x)+(2*sin(x)-sin(2*x))/8))

deriv(g_1,x)
deriv(g_a,x)
deriv(g_c,x)
deriv(g_ci,x)

g3 <- expression((1-4*a-c)^2-9*c^2/2-(1-4*a-c)*(2*a+c/2)-
                   ((1-4*a-c)^2-9/4*c^2)^(3/2))
g3
