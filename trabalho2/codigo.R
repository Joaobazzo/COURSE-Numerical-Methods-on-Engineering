require(lattice)
require(magick)
# -----
#
# solucao bidimensional esquema ADI, analitica
# p.379 livro Nelson
#
# Joao Pedro Bazzo
#
# ----
#
# funcao solucao analitica
analitica <- function(D,L,uo,x,y,tempo){
  t <- tempo*0.01
  exp1 <- uo/(1+4*t*D/L^2)
  exp2 <- exp(-(x^2+y^2)/(L^2+4*D*t))
  return(exp1*exp2)
}
# valores atribuidos
D = 0.020
L = 1
dt = 0.01
delta = 1/35
uo = 1
fo = dt*D/(delta^2)
fo
nt = 1/dt
nd = 1/delta
#nd=20
print(paste0("nd= ",nd))
# matrix sol. numerica
u <- array(data = NA,dim = c(nd,nd,nt))
#
#
# condicoes iniciais
xi <- seq(-delta*nd,+delta*nd,length.out = nd)
yi <- seq(-delta*nd,+delta*nd,length.out = nd)
for(j in (1:nd)){
  for(i in (1:nd)){
    u[i,j,1] <- analitica(D,L,uo,xi[i],yi[j],t = 0)
  }
}
#
# solucao analitica
ana <- u
for(tt in (1:nt)){
  for(i in (1:nd)){
    for(j in (1:nd)){
      ana[i,j,tt]  <- analitica(D = D,L = L,uo = uo,x = xi[i], y=yi[j],tempo = tt-1)
    }
  }
}
#
# condicoes de contorno u(L,y,t) = u(-L,y,t) = u(x,-L,t)= u(x,L,t)
for(tt in (1:nt)){
  u[1,1:nd,tt]  <- analitica(D,L,uo,x = xi[1], y=yi[1:nd],tt-1)
  u[nd,1:nd,tt] <- analitica(D,L,uo,x = xi[nd],y=yi[1:nd],tt-1)
  
  u[1:nd,1,tt]  <- analitica(D,L,uo,x = xi[1:nd], y=yi[1],tt-1)
  u[1:nd,nd,tt] <- analitica(D,L,uo,x = xi[1:nd], y=yi[nd],tt-1)
}
#
# carrega matrix A
A <- matrix(ncol=nd,nrow=nd)
for(i in (1:(nd))){
  if(i==1){
    A[i,]<-c(1+2*fo,-fo,
             rep(0,nd-2))}
  if(i==nd){
    A[i,] <- c(rep(0,nd-2),
               -fo,1+2*fo)}
  if(i>1&&i<nd){
    A[i,]<-c(rep(0,i-2),
             -fo,1+2*fo,-fo,
             rep(0,nd-1-i))}
}
# loop no espaco x,y
v <- u
#break
#j=2
tt=1
while(tt <= (nt-2)){ # tempo
  # -
  # espaço 'x'
  # -
  for(j in (2:(nd-1))){ 
    B <- matrix(ncol=1,nrow=nd)
    # b[1]
    B[1] <- fo*analitica(D,L,uo,x = xi[1]-diff(xi)[1],y=yi[j],tt+1)+
      fo*v[1,j+1,tt]+(1-2*fo)*v[1,j,tt]+fo*v[1,j-1,tt]
    # b[inter]
    B[2:(nd-1)] <- fo*v[2:(nd-1),j+1,tt]+(1-2*fo)*v[2:(nd-1),j,tt]+
      fo*v[2:(nd-1),j-1,tt]
    # b[n]
    B[nd] <- fo*analitica(D,L,uo,x=xi[nd]+diff(xi)[1],y = yi[j],tt+1)+
      fo*v[nd,j+1,tt]+(1-2*fo)*v[nd,j,tt]+fo*v[nd,j-1,tt]
    #
    v[1:nd,j,tt+1] <- solve(A,B)
  }
  # -
  # espaco 'y'
  # -
  # matrix B
  B <- matrix(ncol=1,nrow=nd)
  for(i in (2:(nd-1))){ 
    # b[1]
    B[1] <- fo*analitica(D,L,uo,x=xi[i],y = yi[1]-diff(yi)[1],tt+2)+
      fo*v[i-1,1,tt+1]+(1-2*fo)*v[i,1,tt+1]+fo*v[i+1,1,tt+1]
    # b[inter]
    B[2:(nd-1)] <- fo*v[i-1,2:(nd-1),tt+1]+
      (1-2*fo)*v[i,2:(nd-1),tt+1]+
      fo*v[i+1,2:(nd-1),tt+1]
    # b[n]
    B[nd] <- fo*analitica(D,L,uo,xi[i],y = yi[nd]+diff(yi)[1],tt+2)+
      fo*v[i-1,nd,tt+1]+
      (1-2*fo)*v[i,nd,tt+1]+
      fo*v[i+1,nd,tt+1]
    #
    v[i,1:nd,tt+2] <- solve(A,B)
    #
    #
  }
  tt <- tt+2
}


