require(lattice)
require(magick)
# -----
#
# solucao bidimensional esquema ADI, analitica
# p.379 livro Nelson
#
# ----
analitica <- function(D,L,uo,x,y,tempo){
  t <- tempo*0.01
  exp1 <- uo/(1+4*t*D/L^2)
  exp2 <- exp(-(x^2+y^2)/(L^2+4*D*t))
  return(exp1*exp2)
}
# valores atribuidos
D = 0.1           # difusividade
L = 1             # comprimento L
dt = 0.01         # intervalo de tempo
delta = 1/55      # intervalo no espaco
uo = 1            # velocidade
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
# visualiza solucao analitica
levelplot(t(ana[,,1]))
levelplot(t(ana[,,2]))
levelplot(t(ana[,,3]))
levelplot(t(ana[,,4]))
levelplot(t(ana[,,5]))
levelplot(t(ana[,,6]))
levelplot(t(ana[,,7]))
levelplot(t(ana[,,8]))
levelplot(t(ana[,,9]))
levelplot(t(ana[,,15]))
levelplot(t(ana[,,26]))
levelplot(t(ana[,,37]))
levelplot(t(ana[,,48]))
levelplot(t(ana[,,59]))
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
break
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
  tt <- tt+1
}
    
break
levelplot(t(v[,,1]))
levelplot(t(v[,,2]))
levelplot(t(v[,,3]))
levelplot(t(v[,,4]))
levelplot(t(v[,,5]))
levelplot(t(v[,,6]))
levelplot(t(v[,,7]))
levelplot(t(v[,,18]))
levelplot(t(v[,,29]))
levelplot(t(v[,,35]))
levelplot(t(v[,,46]))
levelplot(t(v[,,67]))
levelplot(t(v[,,88]))
levelplot(t(v[,,99]))


p1 <- levelplot(abs(ana[,,1]-v[,,1]),xlab="x",ylab="y",at=seq(0,0.3,length.out = 100),col.regions=heat.colors(100)[100:1],
                main=paste0("dt = ",1))
p2 <- levelplot(abs(ana[,,12]-v[,,12]),xlab="x",ylab="y",at=seq(0,0.3,length.out = 100),col.regions=heat.colors(100)[100:1],
                main=paste0("dt = ",12))
p3 <- levelplot(abs(ana[,,24]-v[,,24]),xlab="x",ylab="y",at=seq(0,0.3,length.out = 100),col.regions=heat.colors(100)[100:1],
                main=paste0("dt = ",24))
p4 <-levelplot(abs(ana[,,36]-v[,,36]),xlab="x",ylab="y",at=seq(0,0.3,length.out = 100),col.regions=heat.colors(100)[100:1],
               main=paste0("dt = ",36))
p5 <- levelplot(abs(ana[,,48]-v[,,48]),xlab="x",ylab="y",at=seq(0,0.3,length.out = 100),col.regions=heat.colors(100)[100:1],
                main=paste0("dt = ",48))
p6 <-levelplot(abs(ana[,,60]-v[,,60]),xlab="x",ylab="y",at=seq(0,0.3,length.out = 100),col.regions=heat.colors(100)[100:1],
               main=paste0("dt = ",60))
grid.arrange(p1,p2,p3,p4,p5,p6)
levelplot(t(v[,,5]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",5))
levelplot(t(v[,,6]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",6))
levelplot(t(v[,,7]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",7))
levelplot(t(v[,,8]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",8))
levelplot(t(v[,,9]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",9))
levelplot(t(v[,,10]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",10))
levelplot(t(v[,,12]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",12))
levelplot(t(v[,,35]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",35))
levelplot(t(v[,,70]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*70))
levelplot(t(v[,,80]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*80))
levelplot(t(v[,,90]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*90))
levelplot(t(v[,,100]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*100))
levelplot(t(v[,,200]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*200))
levelplot(t(v[,,500]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*500))
levelplot(t(v[,,600]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*600))
levelplot(t(v[,,700]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*700))
levelplot(t(v[,,800]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*800))
levelplot(t(v[,,900]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*900))
levelplot(t(v[,,1000]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*1000))

animation <- image_animate(img, fps = 2)
image_write(animation,"E:/Documents/CICLO/Mestrado/materias/met_mat_eng/trabalho2/difusao2d-in.gif")
fo
