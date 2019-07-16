require(lattice)
require(magick)
# -----
#
# solucao bidimensional esquema ADI, analitica
# p.379 livro Nelson
#
# ----
analitica <- function(D,L,uo,x,y,t){
  exp1 <- uo/(1+4*t*D/L^2)
  exp2 <- exp(-(x^2+y^2)/(L^2+4*D*t))
  return(exp1*exp2)
}
# valores atribuidos
D = 1
L = 1
dt = 0.01
delta = 1/6
uo = 1
fo = dt*D/delta
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
    u[,j,1] <- analitica(D,L,uo,xi[1:nd],yi[j],t = 0)
}
# plot condicoes iniciais
levelplot(u[,,1])
#plot3d(x=0:49,y=0:49,z=u[,,1],type = "l")
#
# condicoes de contorno u(L,y,t) = u(-L,y,t) = u(x,-L,t)= u(x,L,t)
for(tt in (1:nt)){
    u[1,1:nd,tt]  <- analitica(D,L,uo,x = xi[1], y=yi[1:nd],t=tt-1)
    u[nd,1:nd,tt] <- analitica(D,L,uo,x = xi[nd],y=yi[1:nd],t=tt-1)
    
    u[1:nd,1,tt]  <- analitica(D,L,uo,x = xi[1:nd], y=yi[1],t=tt-1)
    u[1:nd,nd,tt] <- analitica(D,L,uo,x = xi[1:nd], y=yi[nd],t=tt-1)
}
#plot(u[1,,35])
levelplot(u[,,1])
levelplot(u[,,2])
levelplot(u[,,3])
levelplot(u[,,4])
levelplot(u[,,5])
levelplot(u[,,10])
levelplot(u[,,15])
levelplot(u[,,20])
levelplot(u[,,25])
levelplot(u[,,40])
levelplot(u[,,50])
#
# carrega matrix A
A <- matrix(ncol=nd,nrow=nd)
for(i in (1:(nd))){
  if(i==1){
    A[i,]<-c(1+1*fo,-fo,
             rep(0,nd-2))}
  if(i==nd){
    A[i,] <- c(rep(0,nd-2),
               -fo,1+1*fo)}
  if(i>1&&i<nd){
    A[i,]<-c(rep(0,i-2),
             -fo,1+2*fo,-fo,
             rep(0,nd-1-i))}
}
A
# loop no espaco x,y
v <- u
#j=2
for(tt in (0:(nt-2))){
  #for(j in (2:(nd-1))){
      for(i in (1:(nd))){
        if(i==1){
          b1 <- as.matrix(c(1-1*fo,fo))
          b2 <- as.matrix(v[,i:(i+1),tt+1])
          B <- t(b1)%*%t(b2)
        }
        if(i==nd){
          b1 <- as.matrix(c(fo,1-2*fo))
          b2 <- as.matrix(v[,(i-1):i,tt+1])
          B <- t(b1)%*%t(b2)
        }
        if(i>1&&i<nd){
          b1 <- as.matrix(c(fo,1-2*fo,fo))
          b2 <- as.matrix(v[,(i-1):(i+1),tt+1])
          B <- t(b1)%*%t(b2)
        }
      
      v[,i,tt+2] <- solve(A,t(B))
   #   }
  }
}

break
#img <- image_graph(600, 600, res = 96)
#contour(v[,,100],col = topo.colors(10))

levelplot(t(v[,,1]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*1))
levelplot(t(v[,,5]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*5))
levelplot(t(v[,,10]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*10))
levelplot(t(v[,,4]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*15))
levelplot(t(v[,,20]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*20))
levelplot(t(v[,,25]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*25))
levelplot(t(v[,,30]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*30))
levelplot(t(v[,,40]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*40))
levelplot(t(v[,,50]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*50))
levelplot(t(v[,,55]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*55))
levelplot(t(v[,,60]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*60))
levelplot(t(v[,,65]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*65))
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
