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
D = 2
L = 1
dt = 0.01
delta = 0.01
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
  u[,j,1] <- analitica(D,L,uo,xi[1:nd],yi[j],t = 0)
}
# plot condicoes iniciais
levelplot(u[,,1])
#plot3d(x=0:49,y=0:49,z=u[,,1],type = "l")
#
# condicoes de contorno u(L,y,t) = u(-L,y,t) = u(x,-L,t)= u(x,L,t)
for(tt in (1:nt)){
  for(i in (1:nd)){
      u[i,1:nd,tt]  <- analitica(D,L,uo,x = xi[i], y=yi[1:nd],t=tt-1)
  }

}
#plot(u[1,,35])
#
break

v <- u
levelplot(t(v[,,1]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*1))
levelplot(t(v[,,2]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*5))
levelplot(t(v[,,3]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*10))
levelplot(t(v[,,4]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*15))
levelplot(t(v[,,5]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*20))
levelplot(t(v[,,6]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*25))
levelplot(t(v[,,7]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*30))
levelplot(t(v[,,8]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*40))
levelplot(t(v[,,9]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*50))
levelplot(t(v[,,10]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*55))
levelplot(t(v[,,12]),xlab="x",ylab="y",at=seq(0,1,by=0.01),col.regions=heat.colors(100)[100:1],main=paste0("dt = ",dt*60))
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
break
par(mfrow=c(2,4))
levelplot(t(v[,,1]))
levelplot(t(v[,,2]))
levelplot(t(v[,,3]))
levelplot(t(v[,,4]))
levelplot(t(v[,,5]))
levelplot(t(v[,,6]))
levelplot(t(v[,,7]))
levelplot(t(v[,,8]))
levelplot(t(v[,,9]))

