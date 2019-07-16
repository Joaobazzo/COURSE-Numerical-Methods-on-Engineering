#
#
# exercicio 2 - sol numerica
#
#
# 1) valores
dx = 1/100
dt = 1/10000
D = 10
L = 1
Fo = D*dt/(dx^2)
Fo
#
# 2) matriz u
nx = round(1/dx)
nt = round(1/dt)
u <- matrix(data = NA,nrow = nx,ncol = nt)
#
# condicoes iniciais
xi <- seq(0,1,length.out=nx)
for(n in (1:nx)){
  u[n,1] <- 2*xi[n]*(1-xi[n])
}
#
# condicoes de contorno
#
u[1,] <- 0 # u(0,t) = 0
u[nx,] <- 0 # u(L,t) = 0
#
#
# carrega matriz
#
for(t in (1:(nt-1))){
  for(i in (2:(nx-1))){
    if(t==1){
      u[i,t+1] = (1-2*Fo)*u[i,t]/(1+2*Fo)+2*Fo*(u[i+1,t]+u[i-1,t])/(1+2*Fo)
      }
    if(t>1){
      aux = (1-2*Fo)*u[i,t-1] + 2*Fo*(u[i+1,t]+u[i-1,t])
      u[i,t+1]=aux/(1+2*Fo)
      }
  }
}
# plota imagem
#
#
#
png(filename = "E:/Documents/CICLO/Mestrado/materias/met_mat_eng/prova1/ex2-3.png",
    width = 750,height = 650,units = "px",pointsize = 20)
par(mfrow=c(4,2),mar=c(2, 4, 2, 2))
plot((1:nx)/nx,u[,1],type="l",ylim=c(0,0.51),xaxs="i",yaxs="i",ylab="u(x,t=1)",xlab="x")
mtext(paste0("t=",1,", Fo=",Fo))
plot((1:nx)/nx,u[,9],type="l",ylim=c(0,0.51),xaxs="i",yaxs="i",ylab="u(x,t=9)",xlab="x")
mtext(paste0("t=",9,", Fo=",Fo))
plot((1:nx)/nx,u[,15],type="l",ylim=c(0,0.51),xaxs="i",yaxs="i",ylab="u(x,t=15)",xlab="x")
mtext(paste0("t=",15,", Fo=",Fo))
plot((1:nx)/nx,u[,17],type="l",ylim=c(0,0.51),xaxs="i",yaxs="i",ylab="u(x,t=17)",xlab="x")
mtext(paste0("t=",17,", Fo=",Fo))
plot((1:nx)/nx,u[,19],type="l",ylim=c(0,0.51),xaxs="i",yaxs="i",ylab="u(x,t=19)",xlab="x")
mtext(paste0("t=",19,", Fo=",Fo))
plot((1:nx)/nx,u[,45],type="l",ylim=c(0,0.51),xaxs="i",yaxs="i",ylab="u(x,t=45)",xlab="x")
mtext(paste0("t=",45,", Fo=",Fo))
plot((1:nx)/nx,u[,65],type="l",ylim=c(0,0.51),xaxs="i",yaxs="i",ylab="u(x,t=65)",xlab="x")
mtext(paste0("t=",65,", Fo=",Fo))
plot((1:nx)/nx,u[,90],type="l",ylim=c(0,0.51),xaxs="i",yaxs="i",ylab="u(x,t=90)",xlab="x")
mtext(paste0("t=",90,", Fo=",Fo))
dev.off()


