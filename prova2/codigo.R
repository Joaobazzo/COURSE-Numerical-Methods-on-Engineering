#
# ===============================================================
#
# Joao Pedro Bazzo
# 
# Equacoes de Saint-Venant unidimensional para um canal
# utilizando o metodo de Lax-Wendroff em duas etapas na marcha
# do tempo
#
# 12/12/2018
#
# ===============================================================
#
#
# a) funcao no contorno do lado esquerdo
#
# ===============================================================
#
f.H <- function(t){
  return(2+sin(0.00014*t+4.36))
}
#
# ===============================================================
#
# b) declaracao de variaveis
#
# ===============================================================
#
L = 1000   # m
tt = 1000  # s
# numero de valores
nx = 100   # 2*25
nt = 500   # 2*250
# valor infinitesimal
dx = L/nx
dt = tt/nt
# vetor com valores de x e t
vx = seq(0,L,length.out=nx)
vt = seq(0,tt,length.out = nt)
# courant
c = 1/f.H(0)+sqrt(9.81*f.H(0))
Cou = c*dt/dx
Cou
# H e Q
H = as.data.frame(matrix(data = NA,nrow = nx,ncol = nt))
Q = as.data.frame(matrix(data = NA,nrow = nx,ncol = nt))
#
# ===============================================================
#
# c) define valores iniciais
#
# ===============================================================
#
H[1,1:nt] = f.H(vt) # H(0,t)
Q[nx,1:nt] = 1      # Q(L,t) = 1
Q[1:nx,1] = 1       # Q(x,0) = 1      (hipotese)
H[2:nx,1] = H[1,1]  # H(x,0) = H(1,0) (hipotese)
#
# ===============================================================
#
# d) variacao no tempo
#
# ===============================================================
#
i = 1
j = 1
while(j <= (nt-1)){
  #
  # ===============================================================
  #
  # d.1) calculo dos invariantes de riemann (lado esquerdo)
  #
  # ===============================================================
  #
  psin <- dx/4
  diffe <- 1
  while(abs(diffe) >= 10^(-6)){
    
    h_psi <- (psin*H[2,j]+(dx/2-psin)*H[1,j])/(dx/2)
    c_delta <- psin/(dt/2)
    c_form <- sqrt(9.81*h_psi)
    diffe <- c_delta-c_form
    
    psi <- psin
    if(diffe < 0){psin <- psi - 0.1*diffe}
    else{psin <- psi +0.1*diffe}
    #print(diffe)
  }
  q_psi <- (psin*Q[2,j]+(dx/2-psin)*Q[1,j])/(dx/2)
  # invariante de riemann
  J <- q_psi/h_psi -2*sqrt(9.81*h_psi)
  #
  #
  # d.2) valor de Q[1,j+1] | Q(0,n+0.5)
  #
  Q[1,j+1] <- (J+2*sqrt(9.81*H[1,j+1]))*H[1,j+1]
  #
  # ===============================================================
  #
  # d.3) Passo 1 - Lax-Wendroff
  #
  # ===============================================================
  #
  for(i in (3:(nx-2))){
  # H
  H[i+1,j+1] <- (H[i+2,j]+H[i,j])*0.5-(0.5*dt/dx)*(Q[i+2,j]-Q[i,j])
  H[i-1,j+1] <- (H[i-2,j]+H[i,j])*0.5-(0.5*dt/dx)*(Q[i,j]-Q[i-2,j])
  
  # Q
  qt2 <- Q[i+2,j]^2/H[i+2,j]+9.81*0.5*H[i+2,j]^2
  qt0 <- Q[i-2,j]^2/H[i-2,j]+9.81*0.5*H[i-2,j]^2
  qt1 <- Q[i,j]^2/H[i,j]+9.81*0.5*H[i,j]^2
  
  Q[i+1,j+1] <- (Q[i+2,j]+Q[i,j])*0.5-(0.5*dt/dx)*(qt2-qt1)
  Q[i-1,j+1] <- (Q[i-2,j]+Q[i,j])*0.5-(0.5*dt/dx)*(qt1-qt0)
  }
  #
  # ===============================================================
  #
  # d.4) Fronteira direita H(nx,n+1/2)
  #
  # ===============================================================
  #
  # calculo do Invariante de Riemann
  #
  psin <- dx/4
  diffe <- 1
  while(abs(diffe) >= 10^(-6)){
    
    h_psi <- (psin*H[nx-1,j]+(dx/2-psin)*H[nx,j])/(dx/2)
    c_delta <- psin/(dt/2)
    c_form <- sqrt(9.81*h_psi)
    diffe <- c_delta-c_form
    
    psi <- psin
    if(diffe < 0){psin <- psi - 0.1*diffe}
    else{psin <- psi +0.1*diffe}
    #print(diffe)
  }
  h_psi <- (psin*H[nx-1,j]+(dx/2-psin)*H[nx,j])/(dx/2)
  q_psi <- (psin*Q[nx-1,j]+(dx/2-psin)*Q[nx,j])/(dx/2)
  J <- q_psi/h_psi +2*sqrt(9.81*h_psi)
  #
  #
  # valor de H[nx,j+1] | H(nx,n+0.5)
  #
  #
  diffe1 <- 1
  hnew <- 0.1
  #dh =1
  while(abs(diffe1) > 10^(-6)){
    Jd <-  Q[nx,j+1]/hnew +2*sqrt(9.81*hnew)
    diffe1 <- J-Jd 
    hold <- hnew
    if(diffe1 >0){hnew <- hold+0.1*diffe1}else{hnew <- hold-0.1*diffe1}
    #dh <- hnew-hold
    #print(diffe1)
    #break
  }
  H[nx,j+1] <- hnew
  #
  # ===============================================================
  #
  # d.5) Passo 2 - Lax-Wendroff
  #
  # ===============================================================
  #
  if(j==(nt-1)){break}else{
  for(i in (2:(nx-1))){
  # H
  H[i,j+2] <- H[i,j] - dt/dx*(Q[i+1,j+1]-Q[i-1,j+1])
  # Q
  q2 <- (Q[i+1,j+1]^2)/H[i+1,j+1]+9.81*0.5*H[i+1,j+1]^2
  q1 <- (Q[i-1,j+1]^2)/H[i-1,j+1]+9.81*0.5*H[i-1,j+1]^2
  Q[i,j+2] <- Q[i,j] - dt/dx*(q2-q1)
  }
  }
  j=j+1
  print(j)
}


