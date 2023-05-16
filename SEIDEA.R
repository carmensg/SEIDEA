library(readxl)
library("lpSolve")
library("MASS")
library(scales)

# A) Read data from file
# FORMAT: 
# DMU; INPUTS; OUTPUTS 
# data<-read_excel("VARIOABLES19.xlsx")

K = 2 # interval DEA: [a-,a+]
N = 12 # DMUs
M = 2 # inputs X
S = 4 # outputs Y

# B) Create directly the data: 
data = t(array(c(62899,77408,2591.8,2591.8,4973988,4973988,53,53,61.3,61.3,19273.02,19273.02,
                     187599,242707,3600.9,3600.9,7765647,7765647,18.5,18.5,21.4,21.4,85851.86,85851.86,
                     607778,791734,21318.8,21318.8,20717244,20717244,249.8,249.8,215.8,215.8,387720.85,387720.85,
                     393107,399659,9553.1,9553.1,15830794,15830794,146.5,146.5,126.6,126.6,221665.62,221665.62,
                     443019,467726,14843.4,14843.4,8439945,8439945,78.5,78.5,67.9,67.9,222775.81,222775.81,
                     616556,677732,11779.4,11779.4,17113029,17113029,68.3,68.3,74.7,74.7,256256.65,256256.65,
                     1080959,1170851,2808.9,2808.9,12219191,12219191,45.4,45.4,34.9,34.9,44339.56,44339.56,
                     794251,794251,29396.5,29396.5,6858772,6858772,257.9,257.9,250.9,250.9,504586.53,504586.53,
                     225166,225166,4662.9,4662.9,4878870,4878870,190.5,190.5,228.5,228.5,156907.85,156907.85,
                     210923,210923,3294.4,3294.4,5802640,5802640,149.7,149.7,179.6,179.6,127939.57,127939.57,
                     90188,90188,3172.1,3172.1,4241271,4241271,19.6,19.6,18.3,18.3,90495.4,90495.4,
                     48096,52667,2149.4,2149.4,3212455,3212455,6.8,6.8,12.7,12.7,65080.97,65080.97),
                   c(N,2*(M+S))))

data = data.frame(misdatos) 
DMUS = c("Attiki","Nisia Aigaiou, Kriti","Cataluña","Comunitat Valenciana","Illes Balears","Provence-Alpes-Côte d'Azur",
         "Jadranska Hrvatska","Veneto","Campania","Sicilia","Cyprus","Malta")

VARIABLES = c("BEDP","REC", "OVER","FEMP","MEMP","GHG")
colnames(misdatos) = c(rbind(paste0(VARIABLES,1),paste0(VARIABLES,2)))
rownames(misdatos) = DMUS

# IDENTIFY INPUTS/OUTPUTS :
idx = c(1,6) # 6th output is undesirable, so it is treated as an input. 
idy = c(2,3,4,5)

# Format of data: 
X = array(c(t(misdatos[1:N,c(1,2,11,12)])),c(K,M,N)) # 1:(2*M)

Y = array(c(t(misdatos[1:N,3:10])),c(K,S,N) ) #(1+2*M):(2*M+2*S)

# Change of units to get more homegenoeus data :  
X[,1,] = X[,1,] / 1e3
X[,2,] = X[,2,] / 1e3
Y[,c(2),] = Y[,c(2),]/1e3

# There are 
NV = N+(M+S)*5 # variables 
NC = 8*(M+S)+1 # constraints 

# VARIABLES:
EI = array(0,N)
SOL = array(0,c(N,NV))
lambda = array(0,c(N,N))
slx = array(0,c(K,M,N)) # same format than INPUTS 
sux = array(0,c(K,M,N)) # same format than INPUTS 
sly = array(0,c(K,S,N)) # same format than OUTPUTS 
suy = array(0,c(K,S,N)) # same format than OUTPUTS 
zx = array(0,c(M,N)) # INPUTS 
zy = array(0,c(S,N)) # OUTPUTS 

# CONSTANT definition :  maximum of the data 
# take L = R 
LR = round(max(X,Y))


# Direction of the constraints
dir  <- c(rep("=",(M+S)*2+1), rep("<=",(M+S)*6)) 
#dir  <- c(rep("<=",2*M), rep(">=",2*S),"=",rep("<=",(M+S)*6)) 


# MATRIX A
# Constraints matrix A, by blocks: (same matrix for each DMU --> outisde the loop)
A = array(0,c(NC,NV))

# INPUTS CONSTRAINTS
for(i in 1:M){
  A[i,1:N] = X[1,i,]
  A[M+i,1:N] = X[2,i,]
} 

A[1:M,c((N+1):(N+M),(N+3*M+1):(N+4*M))] = cbind(diag(1,M),diag(1,M))
A[(M+1):(2*M),c((N+M+1):(N+3*M))] = cbind(diag(1,M),diag(1,M))

# OUTPUTS  CONSTRAINTS
for(r in 1:S){
  A[2*M+r,1:N] = Y[1,r,]
  A[2*M+S+r,1:N] = Y[2,r,]
} 

A[(2*M+1):(2*M+S),
  c((N+4*M+1):(N+4*M+S),(N+4*M+3*S+1):(N+4*M+4*S))] = cbind(diag(-1,S),diag(-1,S))
A[(2*M+S+1):(2*M+2*S),
  c((N+4*M+S+1):(N+4*M+3*S))] = cbind(diag(-1,S),diag(-1,S))

# DEFINIR LA RESTRICCIONES DE LAMBDAS:
A[2*(M+S)+1,1:N] = rep(1,N)

# CONSTRAINTS slacks_low <= slacks_upper
A[(2*M+2*S+2):(3*M+2*S+1),(N+1):(N+2*M)] = cbind(diag(1,M),diag(-1,M))
A[(3*M+2*S+2):(4*M+2*S+1),(N+2*M+1):(N+4*M)] = cbind(diag(1,M),diag(-1,M))

A[(4*M+2*S+2):(4*M+3*S+1),(N+4*M+1):(N+4*M+2*S)] = cbind(diag(1,S),diag(-1,S))
A[(4*M+3*S+2):(4*M+4*S+1),(N+4*M+2*S+1):(N+4*M+4*S)] = cbind(diag(1,S),diag(-1,S))

# BOUNDARY CONSTRAINTS EQUIV. TO NON-LINEAR:
A[(4*M+4*S+2):(5*M+4*S+1),
  c((N+1):(N+M),(N+4*M+4*S+1):(N+5*M+4*S))] = cbind(diag(M),diag(-LR,M))
A[(5*M+4*S+2):(6*M+4*S+1),
  c((N+M+1):(N+2*M),(N+4*M+4*S+1):(N+5*M+4*S))] = cbind(diag(M),diag(-LR,M))

A[(6*M+4*S+2):(7*M+4*S+1),
  c((N+2*M+1):(N+3*M),(N+4*M+4*S+1):(N+5*M+4*S))] = cbind(diag(M),diag(LR,M))
A[(7*M+4*S+2):(8*M+4*S+1),
  c((N+3*M+1):(N+4*M),(N+4*M+4*S+1):(N+5*M+4*S))] = cbind(diag(M),diag(LR,M))

A[(8*M+4*S+2):(8*M+5*S+1),
  c((N+4*M+1):(N+4*M+S),(N+5*M+4*S+1):(N+5*M+5*S))] = cbind(diag(1,S),diag(-LR,S))
A[(8*M+5*S+2):(8*M+6*S+1),
  c((N+4*M+S+1):(N+4*M+2*S),(N+5*M+4*S+1):(N+5*M+5*S))] = cbind(diag(1,S),diag(-LR,S))

A[(8*M+6*S+2):(8*M+7*S+1),
  c((N+4*M+2*S+1):(N+4*M+3*S),(N+5*M+4*S+1):(N+5*M+5*S))] = cbind(diag(1,S),diag(LR,S))
A[(8*M+7*S+2):(8*M+8*S+1),
  c((N+4*M+3*S+1):(N+4*M+4*S),(N+5*M+4*S+1):(N+5*M+5*S))] = cbind(diag(1,S),diag(LR,S))


############################################################################
# We solve PEIMIL : INEFFICIENCY 

for (dmu in 1:N) {
  
    # C VECTOR : objective value  
    C <- rep(0,N)
    
    if(is.null(dim(X[,,dmu]))){
      C = c(C,rep(1/sum(X[,,dmu]),4))
    }else{
      C = c(C,rep(1/apply(X[,,dmu],2,sum),4))
    }
    if(is.null(dim(Y[,,dmu]))){
      C = c(C,rep(1/sum(Y[,,dmu]),4))
    }else{
      C = c(C,rep(1/apply(Y[,,dmu],2,sum),4))
    }
    
    C = c(C,rep(0,M+S))
  
    # Right hand side for the constraints
    B <- c(X[1,,dmu],X[2,,dmu],Y[1,,dmu],Y[2,,dmu],1,
           rep(0,4*M+2*S), rep(LR,2*M),rep(0,2*S),rep(LR,2*S))

    
    # SOLUTION: 
    # Note that every variable is assumed to be >= 0!
    optimum <-  lp(direction="max",
                   objective.in = C,
                   const.mat = A,
                   const.dir = dir,
                   const.rhs = B,
                   #int.vec = (N+4*M+4*M+1):NV
                   binary.vec = (N+4*M+4*S+1):NV
    )
    
    # IIDEA solution: 
    sol.PEIMIL <- optimum$solution
    
    EI[dmu] = optimum$objval
    SOL[dmu,] = sol.PEIMIL
    lambda[dmu,] = sol.PEIMIL[1:N]
    slx[,,dmu] = t(array(sol.PEIMIL[(N+1):(N+2*M)],c(M,2)))
    sux[,,dmu] = t(array(sol.PEIMIL[(N+2*M+1):(N+4*M)],c(M,2)))
    sly[,,dmu] = t(array(sol.PEIMIL[(N+4*M+1):(N+4*M+2*S)],c(S,2)))
    suy[,,dmu] = t(array(sol.PEIMIL[(N+4*M+2*S+1):(N+4*M+4*S)],c(S,2)))
    zx[,dmu] = sol.PEIMIL[(N+4*M+4*S+1):(N+5*M+4*S)]
    zy[,dmu] = sol.PEIMIL[(N+5*M+4*S+1):NV]
    
    
}

## TARGETS:
Targ.X = X
Targ.Y = Y 

for (dmu in 1:N){
  Targ.X[,,dmu] = X[,,dmu] - slx[,,dmu] - sux[2:1,,dmu]
  Targ.Y[,,dmu] = Y[,,dmu] + sly[,,dmu] + suy[2:1,,dmu]
}




############################################################################
## SUPER EFFICIY 

# VARIABLES:
SEI = array(NA,N)
SOL.SEI = array(NA,c(N,NV-1))
lambda.SEI = array(NA,c(N,N-1))
slx.SEI = array(NA,c(K,M,N)) # same format than INPUTS 
sux.SEI = array(NA,c(K,M,N)) # same format than INPUTS 
sly.SEI = array(NA,c(K,S,N)) # same format than OUTPUTS 
suy.SEI = array(NA,c(K,S,N)) # same format than OUTPUTS 
zx.SEI = array(NA,c(M,N)) # INPUTS 
zy.SEI= array(NA,c(S,N)) # OUTPUTS 

A.SEI = array(0,c(NC,NV))

# INPUTS CONSTRAINTS
for(i in 1:M){
  A.SEI[i,1:N] = X[1,i,]
  A.SEI[M+i,1:N] = X[2,i,]
} 

# A.SEI[1:M,c((N+1):(N+M),(N+3*M+1):(N+4*M))] = cbind(diag(1,M),diag(1,M))
# A.SEI[(M+1):(2*M),c((N+M+1):(N+3*M))] = cbind(diag(1,M),diag(1,M))
A.SEI[1:M,c((N+M+1):(N+3*M))] = cbind(diag(-1,M),diag(-1,M))
A.SEI[(M+1):(2*M),c((N+1):(N+M),(N+3*M+1):(N+4*M))] = cbind(diag(-1,M),diag(-1,M))

# OUTPUTS  CONSTRAINTS
for(r in 1:S){
  A.SEI[2*M+r,1:N] = Y[1,r,]
  A.SEI[2*M+S+r,1:N] = Y[2,r,]
} 

# A.SEI[(2*M+1):(2*M+S),c((N+4*M+1):(N+4*M+S),(N+4*M+3*S+1):(N+4*M+4*S))] = cbind(diag(-1,S),diag(-1,S))
# A.SEI[(2*M+S+1):(2*M+2*S),c((N+4*M+S+1):(N+4*M+3*S))] = cbind(diag(-1,S),diag(-1,S))
A.SEI[(2*M+1):(2*M+S),c((N+4*M+S+1):(N+4*M+3*S))] = cbind(diag(1,S),diag(1,S))
A.SEI[(2*M+S+1):(2*M+2*S),c((N+4*M+1):(N+4*M+S),(N+4*M+3*S+1):(N+4*M+4*S))] = cbind(diag(1,S),diag(1,S))

# LAMBDAS CONSTRAINTS:
A.SEI[2*(M+S)+1,1:N] = rep(1,N)

# CONSTRAINTS slacks_low <= slacks_upper
A.SEI[(2*M+2*S+2):(3*M+2*S+1),(N+1):(N+2*M)] = cbind(diag(1,M),diag(-1,M))
A.SEI[(3*M+2*S+2):(4*M+2*S+1),(N+2*M+1):(N+4*M)] = cbind(diag(1,M),diag(-1,M))

A.SEI[(4*M+2*S+2):(4*M+3*S+1),(N+4*M+1):(N+4*M+2*S)] = cbind(diag(1,S),diag(-1,S))
A.SEI[(4*M+3*S+2):(4*M+4*S+1),(N+4*M+2*S+1):(N+4*M+4*S)] = cbind(diag(1,S),diag(-1,S))

# BOUNDARY CONSTRAINTS EQUIV. TO NON LINEAR:
A.SEI[(4*M+4*S+2):(5*M+4*S+1),
  c((N+1):(N+M),(N+4*M+4*S+1):(N+5*M+4*S))] = cbind(diag(M),diag(-LR,M))
A.SEI[(5*M+4*S+2):(6*M+4*S+1),
  c((N+M+1):(N+2*M),(N+4*M+4*S+1):(N+5*M+4*S))] = cbind(diag(M),diag(-LR,M))

A.SEI[(6*M+4*S+2):(7*M+4*S+1),
  c((N+2*M+1):(N+3*M),(N+4*M+4*S+1):(N+5*M+4*S))] = cbind(diag(M),diag(LR,M))
A.SEI[(7*M+4*S+2):(8*M+4*S+1),
  c((N+3*M+1):(N+4*M),(N+4*M+4*S+1):(N+5*M+4*S))] = cbind(diag(M),diag(LR,M))

A.SEI[(8*M+4*S+2):(8*M+5*S+1),
  c((N+4*M+1):(N+4*M+S),(N+5*M+4*S+1):(N+5*M+5*S))] = cbind(diag(1,S),diag(-LR,S))
A.SEI[(8*M+5*S+2):(8*M+6*S+1),
  c((N+4*M+S+1):(N+4*M+2*S),(N+5*M+4*S+1):(N+5*M+5*S))] = cbind(diag(1,S),diag(-LR,S))

A.SEI[(8*M+6*S+2):(8*M+7*S+1),
  c((N+4*M+2*S+1):(N+4*M+3*S),(N+5*M+4*S+1):(N+5*M+5*S))] = cbind(diag(1,S),diag(LR,S))
A.SEI[(8*M+7*S+2):(8*M+8*S+1),
  c((N+4*M+3*S+1):(N+4*M+4*S),(N+5*M+4*S+1):(N+5*M+5*S))] = cbind(diag(1,S),diag(LR,S))


ind.dmu = which(EI == 0)

for (dmu in ind.dmu) { 
  
  C <- rep(0,N-1)
  
  if(is.null(dim(X[,,dmu]))){
    C = c(C,rep(1/sum(X[,,dmu]),4))
  }else{
    C = c(C,rep(1/apply(X[,,dmu],2,sum),4))
  }
  if(is.null(dim(Y[,,dmu]))){
    C = c(C,rep(1/sum(Y[,,dmu]),4),rep(0,M+S))
  }else{
    C = c(C,rep(1/apply(Y[,,dmu],2,sum),4),rep(0,M+S))
  }
  
  # Right hand side for the constraints
  B <-c(X[1,,dmu],X[2,,dmu],Y[1,,dmu],Y[2,,dmu],1,
         rep(0,4*M+2*S), rep(LR,2*M),rep(0,2*S),rep(LR,2*S))
  
  #dir  <- c(rep("=",(M+S)*2+1), rep("<=",(M+S)*6)) 
  dir  <- c(rep("<=",2*M), rep(">=",2*S),"=",rep("<=",(M+S)*6)) 
 
 # SOLUTION: 
  # Note that every variable is assumed to be >= 0!
  optimum <-  lp(direction="min",
                 objective.in = C, #round(C,5),
                 const.mat = A.SEI[,-dmu] , #round(A[,-dmu],5), #round(A.SEI[,-dmu],5),
                 const.dir = dir,
                 const.rhs = B, #round(B,5),
                 binary.vec = (N+4*M+4*S):(NV-1)
  )
  print(paste("DMU:",dmu,"status:",optimum$status))
  #  status = 2 = no feasible solution
  # IIDEA solution:
  SEI[dmu] = optimum$objval
  SOL.SEI[dmu,] = optimum$solution
  lambda.SEI[dmu,] = SOL.SEI[dmu,1:(N-1)]
  slx.SEI[,,dmu] = t(array(SOL.SEI[dmu,(N):(N-1+2*M)],c(M,2)))
  sux.SEI[,,dmu] = t(array(SOL.SEI[dmu,(N+2*M):(N-1+4*M)],c(M,2)))
  sly.SEI[,,dmu] = t(array(SOL.SEI[dmu,(N+4*M):(N-1+4*M+2*S)],c(S,2)))
  suy.SEI[,,dmu] = t(array(SOL.SEI[dmu,(N+4*M+2*S):(N-1+4*M+4*S)],c(S,2)))
  zx.SEI[,dmu] = SOL.SEI[dmu,(N+4*M+4*S):(N-1+5*M+4*S)]
  zy.SEI[,dmu] = SOL.SEI[dmu,(N+5*M+4*S):(NV-1)]
  
}



#########################
## SOME PLOTS 
#########################


library(ggplot2)
library(colorRamps)
library(RColorBrewer)
library(gridExtra)

library(pracma)


#colores = brewer.pal(N,name="Set3")
#colores = rainbow(N)
colores = array(0,N)
colores = c("lightpink1","maroon1","red2","darkorange",
                  "seagreen2","green","green4",
                  "skyblue1","cyan1","blue",
                  "purple1","slateblue1")

ind = which(EI == 0)
ind.ineff = which(EI > 0)

Y.scale = array(0,c(S,N))
TY.scale = array(0,c(S,N))
for(r in 1:S){
  Y.scale[r,] = 0.5*(Y[1,r,] + Y[2,r,]) / max(Y[,r,]) #OUTPUTS
  TY.scale[r,] = 0.5*(Targ.Y[1,r,] + Targ.Y[2,r,]) / max(Targ.Y[,r,]) #OUTPUTS
}

X.scale = array(0,c(M,N))
TX.scale = array(0,c(M,N))
for(i in 1:M){
  X.scale[i,] = 1 - (0.5*(X[1,c(2,1)[i],]+X[2,c(2,1)[i],])/max(X[,c(2,1)[i],])) # primero OUTPUT5
  X.scale[i,] = 1 - (0.5*(X[1,c(2,1)[i],]+X[2,c(2,1)[i],])/max(X[,c(2,1)[i],])) # después INPUT
  TX.scale[i,] = 1 - (0.5*(Targ.X[1,c(2,1)[i],]+Targ.X[2,c(2,1)[i],])/max(Targ.X[,c(2,1)[i],])) # primero OUTPUT5
  TX.scale[i,] = 1 - (0.5*(Targ.X[1,c(2,1)[i],]+Targ.X[2,c(2,1)[i],])/max(Targ.X[,c(2,1)[i],])) # después INPUT
}

leyenda = c()
xdmu = array(0,c(S+M+1,N))
xdmu.targ = array(0,c(S+M+1,N))

ydmu = array(0,c(S+M+1,N))
ydmu.targ = array(0,c(S+M+1,N))

angle = seq(pi/2,pi/2+2*pi,by=(2*pi)/(S+M))

for(j in 1:N){
  for(r in 1:S){
    xdmu[r,j] = (Y.scale[r,j])*cos(angle[r]) #/max(Y[1,r,]
    ydmu[r,j] = (Y.scale[r,j])*sin(angle[r])   #/max(Y[1,r,]
    xdmu.targ[r,j] = (TY.scale[r,j])*cos(angle[r]) #/max(Y[1,r,]
    ydmu.targ[r,j] = (TY.scale[r,j])*sin(angle[r])   #/max(Y[1,r,]
  }
  for(i in 1:M){
    xdmu[S+i,j] = (X.scale[i,j])*cos(angle[i+S]) #/max(Y[1,r,]
    ydmu[S+i,j] = (X.scale[i,j])*sin(angle[i+S])   #/max(Y[1,r,]
    xdmu.targ[S+i,j] = (TX.scale[i,j])*cos(angle[i+S]) #/max(Y[1,r,]
    ydmu.targ[S+i,j] = (TX.scale[i,j])*sin(angle[i+S])   #/max(Y[1,r,]
  }
  # CERRAR PENTAGRAMA: 
  xdmu[M+S+1,j] = xdmu[1,j]
  ydmu[M+S+1,j] = ydmu[1,j]
  xdmu.targ[M+S+1,j] = xdmu.targ[1,j]
  ydmu.targ[M+S+1,j] = ydmu.targ[1,j]
  
  leyenda = c(leyenda,paste0("DMU",j))
}

R = seq(0.1,round((max(xdmu,ydmu)-min(xdmu,ydmu))/2),length=N)

lineas = rep(2,N)
lineas[ind] = 1
width = rep(0.7,N)
width[ind] = 1.2
ptos = rep(16,N)
ptos[ind.ineff] = 17
pts.size = seq(0.5,2.,length=N)





HEXAGON = data.frame(x = Reshape(xdmu,84,1),
                     y = Reshape(ydmu,84,1),
                     Tx = Reshape(xdmu.targ,84,1),
                     Ty = Reshape(ydmu.targ,84,1),
                     dmu = rep(1:N,each = 7))


p = ggplot(HEXAGON, 
           aes(x,y, group = factor(dmu,levels = seq(1,N)),
               color = factor(dmu,levels = seq(1,N)),
               linetype = factor(dmu,levels = seq(1,N))
           )) 

for(j in 1:N){
  p <- p + geom_path(data = HEXAGON[which(HEXAGON$dmu==j),],
                     lwd=1.2) +
    geom_point(data = HEXAGON[which(HEXAGON$dmu==j),],
               cex= 3 )
  
  p <- p + geom_path(data = data.frame(x = R[j]*cos(angle),
                                       y = R[j]*sin(angle),
                                       dmu = j),
                     col="grey60",linetype = 3,lwd=0.5)
  
}

for(r in 1:(M+S)){
  p <- p + geom_path(data = data.frame(x = R*cos(angle[r]),
                                       y = R*sin(angle[r])),
                                       #dmu = 1:N),
                     col="grey50",linetype = 3,lwd=0.5)
}

p = p + geom_text(data = data.frame(label = c(paste("Output", seq(1:4)),"Output* 5", "Input"), 
                                    x = 1.1*cos(angle[-7]), 
                                    y = 1.1*sin(angle[-7])),
                  mapping = aes(x = x, y = y, label = label), 
                  col= c(rep("grey50",6)), #,"grey2"
                  angle = rep(0,6), size = 8) +
  labs(x="",y="") +
  #xlim(c(-1,2))+ylim(c(-2.5,1)) +
  scale_color_manual(name = "DMU",values = colores,
                     breaks = seq(1,N),
                     labels = DMUS
  ) +
  scale_linetype_manual(name = "DMU", 
                        breaks = c(1:N),
                        labels = DMUS,
                        values=lineas) +
  theme_bw() +
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         linestyle=guide_legend(nrow=3, byrow=TRUE)) +
  theme(legend.key.width = unit(2.5, "cm"),
        legend.text=element_text(size=18),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())




##################

pdf("HEXAGON.pdf",15,14)

p

dev.off()


##################

subset1 = sort(c(ind.ineff[1],which(lambda[ind.ineff[1],] >0)))

p1 = ggplot(HEXAGON[which(HEXAGON$dmu %in%subset1),], 
           aes(x,y, 
               group = factor(dmu,levels = subset1),
               color = factor(dmu,levels = subset1),
               linetype = factor(dmu,levels = subset1)
           )) 

for(j in subset1){
  p1 <- p1 + geom_path(data = HEXAGON[which(HEXAGON$dmu==j),],
                     lwd=1.2) +
    geom_point(data = HEXAGON[which(HEXAGON$dmu==j),],
               cex= 3 )+ 
    geom_path(data = data.frame(x = R[j]*cos(angle),
                                y = R[j]*sin(angle)),
                                #dmu = j),
              col="grey60",linetype = 3,lwd=0.5)
  
}

for(r in 1:(M+S)){
  p1 <- p1 + geom_path(data = data.frame(x = R*cos(angle[r]),
                                       y = R*sin(angle[r])),
                     #dmu = 1:N),
                     col="grey50",linetype = 3,lwd=0.5)
}


p1 = p1 + geom_path(data = HEXAGON[which(HEXAGON$dmu==ind.ineff[1]),],
                    aes(Tx,Ty), col = "black",lwd=1,linetype=5) +
  geom_point(data = HEXAGON[which(HEXAGON$dmu==ind.ineff[1]),],
             aes(Tx,Ty), col = "black",pch=22,cex=3) +
  geom_text(data = data.frame(label = c(paste("Out", seq(1:4)),"Out* 5", "In"), 
                              x = 1.1*cos(angle[-7]), 
                              y = 1.1*sin(angle[-7])),
            mapping = aes(x = x, y = y, label = label), 
            col= c(rep("grey50",6)), #,"grey2"
            angle = 0, size = 4) +
  labs(x="",y="") +
  #xlim(c(-1,2))+ylim(c(-2.5,1)) +
  scale_color_manual(name = "DMU",
                     values = colores[subset1],
                     breaks = subset1,
                     labels = DMUS[subset1]
  ) +
  scale_linetype_manual(name = "DMU", 
                        values=lineas[subset1],
                        breaks = subset1,
                        labels = DMUS[subset1]) +
  theme_bw() +
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         linestyle=guide_legend(nrow=3, byrow=TRUE)) +
  theme(legend.key.width = unit(2.5, "cm"),
        legend.text=element_text(size=18),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



##################
subset2 = sort(c(ind.ineff[2],which(lambda[ind.ineff[2],] >0)))

p2 = ggplot(HEXAGON[which(HEXAGON$dmu %in%subset2),], 
            aes(x,y, 
                group = factor(dmu,levels = subset2),
                color = factor(dmu,levels = subset2),
                linetype = factor(dmu,levels = subset2)
            )) 

for(j in subset2){
  p2 <- p2 + geom_path(data = HEXAGON[which(HEXAGON$dmu==j),],
                       lwd=1.2) +
    geom_point(data = HEXAGON[which(HEXAGON$dmu==j),],
               cex= 3 )+ 
    geom_path(data = data.frame(x = R[j]*cos(angle),
                                y = R[j]*sin(angle)),
              #dmu = j),
              col="grey60",linetype = 3,lwd=0.5)
  
}

for(r in 1:(M+S)){
  p2 <- p2 + geom_path(data = data.frame(x = R*cos(angle[r]),
                                         y = R*sin(angle[r])),
                       #dmu = 1:N),
                       col="grey50",linetype = 3,lwd=0.5)
}


p2 = p2 + geom_path(data = HEXAGON[which(HEXAGON$dmu==ind.ineff[2]),],
                    aes(Tx,Ty), col = "black",lwd=1,linetype=5) +
  geom_point(data = HEXAGON[which(HEXAGON$dmu==ind.ineff[2]),],
             aes(Tx,Ty), col = "black",pch=22,cex=3) +
  geom_text(data = data.frame(label = c(paste("Out", seq(1:4)),"Out* 5", "In"), 
                              x = 1.1*cos(angle[-7]), 
                              y = 1.1*sin(angle[-7])),
            mapping = aes(x = x, y = y, label = label), 
            col= c(rep("grey50",6)), #,"grey2"
            angle = 0, size = 4) +
  labs(x="",y="") +
  #xlim(c(-1,2))+ylim(c(-2.5,1)) +
  scale_color_manual(name = "DMU",
                     values = colores[subset2],
                     breaks = subset2,
                     labels = DMUS[subset2]
  ) +
  scale_linetype_manual(name = "DMU", 
                        values=lineas[subset2],
                        breaks = subset2,
                        labels = DMUS[subset2]) +
  theme_bw() +
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         linestyle=guide_legend(nrow=3, byrow=TRUE)) +
  theme(legend.key.width = unit(2.5, "cm"),
        legend.text=element_text(size=18),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



##################
subset3 = sort(c(ind.ineff[3],which(lambda[ind.ineff[3],] >0)))

p3 = ggplot(HEXAGON[which(HEXAGON$dmu %in%subset3),], 
            aes(x,y, 
                group = factor(dmu,levels = subset3),
                color = factor(dmu,levels = subset3),
                linetype = factor(dmu,levels = subset3)
            )) 

for(j in subset3){
  p3 <- p3 + geom_path(data = HEXAGON[which(HEXAGON$dmu==j),],
                       lwd=1.2) +
    geom_point(data = HEXAGON[which(HEXAGON$dmu==j),],
               cex= 3 )+ 
    geom_path(data = data.frame(x = R[j]*cos(angle),
                                y = R[j]*sin(angle)),
              #dmu = j),
              col="grey60",linetype = 3,lwd=0.5)
  
}

for(r in 1:(M+S)){
  p3 <- p3 + geom_path(data = data.frame(x = R*cos(angle[r]),
                                         y = R*sin(angle[r])),
                       #dmu = 1:N),
                       col="grey50",linetype = 3,lwd=0.5)
}


p3 = p3 + geom_path(data = HEXAGON[which(HEXAGON$dmu==ind.ineff[3]),],
                    aes(Tx,Ty), col = "black",lwd=1,linetype=5) +
  geom_point(data = HEXAGON[which(HEXAGON$dmu==ind.ineff[3]),],
             aes(Tx,Ty), col = "black",pch=22,cex=3) +
  geom_text(data = data.frame(label = c(paste("Out", seq(1:4)),"Out* 5", "In"), 
                              x = 1.1*cos(angle[-7]), 
                              y = 1.1*sin(angle[-7])),
            mapping = aes(x = x, y = y, label = label), 
            col= c(rep("grey50",6)), #,"grey2"
            angle = 0, size = 4) +
  labs(x="",y="") +
  #xlim(c(-1,2))+ylim(c(-2.5,1)) +
  scale_color_manual(name = "DMU",
                     values = colores[subset3],
                     breaks = subset3,
                     labels = DMUS[subset3]
  ) +
  scale_linetype_manual(name = "DMU", 
                        values=lineas[subset3],
                        breaks = subset3,
                        labels = DMUS[subset3]) +
  theme_bw() +
  guides(color=guide_legend(nrow=3, byrow=TRUE),
         linestyle=guide_legend(nrow=3, byrow=TRUE)) +
  theme(legend.key.width = unit(2.5, "cm"),
        legend.text=element_text(size=18),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())





pdf("HEXAGON_bis.pdf",17,20)

grid.arrange(p,p1,p2,p3,
             #ncol = 2
             layout_matrix = rbind(c(1,1,1),
                                   c(1,1,1),
                                   c(1,1,1),
                                   c(2,3,4)))


dev.off()



library(ggbreak) 

ind = which(EI == 0)
ind.ineff = which(EI > 0)

INEFF = data.frame(dmu = rep(ind.ineff,each = 6),
                   vars = rep(c("INPUT", paste("OUTPUT",1:5)),3),
                   low = Reshape(rbind(X[1,1,ind.ineff]/max(X[,1,]),
                                       Y[1,,ind.ineff]/max(Y),
                                       X[1,2,ind.ineff]/max(X[,2,])),18,1),
                   upp = Reshape(rbind(X[2,1,ind.ineff]/max(X[,1,]),
                                       Y[2,,ind.ineff]/max(Y),
                                       X[2,2,ind.ineff]/max(X[,2,])),18,1),
                   low.targ = Reshape(rbind(Targ.X[1,1,ind.ineff]/max(X[,1,]),
                                            Targ.Y[1,,ind.ineff]/max(Y),
                                            Targ.X[1,2,ind.ineff]/max(X[,2,])),18,1),
                   upp.targ = Reshape(rbind(X[2,1,ind.ineff]/max(X[,1,]),
                                            Targ.Y[2,,ind.ineff]/max(Y),
                                            Targ.X[2,2,ind.ineff]/max(X[,2,])),18,1)
                   )


INEFF = data.frame(dmu = rep(ind.ineff,each = 6),
                   vars = rep(c("INPUT", paste("OUTPUT",1:5)),3),
                   low = Reshape(rbind(X[1,1,ind.ineff],
                                       Y[1,,ind.ineff],
                                       X[1,2,ind.ineff]),18,1),
                   upp = Reshape(rbind(X[2,1,ind.ineff],
                                       Y[2,,ind.ineff],
                                       X[2,2,ind.ineff]),18,1),
                   low.targ = Reshape(rbind(Targ.X[1,1,ind.ineff],
                                            Targ.Y[1,,ind.ineff],
                                            Targ.X[1,2,ind.ineff]),18,1),
                   upp.targ = Reshape(rbind(X[2,1,ind.ineff],
                                            Targ.Y[2,,ind.ineff],
                                            Targ.X[2,2,ind.ineff]),18,1)
)

INEFF$midp = (INEFF$low + INEFF$upp) / 2
INEFF$sig1 = round(abs(INEFF$midp - INEFF$low),5)
INEFF$sig2 = round(abs(INEFF$midp - INEFF$upp),5)

INEFF$midp.targ = (INEFF$low.targ + INEFF$upp.targ) / 2
INEFF$sig1.targ = round(abs(INEFF$midp.targ - INEFF$low.targ),5)
INEFF$sig2.targ = round(abs(INEFF$midp.targ - INEFF$upp.targ),5)

# INEFF$midp = apply(round(INEFF[,c(3,4)],5),  1, function(x) if(x[1] < x[2] ) return((x[1]+x[2])/2)  else return(0))
# INEFF$midp.targ = apply(round(INEFF[,c(5,6)],5),  1, function(x) if(x[1] < x[2] ) return((x[1]+x[2])/2)  else return(0))

p.dmu2 = ggplot(data = INEFF[which(INEFF$dmu==ind.ineff[1]),])  +
  geom_errorbar(aes(x = vars, y = 1,
                    ymin=1-sig1/midp,
                    ymax=1+sig2/midp),width=.1,lwd=1) +
  geom_errorbar(aes(x = vars, y = midp.targ/midp,
                    ymin=midp.targ/midp-sig1.targ/midp,
                    ymax=midp.targ/midp+sig2.targ/midp),
                col= "red2",width=.1,linetype=1,lwd=0.75) +
  geom_point(aes(x = vars, y = 1),cex=2) + 
  geom_point(aes(x = vars, y = midp.targ/midp),col="red2",cex=4,pch=22) + 
  # geom_text(data = data.frame(label = paste("DMU", ind.ineff[1]), 
  #                             vars = "INPUT", 
  #                             y = 0.45),
  #           mapping = aes(x = vars, y = y, label = label), 
  #           angle = 0, size = 7) +
  #scale_y_continuous(limits = c(-0.1,0.3)) +
  labs(x="", y = "scaled values (x100%)",title = paste("DMU",ind.ineff[1])) + 
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title.y  = element_text(size=18,angle=90),
        axis.title.x = element_text(size = 18),
        title=element_text(size=16)
  )
  
  
p.dmu6 = ggplot(data = INEFF[which(INEFF$dmu==ind.ineff[2]),])  +
  geom_errorbar(aes(x = vars, y = 1,
                    ymin=1-sig1/midp,
                    ymax=1+sig2/midp),width=.1,lwd=1) +
  geom_errorbar(aes(x = vars, y = midp.targ/midp,
                    ymin=midp.targ/midp-sig1.targ/midp,
                    ymax=midp.targ/midp+sig2.targ/midp),
                col= "red2",width=.1,linetype=1,lwd=0.75) +
  geom_point(aes(x = vars, y = 1),cex=2) + 
  geom_point(aes(x = vars, y = midp.targ/midp),col="red2",cex=4,pch=22) + 
  # geom_text(data = data.frame(label = paste("DMU", ind.ineff[1]), 
  #                             vars = "INPUT", 
  #                             y = 0.45),
  #           mapping = aes(x = vars, y = y, label = label), 
  #           angle = 0, size = 7) +
  #scale_y_continuous(limits = c(-0.1,0.3)) +
  labs(x="", y = "scaled values (x100%)",title = paste("DMU",ind.ineff[2])) + 
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title.y  = element_text(size=18,angle=90),
        axis.title.x = element_text(size = 18),
        title=element_text(size=16)
  )


p.dmu11 = ggplot(data = INEFF[which(INEFF$dmu==ind.ineff[3]),])  +
  geom_errorbar(aes(x = vars, y = 1,
                    ymin=1-sig1/midp,
                    ymax=1+sig2/midp),width=.1,lwd=1) +
  geom_errorbar(aes(x = vars, y = midp.targ/midp,
                    ymin=midp.targ/midp-sig1.targ/midp,
                    ymax=midp.targ/midp+sig2.targ/midp),
                col= "red2",width=.1,linetype=1,lwd=0.75) +
  geom_point(aes(x = vars, y = 1),cex=2) + 
  geom_point(aes(x = vars, y = midp.targ/midp),col="red2",cex=4,pch=22) + 
  # geom_text(data = data.frame(label = paste("DMU", ind.ineff[1]), 
  #                             vars = "INPUT", 
  #                             y = 0.45),
  #           mapping = aes(x = vars, y = y, label = label), 
  #           angle = 0, size = 7) +
  #scale_y_continuous(limits = c(-0.1,0.3)) +
  scale_y_break(c(1.5, 2.2)) + 
  labs(x="", y = "scaled values (x100%)",title = paste("DMU",ind.ineff[3])) + 
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title.y  = element_text(size=18,angle=90),
        axis.title.x = element_text(size = 18),
        title=element_text(size=16)
  )
  
  
pdf("INEFF_DMUs.pdf",10,15)

grid.arrange(p.dmu2,p.dmu6,p.dmu11,ncol=1)

dev.off()










