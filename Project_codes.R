library(mvtnorm)
library(coda)
library(tmvtnorm)
require(ggmap)

register_google(key="AIzaSyBT-E0y_XjuxyIiwHp5rX3u_amsZiEYJ2U",write = TRUE)

#my location: 50.443492, -104.510223 (Cosco)  :lat/long
#first point: Co-op Refinery: 50.484095, -104.579288:   lat/long
#second point: Kodiak LED Lighting: 50.5180655, -104.556831: lat/long
#third point: Cornwall: 50.450717, -104.610967: lat/long
#my pos
pos=c(-104.510223, 50.443492)
landmarks=data.frame(lon=c(-104.579288,-104.556831,-104.610967),
                     lat=c(50.484095,50.5180655,50.450717))

file <- system.file("map.png",package = 'imager')
im <- load.image(file)

map <- get_map(pos,zoom=11,maptype = "road")
ggmap(map)
B <- function(P1,P2){
  
  (360+atan((P2[1]-P1[1])/(P2[2]-P1[2]))*180/pi) %% 360
}


alpha <- 312 #bearing between your pos and first point
beta <- 338 #bearing between your pos and second point
gamma <- 276 # bearing between your pos and third point
#lambda follow uniform (-105,-104)
#phi follow uniform (50,51)
#sigma follow exponential(0.05)


d <- seq(0,0.2,0.0001) # Length of the line
line1 <- data.frame(lon=landmarks[1,1] + d*sin(alpha*pi/180+pi),
                    lat=landmarks[1,2] + d*cos(alpha*pi/180+pi))
line2 <- data.frame(lon=landmarks[2,1] + d*sin(beta*pi/180+pi),
                    lat=landmarks[2,2] + d*cos(beta*pi/180+pi))
line3 <- data.frame(lon=landmarks[3,1] + d*sin(gamma*pi/180+pi),
                    lat=landmarks[3,2] + d*cos(gamma*pi/180+pi))

map <- get_map(pos,zoom=11,maptype="watercolor")
mapPlot <- ggmap(map)+
  geom_point(aes(x = lon, y = lat), size = 1, data = landmarks, alpha = .5) +
  geom_line(aes(x=lon,y=lat),data=line1) +
  geom_line(aes(x=lon,y=lat),data=line2) +
  geom_line(aes(x=lon,y=lat),data=line3)
mapPlot


#params=(lambda,phi,sigma)
loglikelihood <- function(data,params){
  dnorm(alpha, B(data[1,],params[1:2])[1,1],params[3],log = T)+
    dnorm(beta, B(data[2,],params[1:2])[1,1],params[3],log=T)+
    dnorm(gamma, B(data[3,],params[1:2])[1,1],params[3],log=T)
}

logprior <- function(params){
  if (params[3]<=0)
    return (0)
  dexp(params[3],20,log = T)
}


prop.cov <- c(4e-6,4e-6,2e-3)*diag(3)
draws <- array(0,dim=c(4000,3,4))
mu <- c(-104.5,50.5,0.05)
draws[4000,,1] <- rtmvnorm(1,mu,prop.cov,c(-Inf,-Inf,0),c(Inf,Inf,Inf),algorithm = "gibbs")
draws[4000,,2] <- rtmvnorm(1,mu,prop.cov,c(-Inf,-Inf,0),c(Inf,Inf,Inf),algorithm = "gibbs")
draws[4000,,3] <- rtmvnorm(1,mu,prop.cov,c(-Inf,-Inf,0),c(Inf,Inf,Inf),algorithm = "gibbs")
draws[4000,,4] <- rtmvnorm(1,mu,prop.cov,c(-Inf,-Inf,0),c(Inf,Inf,Inf),algorithm = "gibbs")


loop <- 1
converged <- FALSE
count <- 0 #we get sample from the second times that chains converge
while(count<2){
  draws[1,,] <- draws[4000,,]
  accepted <- c(1,1,1,1)
  for (step in 2:4000){
    for (chain in 1:4){
      proposed <- c(-104.5,50.5,-1)
      while(proposed[3]<=0){
        proposed <- rmvnorm(1,draws[step-1,,chain],prop.cov)
      }
      
      r <- loglikelihood(landmarks,proposed)+logprior(proposed) + 
        dtmvnorm(draws[step-1,,chain],as.vector(proposed),prop.cov,c(-Inf,-Inf,0),c(Inf,Inf,Inf),log=T)-
        loglikelihood(landmarks,draws[step-1,,chain])-logprior(draws[step-1,,chain])-
        dtmvnorm(as.vector(proposed),draws[step-1,,chain],prop.cov,c(-Inf,-Inf,0),c(Inf,Inf,Inf),log=T)
      
      u <- runif(1)
      if (log(u) < min(0,r)){
        draws[step,,chain] <- proposed
        accepted[chain] <- accepted[chain]+1
      } else {
        draws[step,,chain] <- draws[step-1,,chain]
      }
    }
  }
  chainlist <- mcmc.list(Chain1=mcmc(draws[,,1]),
                         Chain2=mcmc(draws[,,2]),
                         Chain3=mcmc(draws[,,3]),
                         Chain4=mcmc(draws[,,4]))
  converged <- all((gelman.diag(chainlist)$psrf[,2])<1.1)
  if(converged)
    count <- count+1
  cat("ItERATION",loop,"\n")
  cat("Upper Limits of PSRF:",gelman.diag(chainlist)$psrf[,2],"\n")
  cat("Acceptance rate:",accepted/40,"\n")
  loop <- loop+1
}

samples <- rbind(draws[,,1],draws[,,2],draws[,,3],draws[,,4])
dim(samples)

mean(samples[,1]) #lambda: -104.5274
mean(samples[,2]) #phi: 50.44016
mean(samples[,3]) #sigma: 0.68678




################################################################################################################################


prop.cov1 <- c(1e-8,1e-8,1e-6)*diag(3)
draws1 <- array(0,dim=c(4000,3,4))
mu <- c(-104.2,50.2,0.055)
draws1[4000,,1] <- rtmvnorm(1,mu,prop.cov1,c(-Inf,-Inf,0),c(Inf,Inf,Inf),algorithm = "gibbs")
draws1[4000,,2] <- rtmvnorm(1,mu,prop.cov1,c(-Inf,-Inf,0),c(Inf,Inf,Inf),algorithm = "gibbs")
draws1[4000,,3] <- rtmvnorm(1,mu,prop.cov1,c(-Inf,-Inf,0),c(Inf,Inf,Inf),algorithm = "gibbs")
draws1[4000,,4] <- rtmvnorm(1,mu,prop.cov1,c(-Inf,-Inf,0),c(Inf,Inf,Inf),algorithm = "gibbs")


loop1 <- 1 #count number of iterations
converged1 <- FALSE
while(!converged1){
  draws1[1,,] <- draws1[4000,,]
  accepted1 <- c(1,1,1,1)
  for (step in 2:4000){
    for (chain in 1:4){
      proposed <- c(-104.5,50.5,-1)
      while(proposed[3]<=0){
        proposed <- rmvnorm(1,draws1[step-1,,chain],prop.cov1)
      }
      
      r <- loglikelihood(landmarks,proposed)+new_logprior(proposed) + 
        dtmvnorm(draws1[step-1,,chain],as.vector(proposed),prop.cov1,c(-Inf,-Inf,0),c(Inf,Inf,Inf),log=T)-
        loglikelihood(landmarks,draws1[step-1,,chain])-new_logprior(draws1[step-1,,chain])-
        dtmvnorm(as.vector(proposed),draws1[step-1,,chain],prop.cov1,c(-Inf,-Inf,0),c(Inf,Inf,Inf),log=T)
      
      u <- runif(1)
      if (log(u) < min(0,r)){
        draws1[step,,chain] <- proposed
        accepted1[chain] <- accepted1[chain]+1
      } else {
        draws1[step,,chain] <- draws1[step-1,,chain]
      }
    }
  }
  chainlist <- mcmc.list(Chain1=mcmc(draws1[,,1]),
                         Chain2=mcmc(draws1[,,2]),
                         Chain3=mcmc(draws1[,,3]),
                         Chain4=mcmc(draws1[,,4]))
  converged1 <- all((gelman.diag(chainlist)$psrf[,2])<1.1)
  cat("ItERATION",loop1,"\n")
  cat("Upper Limits of PSRF:",gelman.diag(chainlist)$psrf[,2],"\n")
  cat("Acceptance rate:",accepted1/40,"\n")
  loop1 <- loop1+1
}

new_sample <- draws1
