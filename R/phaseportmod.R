# phaseportmod produces a phase portrait based on the best models (differential equations)
# with modeled trajectories for specified entities (e.g. countries)
# for x and y, the phase portrait is the visual simultaneous solution of the two differential equations

# to call: phaseportmod(datap, 30, datap$logGDP, datap$EmanzV, seq(0, 1, by = 0.1), seq(0, 1, by = 0.1), param <- c(0.0012, 0.0071), f <- function(t,Y=c()) rbind(0.0012/Y[1]^2, + 0.0071*Y[1]^3), 1, 2, 4, 5, 7, 9)
# yearnr can be the number of the observations in the data, but in can be also go beyond, to make
# predictions of trajectories of the specified entities

phaseportmod <- function(dataset, yearnr, xv, yv, rangeX, rangeY, param, f, entidx1, entidx2, 
                      entidx3, entidx4, entidx5, entidx6) 
{  
  procdata <- preprocess_data(2, xv, yv)
  xwide <- procdata$xwide
  ywide <- procdata$ywide
  #   
  #   # boundaries of the system, here under assumption that both variables are scaled 0-1
  #   # the visualisation is optimized for this scaling, in case the user would like to use
  #   # the original scale, he/she needs to adjust the code
  xmin=0 
  xmax=1 
  ymin=0 
  ymax=1  
  rgrid = 21
  
  y1 <- linspace(xmin, xmax, rgrid)
  y2 <- linspace(ymin, ymax, rgrid)
  
  # scaling arrows of the vectorfield
  qscale <- range(y1)/range(y2)
  mgrid <- meshgrid(y1, y2)
  x <- mgrid[[1]] 
  y <- mgrid[[2]] 
  
  u <- matrix(0, rgrid, rgrid) 
  v <- matrix(0, rgrid, rgrid)
  
  # compute derivates at each point (every possible intial conditions)
  t = 0
  for (i in 1:length(x))
  {
    Yprime <- f(t, rbind(x[i], y[i]))
    u[i] <- Yprime[1]
    v[i] <- Yprime[2]
  }
  
  # phase portrait with highlighted model trajectories
  dev.set(1)
  postscript("ppmod.eps", horizontal=FALSE, width=5, height=5,
             onefile=FALSE, paper="special", family="ComputerModern")
  
  # setting a plot 
  plot(rangeX, rangeY, col="white", xlab = "X-Variable", ylab = "Y-Variable")
  grid(col="white")
  # velocity plot
  quiver(x, y, u, v, scale=0.5, col = 'gray') 
  
  legend("topleft", bg="white", legend=c('Entity1  ', 'Entity2  ', 'Entity3  ', 'Entity4 ', 'Entity5 ', 'Entity6 '), 
         lwd = 2, col=c('blue', 'darkgreen', 'red', 'cyan', 'magenta', 'black'))
  
  # defining initial conditions/values, 
  # example, if the data has missings and different first time of observation
#   inits1 <- rbind(xwide[entidx1,22], ywide[entidx1,22])
#   inits2 <- rbind(xwide[entidx2,4], ywide[entidx2,4])
#   inits3 <- rbind(xwide[entidx3,15], ywide[entidx3,15])
#   inits4 <- rbind(xwide[entidx4,10], ywide[entidx4,10])
#   inits5 <- rbind(xwide[entidx5,16], ywide[entidx5,16])
#   inits6 <- rbind(xwide[entidx6,1], ywide[entidx6,1]) 
  
  # defining initial conditions/values (default)
  inits1 <- rbind(xwide[entidx1,1], ywide[entidx1,1])
  inits2 <- rbind(xwide[entidx2,1], ywide[entidx2,1])
  inits3 <- rbind(xwide[entidx3,1], ywide[entidx3,1])
  inits4 <- rbind(xwide[entidx4,1], ywide[entidx4,1])
  inits5 <- rbind(xwide[entidx5,1], ywide[entidx5,1])
  inits6 <- rbind(xwide[entidx6,1], ywide[entidx6,1]) 
  
  # defining function for ode (integration of ordinary differential equations)
  f1 <- function(t, state, param)
  { 
    ft <- f(t, rbind(state[1], state[2]))  
    with(as.list(ft), 
    { 
      dX <- ft[1]
      dY <- ft[2]
        list(c(dX, dY))
    })
  } 
  
  # solving/integrating ordinary differential equtions with defined initial values
  times <- seq(0, yearnr, by = 1) # defining the timespan and changes per year
  sol1 <- ode(y = inits1, times = times, func = f1, parms = param, method = "ode45")
  sol2 <- ode(y = inits2, times = times, func = f1, parms = param, method = "ode45")
  sol3 <- ode(y = inits3, times = times, func = f1, parms = param, method = "ode45")
  sol4 <- ode(y = inits4, times = times, func = f1, parms = param, method = "ode45")
  sol5 <- ode(y = inits5, times = times, func = f1, parms = param, method = "ode45")
  sol6 <- ode(y = inits6, times = times, func = f1, parms = param, method = "ode45")
  print(sol1)
  # highlightning modeled trajectories of selected entities
  matplot(sol1[,2], sol1[,3], type='l', col = 'blue', add=TRUE)
  #points(sol1[1,2], sol1[1,3], pch = 20, cex=1, col = 'blue')
  #lines(sol1[nrow(sol1),2], sol1[nrow(sol1),3], type='l', col = 'blue') 
  
  matplot(sol2[,2], sol2[,3], type='l', col = 'darkgreen', add=TRUE)
  #points(sol2[1,2], sol2[1,3], pch = 20, cex=1, col = 'darkgreen')
  #lines(sol2[nrow(sol2),2], sol2[nrow(sol2),3], type='l', col = 'darkgreen')
  
  matplot(sol3[,2], sol3[,3], type='l', col = 'red', add=TRUE)
  #points(sol3[1,2], sol3[1,3], pch = 20, cex=1, col = 'red')
  #lines(sol3[nrow(sol3),2], sol3[nrow(sol3),3], type='l', col = 'red')
  
  matplot(sol4[,2], sol4[,3], type='l', col = 'cyan', add=TRUE)
  #points(sol4[1,2], sol4[1,3], pch = 20, cex=1, col = 'cyan')
  #lines(sol4[nrow(sol4),2], sol4[nrow(sol4),3], type='l', col = 'cyan')
  
  matplot(sol5[,2], sol5[,3], type='l', col = 'magenta', add=TRUE)
  #points(sol5[1,2], sol5[1,3], pch = 20, cex=1, col = 'magenta')
  #lines(sol5[nrow(sol5),2], sol5[nrow(sol5),3], type='l', col = 'magenta')
  
  matplot(sol6[,2], sol6[,3], type='l', col = 'black', add=TRUE)
  #points(sol6[1,2], sol6[1,3], pch = 20, cex=1, col = 'black')
  #lines(sol6[nrow(sol6),2], sol6[nrow(sol6),3], type='l', col = 'black')
  
  dev.off
}