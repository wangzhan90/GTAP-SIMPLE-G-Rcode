# Calibrate key parameters with Nelder–Mead method
# Zhan Wang (zhanwang@purdue.edu)

# Setup ----
rm(list=ls())
#library('devtools')
#devtools::install_git('https://github.com/USDA-ERS/MTED-HARr.git')
#require(HARr)
library(HARr)
library(plyr)
library(ggplot2)

# Functions
# Add elements to list
addElement = function(list, variable, name, longname)
{
  attr(variable, 'description') = longname
  list = c(list, variable = list(variable))
  names(list)[length(list)] = name
  return(list)
}

# Run simulation
simulation = function(df.input, subinterval)
{
  print(paste0("Run simulation with subinterval = ",subinterval))
  # Simulate model
  # Use Shell to call batch file within R, which allows "wait"
  # use wait = T (default) to let R code wait until the batch file is finished
  # Use intern = T to capture console output as a R vector, and then turn it off with invisible()
  invisible(shell(paste0("do_",subinterval,".bat"),wait = T, intern = T))
  
  # Check success or not
  successFlag = F
  
  outputlog = "out/Calibration.log"
  if(file.exists(outputlog))
  {
    successNote = grep("The program has completed without error", readLines(outputlog), value = TRUE)
    if(length(successNote) > 0)
    {
      return(T)
    }
  }
  return(successFlag)
}


# Run simulation with # subinterval = 5, 10, 20, 30, 50, 70, 100
# If still fail, report a highest cost
cost = function(a.input)
{
  if (any(a.input >= upper)) # Elasticitiy of transformation must less than 0
  {
    cost = Inf
  } else {
    df.input = data.frame(ETCP = a.input, 
                          ETMC.irri = df.ETMC$ETMC, 
                          ETMC.rain = df.ETMC$ETMC)
    # To speed up, read the parameter from an independent file
    list.gridparm.cali = list()
    header.etcp = array(df.input$ETCP,
                        dim = c(length(set.gid.g)),
                        dimnames = list(GID = set.gid.g))
    list.gridparm.cali = addElement(list.gridparm.cali, header.etcp, "ETCP", "Elas. of trans. between cropland and pasture by grid")
    
    header.etmc = array(c(df.input$ETMC.irri, df.input$ETMC.rain),
                        dim = c(length(set.gid.g),2),
                        dimnames = list(GID = set.gid.g, LTYPE = set.ltyp))
    list.gridparm.cali = addElement(list.gridparm.cali, header.etmc , "ETMC", "Elas. of trans. between crop types by grid")
    write_har(list.gridparm.cali, 'in/CALIPARM.har', maxSize = 1e6)
    
    simFlag = simulation(df.input, 5)
    if(simFlag == F)
    {
      simFlag = simulation(df.input, 10)
      if(simFlag == F)
      {
        simFlag = simulation(df.input, 20)
        if(simFlag == F)
        {
          simFlag = simulation(df.input, 30)
          if(simFlag == F)
          {
            simFlag = simulation(df.input, 50)
            if(simFlag == F)
            {
              simFlag = simulation(df.input, 70)
              if(simFlag == F)
              {
                simFlag = simulation(df.input, 100)
              }
            }
          }
        }
      }
    }
    
    
    if(simFlag == T)
    {
      result = read_SL4("out/Calibration.sl4")
      
      df.sim = data.frame(GID = set.gid.g)
      df.sim$QLND.irri =  GRIDDATA$qlnd[,1]
      df.sim$QLND.rain =  GRIDDATA$qlnd[,2]
      
      df.sim$pQLND.irri = result$p_qlandgl[,1,1]
      df.sim$pQLND.rain = result$p_qlandgl[,2,1]
      
      df.sim$QLND.sim = df.sim$QLND.irri * (1+df.sim$pQLND.irri/100) + 
        df.sim$QLND.rain * (1+df.sim$pQLND.rain/100)
      
      df.sim$QPLD = GRIDDATA$qpld
      df.sim$pQPLD = result$p_qplandg[,1]
      df.sim$QPLD.sim = df.sim$QPLD*(1+df.sim$pQPLD/100)
      
      # Calculate cost function, adjust by the total cropland and pasture area
      cost = (sum((df.target$QLND.target - df.sim$QLND.sim)^2)*(total.pasture / total.cropland)^2 +
                sum((df.target$QPLD.target - df.sim$QPLD.sim)^2)) / num.gid.g
      
    } else {
      print("Model fails to run with maximum of subintervals")
      cost = Inf
    }
  }
  return(cost)
}

# Main script ----
GRIDSETS = read_har('in/GRIDSETS.har')
GRIDDATA = read_har('in/GRIDDATA.har')

set.gid.g = toupper(GRIDSETS$gid)
num.gid.g = length(set.gid.g)
set.ltyp = c("Irrigated","Rainfed")

# Read in shock from raw data
df.product = read.csv("calibrate/in/CropOutput.csv", header = T)
df.LULC = read.csv("calibrate/in/LULC.csv", header = T)

df.target = read.csv("calibrate/in/state.csv", header = T)
df.target$index = 1:nrow(df.target)
df.target$GID = paste0("G",sprintf("%06d", df.target$index))
# Recalculate value change from shock (in case the baseline data does not match)
df.target$QLND.irri =  GRIDDATA$qlnd[,1]
df.target$QLND.rain =  GRIDDATA$qlnd[,2]
df.target$pQLND = df.LULC$cropland_shock
df.target$QLND.target = (df.target$QLND.irri + df.target$QLND.rain)*(1+df.target$pQLND/100)

df.target$QPLD =  GRIDDATA$qpld
df.target$pQPLD = df.LULC$pasture_shock
df.target$QPLD.target = (df.target$QPLD)*(1+df.target$pQPLD/100)

total.cropland = sum(df.target$QLND.irri) + sum(df.target$QLND.rain) 
total.pasture = sum(df.target$QPLD)

# Read in previous calibrated values
df.ETCP = read.csv("calibrate/previous validation/ETCP.csv", header = T)
df.ETMC = read.csv("calibrate/previous validation/ETMC.csv", header = T)

# Algorithm starts ----
num = num.gid.g
upper = 0
lower = -1

alpha = 1
gamma = 2
rho = 1/2
sigma = 1/2

iter.max = 40
cost.tolerance = 1

# Step 1 initial parameters and their cost
l = list()

# For the first time: create n random point plus 1 calibrated point
for(i in 1:(num))
{
  print(paste0("Generate initial candidate: ",i))
  df = data.frame(x = runif(num, min = lower, max = upper))
  # The advantage of using a df instead of an array is that:
  # if I combine an array with a number in c(), they will be regarded as in the same array
  # using df helps me to separate the df for parameter and the cost function later
  l[[i]] = c(df, cost(df$x))
}

df.calib = data.frame(x = df.ETCP$ETCP)
l[[num+1]] = c(df.calib, cost(df.calib$x))

df.vis = data.frame(iter = 0, cost_min = 0, cost_max = 0)
# For saving purpose
saveRDS(l, file="calibrate/temp/NM.RData")
saveRDS(l, file="calibrate/temp/NM_initial.RData")
write.csv(df.vis,"calibrate/temp/df_vis.csv", row.names = F)
# For repeated case
if(F)
{
  l = readRDS(file="calibrate/temp/NM.RData")
  df.vis = read.csv("calibrate/temp/df_vis.csv", header = T)
}

iter = df.vis[nrow(df.vis),]$iter + 1
# Loop starts here
while (iter <= iter.max)
{
  print(paste0("Iteration: ",iter))
  
  # step 2: reorder parameters with the cost in ascending order
  l = l[order(as.numeric(sapply(l, `[`, 2)))]
  df.vis = rbind(df.vis, data.frame(iter = iter, cost_min = l[[1]][[2]], cost_max = l[[num+1]][[2]]))
  
  if(l[[1]][[2]] < cost.tolerance)
  {
    print("Optimum found")
    break
  }
  
  # step 3: calculate the centroid except n+1
  l.n = l[1:num]
  # a.o: the centroid
  a.o = array(rep(0, num))
  for (i in 1:length(l.n))
  {
    a.o = a.o + l.n[i][[1]]$x
  }
  a.o = a.o / num
  
  # step 4: calculate the reflection point
  a.r = a.o + alpha*(a.o - l[num+1][[1]]$x)
  cost.r = cost(a.r)
  
  # step 5: compare the cost of reflection point
  if(cost.r < l[[num]][[2]] & cost.r >= l[[1]][[2]]) # reflection
  {
    print("perform reflection")
    l[[num+1]] = c(data.frame(x = a.r), cost.r)
  } else if (cost.r < l[[1]][[2]]) { # expansion
    print("perform expansion")
    a.e = a.o + gamma*(a.r - a.o)
    cost.e = cost(a.e)
    if (cost.e < cost.r)
    {
      l[[num+1]] = c(data.frame(x = a.e), cost.e)
    } else {
      l[[num+1]] = c(data.frame(x = a.r), cost.r)
    }
  } else { # contraction
    if (cost.r < l[[num+1]][[2]]) # outside contraction or shrink
    {
      a.c = a.o + rho*(a.r - a.o)
      cost.c = cost(a.c)
      if (cost.c < cost.r)
      {
        print("perform outside contraction")
        l[[num+1]] = c(data.frame(x = a.c), cost.c)
      } else { # shrink
        print("perform shrink")
        for(i in 2:(num+1))
        {
          a.i = l[1][[1]]$x + sigma*(l[i][[1]]$x - l[1][[1]]$x)
          l[[i]] = c(data.frame(x = a.i), cost(a.i))
        }
      }
    } else { # inside contraction or shrink
      a.c = a.o + rho*(l[num+1][[1]]$x - a.o)
      cost.c = cost(a.c)
      if (cost.c < l[[num+1]][[2]])
      {
        print("perform inside contraction")
        l[[num+1]] = c(data.frame(x = a.c), cost.c)
      } else { # shrink
        print("perform shrink")
        for(i in 2:(num+1))
        {
          a.i = l[1][[1]]$x + sigma*(l[i][[1]]$x - l[1][[1]]$x)
          l[[i]] = c(data.frame(x = a.i), cost(a.i))
        }
      }
    }
  }
  
  # For saving purpose
  saveRDS(l, file="calibrate/temp/NM.RData")
  saveRDS(l, file=paste0("calibrate/temp/NM_", iter, ".RData"))
  write.csv(df.vis,"calibrate/temp/df_vis.csv", row.names = F)
  iter = iter + 1
  
  if (iter > iter.max)
  {
    print("Optimum does not found with given iterations")
  }
}

l = l[order(as.numeric(sapply(l, `[`, 2)))]

df.vis = read.csv("calibrate/temp/df_vis.csv", header = T)
ggplot()+            
  geom_point(data = df.vis[2:nrow(df.vis),], aes(iter, cost_min)) +                                      
  geom_line(data = df.vis[2:nrow(df.vis),], aes(iter, cost_min),
            linetype="dashed")+
  ylab("Error")+
  xlab("Iterations")+
  theme(
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black"),
    panel.background = element_blank(),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 15))
ggsave(paste0("calibrate/out/", "ETCP_Error_min.png"))


ggplot()+            
  geom_point(data = df.vis[2:nrow(df.vis),], aes(iter, cost_max)) +                                      
  geom_line(data = df.vis[2:nrow(df.vis),], aes(iter, cost_max),
            linetype="dashed")+
  ylab("Error")+
  xlab("Iterations")+
  theme(
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black"),
    panel.background = element_blank(),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 15))
ggsave(paste0("calibrate/out/", "ETCP_Error_max.png"))

write.csv(data.frame(GID = set.gid.g, ETCP = l[[1]]$x), "calibrate/out/ETCP.csv", row.names = F)








