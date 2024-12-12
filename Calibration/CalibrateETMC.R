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


# Run simulation with # subinterval = 30, 50, 70, 100
# If still fail, report a highest cost
cost = function(a.input)
{
  if (any(a.input >= upper)) # Elasticitiy of transformation must less than 0
  {
    cost = Inf
  } else {
    df.input = data.frame(ETCP = df.ETCP$ETCP, 
                          ETMC.irri = a.input, 
                          ETMC.rain = a.input)
    
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
      
      df.sim$Qosd.irri =  GRIDDATA$qmcp[,1,"oilseed"]
      df.sim$Qosd.rain =  GRIDDATA$qmcp[,2,"oilseed"]
      
      df.sim$Qgro.irri =  GRIDDATA$qmcp[,1,"othergrains"]
      df.sim$Qgro.rain =  GRIDDATA$qmcp[,2,"othergrains"]
      
      df.sim$Qc_b.irri =  GRIDDATA$qmcp[,1,"sugarcrops"]
      df.sim$Qc_b.rain =  GRIDDATA$qmcp[,2,"sugarcrops"]
      
      
      df.sim$pQosd.irri = result$p_qmcropglc[,1,"oilseed",1]
      df.sim$pQosd.rain = result$p_qmcropglc[,2,"oilseed",1]
      
      df.sim$pQgro.irri = result$p_qmcropglc[,1,"othergrains",1]
      df.sim$pQgro.rain = result$p_qmcropglc[,2,"othergrains",1]
      
      df.sim$pQc_b.irri = result$p_qmcropglc[,1,"sugarcrops",1]
      df.sim$pQc_b.rain = result$p_qmcropglc[,2,"sugarcrops",1]
      
      df.sim$Qosd.sim = df.sim$Qosd.irri * (1+df.sim$pQosd.irri/100) + 
        df.sim$Qosd.rain * (1+df.sim$pQosd.rain/100)
      
      df.sim$Qgro.sim = df.sim$Qgro.irri * (1+df.sim$pQgro.irri/100) + 
        df.sim$Qgro.rain * (1+df.sim$pQgro.rain/100)
      
      df.sim$Qc_b.sim = df.sim$Qc_b.irri * (1+df.sim$pQc_b.irri/100) + 
        df.sim$Qc_b.rain * (1+df.sim$pQc_b.rain/100)
      
      # Calculate cost function, adjust by the total cropland and pasture area
      cost = (sum((df.target$Qosd.target - df.sim$Qosd.sim)^2) +
                sum((df.target$Qgro.target - df.sim$Qgro.sim)^2)*(price.corn / price.soy)^2 +
                sum((df.target$Qc_b.target - df.sim$Qc_b.sim)^2)*(price.cane / price.soy)^2) / num.gid.g
      
    } else {
      print("Model fails to run with maximum of subintervals")
      cost = Inf
    }
  }
  return(cost)
}

# Main script ----
# Remember to change the cost function with the weight of crop price square!

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
if(F)
{
  df.target$QLND.irri =  GRIDDATA$qlnd[,1]
  df.target$QLND.rain =  GRIDDATA$qlnd[,2]
  df.target$pQLND = df.LULC$cropland_shock
  df.target$QLND.target = (df.target$QLND.irri + df.target$QLND.rain)*(1+df.target$pQLND/100)
  
  df.target$QPLD =  GRIDDATA$qpld
  df.target$pQPLD = df.LULC$pasture_shock
  df.target$QPLD.target = (df.target$QPLD)*(1+df.target$pQPLD/100)
  
  total.cropland = sum(df.target$QLND.irri) + sum(df.target$QLND.rain) 
  total.pasture = sum(df.target$QPLD)
}

df.target$Qosd =  GRIDDATA$qmcp[,1,"oilseed"] + GRIDDATA$qmcp[,2,"oilseed"]
df.target$Qgro =  GRIDDATA$qmcp[,1,"othergrains"] + GRIDDATA$qmcp[,2,"othergrains"]
df.target$Qc_b =  GRIDDATA$qmcp[,1,"sugarcrops"] + GRIDDATA$qmcp[,2,"sugarcrops"]

df.target$pQosd = df.product$soy_shock
df.target$pQgro = df.product$corn_shock
df.target$pQc_b = df.product$cane_shock

df.target$Qosd.target =  df.target$Qosd * (1+df.target$pQosd/100)
df.target$Qgro.target =  df.target$Qgro * (1+df.target$pQgro/100)
df.target$Qc_b.target =  df.target$Qc_b * (1+df.target$pQc_b/100)

total.osd = sum(df.target$Qosd)
total.gro = sum(df.target$Qgro)
total.c_b = sum(df.target$Qc_b)
# Read in previous calibrated ETCP
df.ETCP = read.csv("calibrate/previous validation/ETCP.csv", header = T)
df.ETMC = read.csv("calibrate/previous validation/ETMC.csv", header = T)

df.price = read.csv("calibrate/in/BrazilCropPrice.csv", header = T)
price.soy = df.price[df.price$Item == "Soya beans",]$Price
price.corn = df.price[df.price$Item == "Maize (corn)",]$Price
price.cane = df.price[df.price$Item == "Sugar cane",]$Price
# Algorithm starts ----
num = num.gid.g
upper = 0
lower = -2

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

df.calib = data.frame(x = df.ETMC$ETMC)
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
ggsave(paste0("calibrate/out/", "ETMC_Error_min.png"))


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
ggsave(paste0("calibrate/out/", "ETMC_Error_max.png"))

write.csv(data.frame(GID = set.gid.g, ETMC = l[[1]]$x), "calibrate/out/ETMC.csv", row.names = F)

if(F)
{
  l.init = readRDS("calibrate/temp/NM_initial.RData")
  l.init = l.init[order(as.numeric(sapply(l.init, `[`, 2)))]
  
  l.init[[1]]
  l[[1]]
  
}


