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
require(reshape2)

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

# Read in simulation results
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

# Get data for visualization

df.vis = subset(df.target, select = -c(QLND.irri,QLND.rain,pQLND, QPLD,pQPLD,
                                       Qosd, Qgro,Qc_b,pQosd,pQgro,pQc_b))
df.sim.refine = subset(df.sim, select = c(GID,QLND.sim,QPLD.sim,Qosd.sim,Qgro.sim,Qc_b.sim))

df.vis = join(df.vis, df.sim.refine)


visResult = function(df.vis.fig, var.name, var.unit)
{
  colnames(df.vis.fig) = c("ShortName","sim","target")
  df.vis.fig.melt = as.data.frame.array(melt(df.vis.fig, id = "ShortName"))
  df.vis.fig.melt$ShortName = factor(df.vis.fig.melt$ShortName, levels = df.vis.fig$ShortName)
  
  ggplot() + 
    geom_bar(data = df.vis.fig.melt, aes(x = ShortName, y = value, fill = variable), 
             position = "dodge", stat = "identity")+
    ylab(var.unit) +
    xlab("State") + 
    ggtitle(var.name)+
    scale_fill_manual(name = "Results", 
                      labels = c("Simulation", "Observation"), 
                      breaks = c("sim", "target"),
                      values=c("cornflowerblue","darkgoldenrod1"))+
    geom_hline(yintercept=0) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line.y = element_line(colour = "black"),
          text = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.title.y = element_text(size = 15))
  ggsave(paste0("visualization/",var.name,".png"),width=12, height=6)
}

df.vis.crop = subset(df.vis, select = c(ShortName, QLND.sim, QLND.target))
visResult(df.vis.crop, var.name = "Cropland area", var.unit = "1000 ha")

df.vis.pasture = subset(df.vis, select = c(ShortName, QPLD.sim, QPLD.target))
visResult(df.vis.pasture, var.name = "Pasture area", var.unit = "1000 ha")

df.vis.osd = subset(df.vis, select = c(ShortName, Qosd.sim, Qosd.target))
visResult(df.vis.osd, var.name = "Oilseeds production", var.unit = "1000 t")

df.vis.gro = subset(df.vis, select = c(ShortName, Qgro.sim, Qgro.target))
visResult(df.vis.gro, var.name = "Other grains production", var.unit = "1000 t")

df.vis.c_b = subset(df.vis, select = c(ShortName, Qc_b.sim, Qc_b.target))
visResult(df.vis.c_b, var.name = "Sugar crops production", var.unit = "1000 t")

write.csv(df.vis,"visualization/vis.csv", row.names = F)
