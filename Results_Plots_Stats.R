##### Plots done with R Studio

# Librares
install.packages("devtools")
install.packages("HH")
library("devtools")
install_github("ndphillips/yarrr")
library("yarrr")
install.packages("survminer")
library("survminer")
library(ggplot2)
library(ggExtra)
library(dplyr)
library(magrittr)
library(data.table)


####################################################################
#
# Figure 3 (ABC): Pirate plots of heteroclonal xenograft composition
#
####################################################################

# Import data
Samples <-read_excel("Xenograft_composition.xlsx")

# Divide in quartiles by percentage of tumor
Samples$quartile <- smart_cut(Samples$Percentage, 4, "g", output = "numeric")

# Pirate Plots

#### Proliferation
plot3A<-pirateplot(formula = Alpha ~ quartile, 
                  data = Samples, 
                  xlab = "Quartile", 
                  ylab = "",
                  main = "Proliferation",
                  theme = 1,
                  # Choose your color palette, or give common color vector
                  pal = "southpark",
                  
                  # Set transparency of the elements:
                  avg.line.o = 0.9,
                  bar.f.o = .5,
                  bean.b.o = .9,
                  point.o = .5,
                  
                  # Shape of point
                  #point.pch = 2,
                  
                  #Background color
                  back.col = gray(.98),
                  quant = c(.25, .75), # 10th and 90th quantiles
                  quant.col = "black" # Black quantile lines
)
plot3A 


##### APNG

plot3B<-pirateplot(formula = APNG ~ quartile, 
                  data = Samples, 
                  xlab = "Quartile", 
                  main = "APNG", 
                  lab = "",
                  theme = 1,
                  # Choose your color palette, or give common color vector
                  pal = "southpark",
                  
                  # Set transparency of the elements:
                  avg.line.o = 0.9,
                  bar.f.o = .5,
                  bean.b.o = .9,
                  point.o = .5,
                  
                  # Shape of point
                  #point.pch = 2,
                  
                  #Background color
                  back.col = gray(.98),
                  quant = c(.25, .75), # 10th and 90th quantiles
                  quant.col = "black" # Black quantile lines
)
plot3B 


#### MGMT

plot3C<-pirateplot(formula = MGMT ~ quartile, 
                  data = Samples,
                  xlab = "Quartile", 
                  main = "MGMT", 
                  ylab = "",
                  theme = 1,
                  # Choose your color palette, or give common color vector
                  pal = "southpark",
                  
                  # Set transparency of the elements:
                  avg.line.o = 0.9,
                  bar.f.o = .5,
                  bean.b.o = .9,
                  point.o = .5,
                  
                  # Shape of point
                  #point.pch = 2,
                  
                  #Background color
                  back.col = gray(.98),
                  quant = c(.25, .75), # 10th and 90th quantiles
                  quant.col = "black" # Black quantile lines
)
plot3C




####################################################################
#
# Figure 3 (DE): Scatter plots of heteroclonal xenograft composition
#
####################################################################

# Import Data
# Data trimed to only contain most and least resistant
DataSet1 <- read_excel("Most_Least_Resistant.xlsx")

#Separate Data
AlphaW <- DataSet1$alphaWinner
APNGw <- DataSet1$APNGwinner
MGMTw <- DataSet1$MGMTwinner

AlphaL <- DataSet1$alphaLoser
APNGl <- DataSet1$APNGloser
MGMTl <- DataSet1$MGMTloser


# Plot with marginal histogram

plot3D<-ggplot(loserRepair, aes(x=APNGl, y=MGMTl, colour=AlphaL)) + 
  geom_point()+geom_rug(col="steelblue",alpha=0.1, size=1.5) + 
  geom_smooth(method=lm , color="gold", se=TRUE)+ 
  labs(title = "Most Sensitive",x = "APNG", y = "MGMT") + 
  theme(legend.position="none", plot.title = element_text(hjust=1))
ggMarginal(plot4b, type="histogram",col="navajowhite4",fill = "navajowhite")


plot3E<-ggplot(winnerRepair, aes(x=APNGw, y=MGMTw, colour=AlphaW)) + 
  geom_point()+geom_rug(col="steelblue",alpha=0.1, size=1.5) + 
  geom_smooth(method=lm , color="gold", se=TRUE) + 
  labs(title = "Most Resistant",x = "APNG", y = "MGMT") + 
  theme(legend.position="none", plot.title = element_text(hjust=1))
ggMarginal(plot1b, type="histogram",col="navajowhite4",fill = "navajowhite")




####################################################################
#
# Figure 4: Monoclonal time course
#
####################################################################

# Import Data

Null <- read_excel("Data_Null_Monoclonal.xlsx")
APNG <- read_excel("Data_APNG_Monoclonal.xlsx")
MGMT <- read_excel("Data_MGMT_Monoclonal.xlsx")
Both <- read_excel("Data_Both_Monoclonal.xlsx")



# For each treatment we overlap time courses of mean tumor sice and confidence intervals

MakePlot = function(x,Con,ConSD,TMZ,TMZSD,Opt,OptSD){
  df <- data.frame(x, Con,TMZ,Opt)
  
  eb1 <- aes(ymax = (Con + ConSD), ymin = pmax((Con - ConSD),0))
  eb2 <- aes(ymax = (TMZ + TMZSD), ymin = pmax((TMZ - TMZSD),0))
  eb3 <- aes(ymax = (Opt + OptSD), ymin = pmax((Opt - OptSD),0))
  
  plot1<-ggplot(df, aes(x)) +                    # basic graphical object
    geom_line(aes(y=Con), size = 1, colour="red4") + 
    geom_ribbon(eb1, alpha = 0.3, fill = "red")+
    geom_line(aes(y=TMZ), size = 1, colour="blue4",linetype = "dotted") + 
    geom_ribbon(eb2, alpha = 0.3, fill = "blue")+
    geom_line(aes(y=Opt), size = 1, colour="orange",linetype = "dashed") + 
    geom_ribbon(eb3, alpha = 0.3, fill = "gold") + labs(x = "Time", y = "Tumor Cells") + 
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
    theme(text = element_text(size=15)) 
  return(plot1)
  
}


# Plots

plot1a <- MakePlot(Null$Time, Null$Conventional_mean, Null$Conventional_SD ,Null$TMZ_mean, Null$TMZ_SD, Null$Optimal_mean, Null$Optimal_SD)
plot1b <- MakePlot(APNG$Time, APNG$Conventional_mean, APNG$Conventional_SD ,APNG$TMZ_mean, APNG$TMZ_SD, APNG$Optimal_mean, APNG$Optimal_SD)
plot1c <- MakePlot(MGMT$Time, MGMT$Conventional_mean, MGMT$Conventional_SD ,MGMT$TMZ_mean, MGMT$TMZ_SD, MGMT$Optimal_mean, MGMT$Optimal_SD)
plot1d <- MakePlot(Both$Time, Both$Conventional_mean, Both$Conventional_SD ,Both$TMZ_mean, Both$TMZ_SD, Both$Optimal_mean, Both$Optimal_SD)



# Inlets

plot2a <- plot1a + xlim(6, 7) +ylim(0,10000)
plot2b <- plot1b + xlim(6, 7) +ylim(0,20000)
plot2c <- plot1c + xlim(6, 7) +ylim(0,700000)
plot2d <- plot1d + xlim(6, 7) +ylim(0,1500000)





####################################################################
#
# Figure 5: Heteroclonal comparison of treatment strategies
#
####################################################################



# Import data

# For visualization this data set contains empty columns to add space between days
DataH <- read_excel("Heteroclonal_treatments.xlsx")
dodge <- position_dodge(width = 0.9)

# Set limits
limits <- aes(ymax = log10(DataH$Mean/10^5 + DataH$SD/10^5),
              ymin = log10(pmax(DataH$Mean/10^5 - DataH$SD/10^5,1e-5)))

# Create plot
p <- ggplot(DataH, aes(x=factor(DataH$Day), y=log10(DataH$Mean/10^5), fill=factor(DataH$Treatment)))

# Add bars, error bars, and labels
plot1 <-p + geom_bar(stat = "identity", width = 1, position = position_dodge(1)) +
            geom_errorbar(limits, position = position_dodge(1),width = 0.5,
            # Only color columns containing data
            colour=c("white","black","black","black","black","white","white","black","black","black","black","white","white","black","black","black","black","white","white","black","black","black","black","white")) +
            labs(x = "Day", y = "Fold Change") + 
            # Change y label
            scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6),labels=c("1e-6","1e-4", "1e-2", "1","1e2","1e4","1e6"))+ 
            theme(axis.text.x = element_text(color = "black", size = 14),
                  axis.text.y = element_text(color = "black", size = 14),
                  axis.text=element_text(size=14),
                  axis.title=element_text(size=14),
                  legend.text=element_text(size=14),
                  legend.title=element_text(size=14)) 

# Set colors manually
plot1 <- plot1+scale_fill_manual("Treatment",values = c("#A6CEE3","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99"), breaks = c("A", "B", "C","D"))

plot1
ggsave("Heteroclonal_Plot.png", plot = plot1)





