library(ggplot2)
library(cowplot)
library(openxlsx)

setwd('C:/Users/Nathan/Downloads/PerturbationMethods/Model2')

pdata = t(fread('HistoricalEndogenousP.csv'))
historical_p = fread('Historical_True_P_Path.csv')
quants <- c(0.005,0.025,0.05,0.1,0.25,0.50,0.75,0.90,0.95,0.975,0.995)
widetable_with_pctiles = data.table(year=1985:2017,trueP = historical_p$V1,t(apply( pdata , 2 , quantile , probs = quants)))
widetable_with_pctiles[,P_label:='Historical P'
                     ][,P_label2:='Model Median']

ggplot(widetable_with_pctiles) +
  geom_ribbon(aes(x = year,ymin = `25%`, ymax = `75%`,colour = '50%'), alpha = 0.36,fill='blue', linetype = 0) +
  geom_ribbon(aes(x = year,ymin = `10%`, ymax = `90%`,colour = '80%'), alpha = 0.25,fill='blue', linetype = 0) +
  geom_ribbon(aes(x = year,ymin = `5%`, ymax = `95%`,colour = '90%'), alpha = 0.15,fill='blue', linetype = 0) +
  geom_ribbon(aes(x = year,ymin = `2.5%`, ymax = `97.5%`,colour = '95%'), alpha = 0.08,fill='blue', linetype = 0) +
  geom_ribbon(aes(x = year,ymin = `0.5%`, ymax = `99.5%`,colour = '99%'), alpha = 0.04,fill='blue', linetype = 0) +
  geom_line(aes(x = year,y = trueP,linetype = P_label),size = 1.2,colour='#c05020') +
  geom_line(aes(x = year,y = `50%`,linetype = P_label2),size = 1.2,colour = 'black')+
  # scale_fill_manual("",values="black") +
  scale_x_continuous(limits = c(1985,2017), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_cowplot() +
  labs(y= "P", x = "Year") +
  theme(axis.title.y = element_text(angle=0, vjust = 0.5),
        legend.title = element_blank(),
        legend.justification=c(-0.2,1),legend.position=c(0,1)) +
  guides(colour = guide_legend(override.aes = list(alpha = c(0.18,0.12,0.07,0.04,0.02))),
         linetype = guide_legend(override.aes = list(colour = c('#c05020','black'))))

ggsave('Historical_P_Plot.png')

wb=createWorkbook()
addWorksheet(wb, 'Percentiles')
addWorksheet(wb, 'Raw Data')
writeData(wb, 'Percentiles',widetable_with_pctiles[,`:=`(P_label=NULL,P_label2=NULL)])
writeData(wb, 'Raw Data',pdata)
saveWorkbook(wb, 'HistoricalP.xlsx',overwrite = T)