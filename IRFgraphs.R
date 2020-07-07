setwd('C:/Users/Nathan/Downloads/PerturbationMethods/Model2')
library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)
library(foreach)
library(iterators)

#I change the Ps to \U0001d4ab with MS Paint.

irfs = fread('IRFs/IRFstoGraph.csv')
irfs[,t := .I]
irfslong = melt(irfs,id.vars = 't',variable.name = 'irf',variable.factor = F)
irfslong[,shock := substr(irf,1,1)]
irfslong[,LambdaP := substr(irf,2,nchar(irf)-1)
       ][LambdaP == 0.8, policy := 'Active Policy'
       ][LambdaP == 0.95, policy := 'Inactive Policy']
irfslong[,var := substr(irf,nchar(irf),nchar(irf))]
graphorder = data.table(var = c('P','P','k','k','y','y','r','r'),
                        LambdaP = rep(c('0.8','0.95'), times = 4))
graphorder[,order := .I]
irfslong[graphorder, on = c('var','LambdaP'), order := i.order]
setkey(irfslong,order)
irfslong[,graphname := paste0(var,LambdaP)
       ][,graphname := factor(graphname,levels = unique(graphname))]

irfslong[shock == 'P', shockname := '1 standard deviation shock to \U0001d4ab'
       ][shock == 'Z', shockname := '1 standard deviation shock to \U03B6']

firstpart = function(x) if(substr(x,4,4)=='8'){
  paste0(substr(x,1,1),'     Active Policy: ')
} else {
  paste0(substr(x,1,1),'     Inactive Policy: ')
}
secondpart = function(x)paste0(' = ',substr(x,2,nchar(as.character(x))))

p_and_k = ggplot(irfslong[var %in% c('P','k') & t <= 60]) +
  geom_line(aes(x = t, y = value, linetype = shockname, colour = shockname),size = 0.8) +
  facet_wrap(~graphname,nrow = 2,scales='free',
    labeller = label_bquote(.(firstpart(graphname))*lambda[P]*.(secondpart(graphname)))) +
  # facet_wrap(~graphname,nrow = 2,scales='free', labeller = as_labeller(scriptP)) +
  labs(title = 'Figure 2.1a: Impulse Response Functions of \U0001d4ab and k\n')
  
p_and_k + theme_excel_new() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.2, color='grey' ),
        strip.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black'),
        legend.text = element_text(color = 'black',size = 11),
        panel.spacing = unit(2,'lines')) +
  scale_linetype_manual(values=c('solid', 'dashed')) +
  scale_color_manual(values=c('#777777', 'black')) +
  guides(linetype = guide_legend(nrow = 2, override.aes = list(size = 0.5)))
ggsave('IRFs/P_and_k_plot.png',dpi = 400)

r_and_y = ggplot(irfslong[var %in% c('r','y') & t <= 60]) +
  geom_line(aes(x = t, y = value, linetype = shockname, colour = shockname),size = 0.8) +
  facet_wrap(~graphname,nrow = 2,scales='free',
             labeller = label_bquote(.(firstpart(graphname))*lambda[P]*.(secondpart(graphname)))) +
  labs(title = 'Figure 2.1b: Impulse Response Functions of y and r\n')
r_and_y + theme_excel_new() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.2, color='grey' ),
        strip.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black'),
        legend.text = element_text(color = 'black'),
        panel.spacing = unit(2,'lines')) +
  scale_linetype_manual(values=c('solid', 'dashed')) +
  scale_color_manual(values=c('#777777', 'black')) +
  guides(linetype = guide_legend(nrow = 2, override.aes = list(size = 0.5)))
ggsave('IRFs/y_and_r_plot.png',dpi = 400)

# foreach(i = iter(unique(irfslong[,.(var,LambdaP)]),by='row'))%do%{
# ggplot(irfslong[var == i$var & LambdaP == i$LambdaP & t <= 60]) +
#   geom_line(aes(x = t, y = value, linetype = shockname, colour = shockname),size = 0.8) +
#   theme_excel_new() +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line( size=.3, color='#505050'),
#         strip.text = element_text(size = 12),
#         legend.position = 'none',
#         axis.text.x = element_text(size = 28),
#         axis.text.y = element_text(size = 28),
#         panel.background = element_rect(fill = 'white', colour = 'white')) +
#   scale_linetype_manual(values=c('solid', 'dashed')) +
#   scale_color_manual(values=c('#777777', 'black'))
# ggsave(paste0('IRFs/',i$var,as.character(100*as.numeric(i$LambdaP)),'.png'),dpi = 600)
# }