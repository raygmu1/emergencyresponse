#Creating dataframe from output
library(jtools)
colnames(res_graph) <- c("focus_nodes_num", "focusresult", "baselinresult", "degree_cen", "eig_cen", "hub_cen", "authority_cen", "close_cen", "reach_cen", "between_cen", "coleman_s", "coleman_t", "eiindex", "skeptics_num", "trusting_num")

res_graph <- as.data.frame(res_graph)

#Running linear regression
focdata <- res_graph[,c(1,2,4:15)]
basedata <- res_graph[,c(3,4:15)]

#Scaled
sfocd <- as.data.frame(scale(focdata))
sbasd <- as.data.frame(scale(basedata))

#Removing any NaN columns
sfocd <- sfocd[ , colSums(is.na(sfocd)) == 0]
sbasd <- sbasd[ , colSums(is.na(sbasd)) == 0]

modfoc <- lm(focusresult ~ ., data = sfocd, na.action=na.omit)
modbas <- lm(baselinresult ~ ., data = sbasd, na.action=na.exclude)

summ(modfoc)

summ(modbas)

library(stargazer)

stargazer(modfoc, modbas, title="Results Barabasi Albert", align=TRUE)

library(edgynode)



export_summs(modfoc, modbas, scale = TRUE)