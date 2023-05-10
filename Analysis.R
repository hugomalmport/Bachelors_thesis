

m1 = lm(LW_ratio ~ knot, data = data_measurements)
summary(m1)

plot(data_measurements$knot, data_measurements$LW_ratio)
library(ggpubr)
ggscatter(x = "knot", y = "LW_ratio", color = "presumed_species", data = data_measurements,
          add = "reg.line")
m1 = lm(LP_ratio ~ length_mm, data = data_measurements[data_measurements$presumed_species!="A. nodosum",])
summary(m1)


ggboxplot(x = "condition", y = "LP_ratio", data = data_measurements, color = "presumed_species")
