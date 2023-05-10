
library(googlesheets4)
#data=read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0")
#str(data)

#library(tidyverse)
#data_test = data %>% filter(transect_label=="A2")
#plot(data_test$Knot, data_test$`L/W ratio`)

#Library
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(readxl)
library(tidyverse)
library(ggfortify)
library(vegan)

#Read data
data_depth = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "depth")
data_measurements = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "measurements")

data_measurements$germling_label = as.character(data_measurements$germling_label)
data_measurements$germling_label[data_measurements$germling_label == "NULL"] = NA

# filter depth data and combine with other data
data_depth = data_depth %>% select(germling_label, depth_calculated)

data = left_join(data_measurements, data_depth, by = "germling_label")


#Plot
#ggplot(data, aes(x = TDMC, y = width_mm, color = presumed_species, size = depth(calculated)) +
 #  geom_point() +
  #theme_ipsum()


#get table with mean and sd per species
data_summary = data %>% group_by(presumed_species) %>% summarize(mean_depth = mean(depth_calculated, na.rm = TRUE),
                                                                 sd_depth = sd(depth_calculated, na.rm = TRUE))


ggplot() +
  geom_point(data = data_summary, 
             mapping = aes(x = presumed_species,
                           y = mean_depth, colour = presumed_species), 
             size = 2, shape = 19) +

geom_errorbar(data = data_summary,
                mapping = aes(x = presumed_species, colour = presumed_species, 
                              ymin = mean_depth - sd_depth,
                              ymax = mean_depth + sd_depth),
                width = 0.05)
 

 geom_quasirandom(data = all_depth,
                   mapping = aes(x = binomial_code, y = depth_correct, colour = binomial_code),
                   alpha = 0.3, shape = 1, width = 0.2)
  geom_point(data = all_depth_summary, 
             mapping = aes(x = binomial_code,
                           y = m_depth_correct, colour = binomial_code), 
             size = 2.75, shape = 18)


  
  
  ggplot(data = data, 
               mapping = aes(x = LW_ratio,
                             y = depth_calculated, color = presumed_species, size = length_mm), 
               size = 2, shape = 19) +
  geom_point() +
    xlim(c(0, 20))
  
  
plot(data %>% select(LW_ratio, depth_calculated, length_mm, `area(mm)`,`perimeter(mm)`,TDMC,LP_ratio,`SA:P`,`STA(mm(^2)g^-1))`))
round(cor(data %>% select(LW_ratio, depth_calculated, length_mm, `area(mm)`,`perimeter(mm)`,TDMC,LP_ratio,`SA:P`,`STA(mm(^2)g^-1))`), use = "pairwise.complete.obs"), 3)


data$presumed_species[is.na(data$presumed_species)] = "not yet defined"
data_pca = na.omit(data %>% select(TDMC, LP_ratio, `STA(mm(^2)g^-1))`, LW_ratio, presumed_species, length_mm, depth_calculated))



pca = prcomp(data_pca[c(-5, -6, -7)], center = TRUE, scale = TRUE)
screeplot(pca)
autoplot(pca, loadings = TRUE, loadings.label = TRUE, data = data_pca, color = "presumed_species", size = "depth_calculated")
adonis2(scale(data_pca[c(-5, -6, -7)])~depth_calculated , data = data_pca, permutations = 9999 , method = "euclidian")

        