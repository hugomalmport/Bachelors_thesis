
#Library
library(ggplot2)
library(ggpubr)
library(dplyr)
library(googlesheets4)
library(lme4)
library(lmerTest)


data_measurements = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "measurements")
data_measurements$transect_pos_cm = data_measurements$knot*10-10


#Add missing depths
data_depth = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "depth")


# here i recalculate the depth 
data_depth$depth_corrected = data_depth$`viva(cm)` + data_depth$`depth_measured (cm)`

# here i remove the depth corrected when there is no measurement
data_depth$depth_corrected[is.na(data_depth$`depth_measured (cm)`)]=NA
data_depth$transect_pos_cm = data_depth$knot*10-10


library(tidyr)
library(dplyr)


# Step 1: Interpolation within known points
data_depth <- data_depth %>%
  group_by(transectlabel) %>%
  arrange(transect_pos_cm) %>%
  mutate(depth_interpolated = approx(transect_pos_cm, depth_corrected, transect_pos_cm, rule = 1)$y)

# Step 1.5: Calculate slope between each pair of known points
data_depth <- data_depth %>%
  mutate(slope_interpolated = (depth_interpolated - lag(depth_interpolated)) / (transect_pos_cm - lag(transect_pos_cm)))

# Step 2: Extrapolation for the edges
data_depth <- data_depth %>%
  ungroup() %>%
  arrange(transectlabel, transect_pos_cm)

# For-loop for start points
for(i in unique(data_depth$transectlabel)) {
  # Calculate slope of the first non-NA pair
  slope_start <- data_depth$slope_interpolated[data_depth$transectlabel == i & !is.na(data_depth$slope_interpolated)][1]
  
  # Find the indices to be replaced
  indices <- which(data_depth$transectlabel == i & is.na(data_depth$depth_interpolated) & data_depth$transect_pos_cm < min(data_depth$transect_pos_cm[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)]))
  
  # Replace the values
  data_depth$depth_interpolated[indices] <- data_depth$depth_interpolated[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)][1] - slope_start * (data_depth$transect_pos_cm[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)][1] - data_depth$transect_pos_cm[indices])
}

# For-loop for end points
for(i in unique(data_depth$transectlabel)) {
  # Calculate slope of the last non-NA pair
  slope_end <- data_depth$slope_interpolated[data_depth$transectlabel == i & !is.na(data_depth$slope_interpolated)][length(data_depth$slope_interpolated[data_depth$transectlabel == i & !is.na(data_depth$slope_interpolated)])]
  
  # Find the indices to be replaced
  indices <- which(data_depth$transectlabel == i & is.na(data_depth$depth_interpolated) & data_depth$transect_pos_cm > max(data_depth$transect_pos_cm[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)]))
  
  # Replace the values
  data_depth$depth_interpolated[indices] <- data_depth$depth_interpolated[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)][length(data_depth$depth_interpolated[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)])] + slope_end * (data_depth$transect_pos_cm[indices] - data_depth$transect_pos_cm[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)][length(data_depth$depth_interpolated[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)])])
}

# In these for-loops, we iterate over each unique transect. For each transect, we first compute the slope at the start (or end) using the first (or last) non-NA interpolated depth value. We then identify the indices of the depth_interpolated values to be replaced, and replace them with the extrapol

data_measurements = left_join(data_measurements,data_depth, by =c("transectlabel","transect_pos_cm"))
data_measurements$depth_interpolated = round(data_measurements$depth_interpolated, 2)


######


data_measurements$depth_interpolated = data_measurements$depth_interpolated*-1



data_fucus = data_measurements %>% filter(defining_morphology != "A. nodosum")

#Do traits of fucus change on a depth gradient
m1 = lm(TDMC ~ depth_interpolated, data = data_fucus)
summary(m1)
ggscatter(x = "depth_interpolated", y = "TDMC",
          title = "",
          xlab = "Depth (cm)", ylab = "LW_ratio",
          color = "defining_morphology", palette = c("darkgreen", "orange", "lightgreen", "brown" ),
          data = data_fucus, add = "reg.line")


#Is there a difference between exposed and sheltered
m1 = lm(LW_ratio ~ Condition   * length_mm, data = data_fucus)
summary(m1)


ggboxplot(x = "Condition  ", y = "LW_ratio", data = data_fucus)
ggscatter(x = "length_mm", y = "LW_ratio", color = "Condition", data = data_fucus)


m1 = lmer(TDMC ~ Condition + length_mm + (1 | transectlabel), data = data_fucus)
summary(m1)
ggscatter(x = "depth_interpolated", y = "TDMC", 
          color = "Condition", add = "reg.line", 
          data = data_fucus, 
          palette = c("red", "darkgreen", "grey", "blue"),
          ellipse = TRUE, mean.point = TRUE,
          star.plot = TRUE)

ggboxplot(x = "transectlabel", y = "TDMC", color = "Condition ", data = data_fucus)

#Is it possible to identify germlings based on their traits
#Is serratus distinguishable by traits - most likely with area

ggboxplot(x = "Condition", y = "depth_interpolated", data = data_measurements)

#pca
library(car)
library(ggfortify)
library(MuMIn)

m1 = lmer(TDMC ~ scale(length_mm) + scale(depth_interpolated) + (1 | transectlabel), data = data_fucus)
summary(m1)
r.squaredGLMM(m1)
qqp(resid(m1))


data_pca = data_fucus %>% select(TDMC, SA_P, LW_ratio, LP_ratio, STAmm2_g, wetweight_g, depth_interpolated, Condition)
data_pca = na.omit(data_pca)
data_pca_result = prcomp(data_pca[c(-6, -7, -8)], scale = TRUE)
screeplot(data_pca_result)
data_pca_x = data.frame(data_pca_result$x)

ggscatter(x = "PC1", y = "PC2", color = "PC3", add = "reg.line", data = data_pca_x)

autoplot(data_pca_result, loadings = TRUE, loadings.label = TRUE, color = "Condition", size = "depth_interpolated", data = data_pca)


#Does the distribution of germlings look different from the adult distribution, filter for sheltered



######


m1 = lm(LW_ratio ~ knot, data = data_measurements)
summary(m1)

plot(data_measurements$knot, data_measurements$LW_ratio)
library(ggpubr)
ggscatter(x = "knot", y = "LW_ratio", color = "defining_morphology", data = data_measurements,
          add = "reg.line")
m1 = lm(TDMC ~ length_mm, data = data_measurements[data_measurements$defining_morphology!="A. nodosum",])
summary(m1)


ggboxplot(x = "Condition ", y = "STAmm2_g", data = data_measurements, color = "defining_morphology")
