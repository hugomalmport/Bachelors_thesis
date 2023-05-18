
#Library
library(ggplot2)
library(ggpubr)
library(dplyr)
library(googlesheets4)
library(lme4)
library(lmerTest)
library(directlabels)


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

#####Do traits of fucus change on a depth gradient

#SA:P
m1 = lm(TDMC ~ depth_interpolated, data = data_fucus)
summary(m1)
ggscatter(x = "depth_interpolated", y = "TDMC",
          title = "",
          xlab = "Depth (cm)", ylab = "TDMC",
          color = "defining_morphology", palette = c("darkgreen", "orange", "lightgreen", "brown" ),
          data = data_fucus, add = "reg.line")
#TDMC
m2 = lm(SA_P ~ depth_interpolated, data = data_fucus)
summary(m2)
ggscatter(x = "depth_interpolated", y = "SA_P",
          title = "",
          xlab = "Depth (cm)", ylab = "SA_P",
          color = "defining_morphology", palette = c("darkgreen", "orange", "lightgreen", "brown" ),
          data = data_fucus, add = "reg.line")
#LW ratio
m3 = lm(LW_ratio ~ depth_interpolated, data = data_fucus)
summary(m3)
ggscatter(x = "depth_interpolated", y = "LW_ratio",
          title = "",
          xlab = "Depth (cm)", ylab = "Length width ratio",
          color = "defining_morphology", palette = c("darkgreen", "orange", "lightgreen", "brown" ),
          data = data_fucus, add = "reg.line")

#LP ratio
m4 = lm(LP_ratio ~ depth_interpolated, data = data_fucus)
summary(m4)
ggscatter(x = "depth_interpolated", y = "LP_ratio",
          title = "",
          xlab = "Depth (cm)", ylab = "Length perimeter ratio",
          color = "defining_morphology", palette = c("darkgreen", "orange", "lightgreen", "brown" ),
          data = data_fucus, add = "reg.line")

#STA(mm2/g-1)
m5 = lm(TDMC ~ depth_interpolated, data = data_fucus)
summary(m5)
ggscatter(x = "depth_interpolated", y = "STAmm2_g",
          title = "",
          xlab = "Depth (cm)", ylab = "STA(mm2/g)",
          color = "defining_morphology", palette = c("darkgreen", "orange", "lightgreen", "brown" ),
          data = data_fucus, add = "reg.line")

####Is there a difference between exposed and sheltered
m6 = lm(LW_ratio ~ Condition + length_mm, data = data_fucus)
summary(m6)

#LW_ratio
ggboxplot(x = "Condition", y = "LW_ratio", color = "Condition", 
          xlab = "Condition", ylab = "Length width ratio",
          palette = c("#FF6699", "#0099FF"),
          data = data_fucus)

#TDMC
m7 = lmer(TDMC ~ Condition + length_mm + (1 | transectlabel), data = data_fucus)
summary(m7)
ggscatter(x = "depth_interpolated", y = "TDMC", 
          color = "Condition", add = "reg.line",
          xlab = "Depth (cm)",
          data = data_fucus, 
          palette = c("#FF3333", "#66CC00"),
          ellipse = TRUE, mean.point = TRUE,
          star.plot = TRUE)

ggboxplot(x = "transectlabel", y = "TDMC", color = "Condition", data = data_fucus)

#Is it possible to identify germlings based on their traits
#Is serratus distinguishable by traits - most likely with area

ggboxplot(x = "Condition", y = "depth_interpolated", data = data_measurements)

#pca
library(car)
library(ggfortify)
library(MuMIn)

m8 = lmer(TDMC ~ scale(length_mm) + scale(depth_interpolated) + (1 | transectlabel), data = data_fucus)
summary(m8)
r.squaredGLMM(m1)
qqp(resid(m1))

data_pca = data_fucus %>% select(TDMC, SA_P, LW_ratio, LP_ratio, STAmm2_g, wetweight_g, depth_interpolated, Condition)
data_pca = na.omit(data_pca)
data_pca_result = prcomp(data_pca[c(-6, -7, -8)], scale = TRUE)

#PCA screeplot
eigenvalues <- data_pca_result$sdev^2
pcs <- 1:length(eigenvalues)
regression_line <- lm(eigenvalues ~ pcs)

screeplot(data_pca_result, main = "PCA", ylab = "Eigenvalues", xlab = "PCs")

variance_percentage <- eigenvalues / sum(eigenvalues) * 100
max_height <- max(data_pca_result$sdev)
label_position <- eigenvalues + 0.02 * max_height

labels <- c("38%", "29%", "17%", "12%", "4%")
label_pcs <- c(1, 2, 3, 4, 5)

text(label_pcs, label_position - 0.05 * max_height, labels, pos = 3, adj = c(0.5, 0.5))

abline(regression_line, col = "black")


data_pca_x = data.frame(data_pca_result$x)

ggscatter(x = "PC1", y = "PC2", color = "PC3", add = "reg.line", data = data_pca_x)

autoplot(data_pca_result, loadings = TRUE, loadings.label = TRUE, color = "Condition", size = "depth_interpolated", data = data_pca)


#Does the distribution of germlings look different from the adult distribution, filter for sheltered



