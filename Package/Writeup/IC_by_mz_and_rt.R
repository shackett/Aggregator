
library(dplyr)

MZ_range <- c(150, 1200)
MZ_cv <- 0.002
MZ_sd <- 1.5
RT_range <- c(10, 60)
RT_cv <- 0.01
RT_sd <- 0.4
ICmean <- 1e6

nions <- 500

sample_ions <- data.frame(
  MZ = runif(nions, MZ_range[1], MZ_range[2]),
  RT = runif(nions, RT_range[1], RT_range[2]),
  IC = rgamma(nions, shape = 0.2, scale = ICmean)
) %>% tbl_df()

sample_ions <- sample_ions %>% mutate(MZsd = MZ_sd + MZ*MZ_cv,
                                      RTsd = RT_sd + RT*RT_cv)

#sample_ions <- sample_ions %>% mutate(MZsd = MZ*MZ_cv,
#                                      RTsd = RT*RT_cv)

sample_ions <- sample_ions %>% filter(IC > 1000) %>%
  mutate(ionID = 1:n())

# Place m/s dist on a grid
MZ_plot <- seq(MZ_range[1], MZ_range[2], by = 0.3)
RT_plot <- seq(RT_range[1], RT_range[2], by = 0.3)

MZ_RT_grid <- expand.grid(testMZ = MZ_plot, testRT = RT_plot, ionID = sample_ions$ionID) %>% tbl_df()

MZ_RT_grid <- MZ_RT_grid %>% left_join(sample_ions, by = "ionID")

MZ_RT_grid <- MZ_RT_grid %>% filter(abs(MZ - testMZ)/MZsd < 5,
                                    abs(RT - testRT)/RTsd < 5)

MZ_RT_grid <- MZ_RT_grid %>% mutate(ionDensity = dnorm(testMZ, mean = MZ, sd = MZsd) *
                        dnorm(testRT, mean = RT, sd = RTsd) * IC)

MZreduced <- MZ_RT_grid %>% dplyr::select(testMZ, testRT, ionDensity) %>%
  group_by(testMZ, testRT) %>% dplyr::summarize(IC = sum(ionDensity))

# fill out grid
possible_grid_vals <- expand.grid(testMZ = MZ_plot, testRT = RT_plot) %>% tbl_df()

MZreduced <- rbind(MZreduced,
                   possible_grid_vals %>% anti_join(MZreduced, by = c("testMZ", "testRT")) %>% mutate(IC = 0))

library(tidyr)
library(rgl)

MZRTgrid <- MZreduced %>% spread(testMZ, IC)
RTvals <- MZRTgrid$testRT
MZRTgrid <- MZRTgrid %>% dplyr::select(-testRT) %>% as.data.frame() %>% as.matrix()
MZvals <- as.numeric(colnames(MZRTgrid))

nbcol = 100
color = colorRampPalette(c("blue", "red"))(nbcol)
color = colorRampPalette(c("gray50", "dodgerblue4", "deepskyblue"))(nbcol)
color = colorRampPalette(c("gray50", "dodgerblue4", "deepskyblue", "firebrick2"))(nbcol)
zcol  = cut(MZRTgrid, nbcol)

persp3d(x = RTvals, y = MZvals, z = MZRTgrid, col= color[zcol])

# save
rgl.postscript(filename = "~/Desktop/Rabinowitz/Aggregator/Package/Writeup/Figures/wireframe.pdf", fmt = "pdf")
