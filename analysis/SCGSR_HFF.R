## ----setup, include=FALSE------------------------------------------------------------------
library(formatR)
knitr::opts_chunk$set(echo = TRUE, fig.align = "center", fig.show = "hold", time_it = TRUE, dpi = 100)
knitr::opts_chunk$set(tidy.opts = list(width.cutoff = 60), tidy = T, collapse = TRUE, echo = FALSE, message = FALSE, warning = FALSE)
knitr::opts_knit$set(root.dir = '/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/')
library(terra)
library(lme4)
library(MASS)
library(lmerTest)
library(tidyverse)
library(sf)
library(ggplot2)
library(stats)
library(ggcorrplot)
library(RColorBrewer)
library(webshot)
library(kableExtra)
library(formatR)
library(stringr)
library(data.table)
library(openxlsx)


#knitr::knit_hooks$set(webgl = hook_webgl)
rgl::setupKnitr(autoprint = TRUE)

terraOptions(
    memfrac = 0.1
)

setGDALconfig("GDAL_PAM_ENABLED", "FALSE")


## ------------------------------------------------------------------------------------------
HLEF_Full <- openxlsx::read.xlsx("/Users/Anthony/OneDrive - UW/University of Washington/Lab Work/HLEF2022/HLEF2022LAB.xlsx", sheet = "HLEF_FullData")

HLEF_HFF <- HLEF_Full |> filter(!is.na(HFF_percC)) |> 
    mutate(HFF_mgC_gSoil = (HFF_percC/100)*(HFF/df_start_wt)*1000,
           HFF_percC_of_soilC = ((HFF_percC/100)*HFF)/((carbon_perc/100)*df_start_wt),
           siltclay = silt+clay) 

readr::write_csv(HLEF_HFF, "SCGSR/data/dataframes/HlEF_HFF.csv")


## ------------------------------------------------------------------------------------------
HLEF_WIP <- 1-rast("HLEF/SEAK_WIP_9-10-2021_output.tif")
plot(HLEF_WIP)
HLEF_geo <- vect("HLEF/Geology/Southeast_Alaska_Geology_(TNC)_/Southeast_Alaska_Geology_(TNC)_.shp") |> terra::crop(ext(HLEF_WIP)) |> 
    tidyterra::mutate(label = as.factor(label), 
                      LITH = case_when(str_detect(label, "KJs") ~ "Slate",
                                       str_detect(label, "KP") ~ "Phyllite",
                                       str_detect(label, "TK") ~ "Tonalite",
                                       str_detect(label, "Qs") ~ "Quaternary",
                                       str_detect(label, "KJv") ~ "Metavolcanic",
                                       str_detect(label, "g") ~ "Glacier",
                                       .default = "Quaternary"))

HLEF_HFF_pts <- HLEF_HFF |> vect(geom = c("lon", "lat")) |> 
    terra::extract(x = HLEF_WIP, bind = TRUE) |> 
    tidyterra::rename(WIP = SEAK_WIP_9.10.2021_output) 

HLEF_geo_ext <- terra::extract(HLEF_HFF_pts, x = HLEF_geo)

HLEF_HFF_pts_df <- HLEF_HFF_pts |> data.frame() |> 
    cbind(HLEF_geo_ext$LITH)  |>
    dplyr::rename("LITH" = "HLEF_geo_ext$LITH") |> 
    mutate(LITH = case_when(is.na(LITH) ~ "Quaternary",
                           .default = LITH)) |> 
    filter(HFF_percC_of_soilC < 1) 

readr::write_csv(HLEF_HFF_pts_df, file = "HLEF/Dataframes/HLEF_HFF_pts_df.csv")

HLEF_HFF_pts_df$LITH <- factor(HLEF_HFF_pts_df$LITH)

ggplot(HLEF_HFF_pts_df, aes(x = LITH, y = HFF_mgC_gSoil)) + 
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
    geom_jitter(width = 0.1) +
    theme(text = element_text(size = 14),
        panel.background = element_blank(),
            panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank())


## ------------------------------------------------------------------------------------------
HLEF_HFF_pts_sum <- HLEF_HFF_pts |> data.frame() |>  dplyr::group_by(sample_ID) |> dplyr::summarise(sum = sum(HFF_mgC_gSoil),
                                                                                wip = mean(WIP)) 

hist(HLEF_HFF_pts_sum$wip, breaks = 20)




## ------------------------------------------------------------------------------------------

hff_lm <- lm(HFF_mgC_gSoil ~ siltclay + ph + cn +WIP, 
             data = HLEF_HFF_pts |> 
                 filter(HFF_percC_of_soilC < 1))
summary(hff_lm)

hff_lmer <- lmer((HFF_mgC_gSoil) ~ siltclay + cn + WIP +(1|sample_ID), 
             data = HLEF_HFF_pts_df, REML = F)
summary(hff_lmer)

HLEF_HFF |> 
    filter(HFF_percC_of_soilC < 1) |> 
ggplot(aes(x = siltclay, y = HFF_mgC_gSoil, colour = ph)) +
    geom_point( size = 3.5) + 
    geom_smooth( method = "lm", se = F, fullrange = TRUE) +
    geom_abline(slope = 0.86, intercept = 0, colour = "darkgreen", lty = 2, size = 1) +
    geom_abline(slope = 0.48, intercept = 0, colour = "black", lty = 2, size = 1) +
    scale_x_continuous(limits = c(0, 101), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 101), expand = c(0,0), name = "Observed HFF SOC mgC/gSoil") +
    ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "adj.R2", "P")), size = 4) +
    theme(text = element_text(size = 14),
          panel.background = element_blank(),
            panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank())

HLEF_HFF |> 
    filter(HFF_percC_of_soilC < 1) |> 
ggplot(aes(x = siltclay, y = HFF_mgC_gSoil, colour = cn)) +
    geom_point( size = 3.5) + 
    geom_smooth( method = "lm", se = F, fullrange = TRUE) +
    geom_abline(slope = 0.86, intercept = 0, colour = "darkgreen", lty = 2, size = 1) +
    geom_abline(slope = 0.48, intercept = 0, colour = "black", lty = 2, size = 1) +
    scale_x_continuous(limits = c(0, 101), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 101), expand = c(0,0), name = "Observed HFF SOC mgC/gSoil") +
    ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "adj.R2", "P")), size = 4) +
    theme(text = element_text(size = 14),
          panel.background = element_blank(),
            panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank())

HLEF_HFF_pts_df |> data.frame() |> 
    filter(HFF_percC_of_soilC < 1) |> 
    ggplot(aes(x = siltclay, y = HFF_mgC_gSoil, colour = LITH)) +
    geom_point( size = 3.5) + 
    geom_smooth( method = "lm", se = F, fullrange = FALSE) +
    geom_abline(aes( slope = 0.86, intercept = 0, colour = "darkgreen"), lty = 2, size = 1) +
    geom_abline(aes(slope = 0.48, intercept = 0, colour = "darkblue"), lty = 2, size = 1) +
    scale_colour_manual(name = "Lithology", 
                          values = c("steelblue", "brown", "gray25", "darkgreen", "darkblue")) +
    scale_x_continuous(limits = c(0, 101), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 101), expand = c(0,0)) +
    #ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "adj.R2", "P")), size = 4) +
    
    theme(text = element_text(size = 14),
          panel.background = element_blank(),
            panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank())

HLEF_HFF_pts_df |> data.frame() |> 
    filter(HFF_percC_of_soilC < 1) |> 
    ggplot(aes(x = WIP, y = HFF_mgC_gSoil, colour = siltclay)) +
    geom_point( size = 3.5) + 
    geom_smooth( method = "lm", se = F, fullrange = FALSE) +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 45),expand = c(0,0)) +
    #ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "adj.R2", "P")), size = 4) +
    
    theme(text = element_text(size = 14),
          panel.background = element_blank(),
            panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank())


ggplot(aes(x = HFF_mgC_gSoil, y = predict(hff_lmer, newdata = HLEF_HFF_pts_df, allow.new.levels = T), colour = siltclay), data = HLEF_HFF_pts_df) +
    geom_point( size = 3.5) + 
    geom_smooth( method = "lm", se = F, fullrange = TRUE, color = "red") +
    geom_abline(slope = 1, intercept = 0, colour = "black", lty = 2, size = 1) +
    #geom_abline(slope = 0.48, intercept = 0, colour = "black", lty = 2, size = 1) +
    scale_x_continuous(limits = c(0, 50), expand = c(0,0), name = "Predicted HFF SOC mgC/gSoil") +
    scale_y_continuous(limits = c(0, 50), expand = c(0,0), name = "Observed HFF SOC mgC/gSoil") +
    scale_colour_viridis_c(option = "mako", name = "Silt + Clay %") +
    ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "adj.R2", "P")), size = 4) +
    theme(text = element_text(size = 14),
          panel.background = element_blank(),
            panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank())


## ------------------------------------------------------------------------------------------
ICP_DC <- openxlsx::read.xlsx("/Users/Anthony/OneDrive - UW/University of Washington/Lab Work/SCGSR/data/ICP_MS/20241107 Organic Samples for Keith.xlsx", sheet = 3) |> 
    dplyr::filter(HLEF_ID != "Blank")

CAMS_14C <- openxlsx::read.xlsx("/Users/Anthony/OneDrive - UW/University of Washington/Lab Work/SCGSR/data/14C/14C_Dataframe.xlsx", sheet = 1) |> 
    dplyr::rename("delta_14C" = "∆14C")

HLEF_HFF_pts_df |> left_join(y = CAMS_14C, by = join_by("HLEF_ID")) |> 
    dplyr::left_join(ICP_DC, by = join_by("HLEF_ID")) |> 
    filter(!is.na(delta_14C)) |> 
    mutate(wetupl = case_when(WIP >= 0.50 ~ "wet",
                              .default = "upl")) |> 
    ggplot(aes(x = (Si), y = delta_14C, colour = ph)) +
    geom_point() + 
    #geom_violin() #+
    geom_smooth(method = "lm")





## ------------------------------------------------------------------------------------------
cor_dat <- HLEF_HFF_pts_df |> left_join(y = CAMS_14C, by = join_by("HLEF_ID")) |> 
    #dplyr::left_join(ICP_DC, by = join_by("HLEF_ID")) |> 
    filter(!is.na(delta_14C)) |> 
    dplyr::select(where(is.numeric)) |> as.matrix()

corplot <- ggcorrplot(cor(cor_dat), method = "square", type = "full", lab = T, lab_size = 1.5)

corplot

