# ---- Notes & Reference Calculations ----
#demo groups income, % black, % hispanic, % white
#weight by county pop and number of observations - alpha and/or size OR make it continuous variable
#cutoffs: all march/apr 2019 vs 2020 and by demo in each group
#half of each demo below 50% in march/apr 2019 and half of each demo above 50% in march/apr 2020 - RD but artif. restr. data
#variance among all samples of each demo group between march/apr 2019 and 2020 - good if relevant
#number of unique sites per county / county pop (sum unique block group pops)

#reference for observations by data set and geog delineation
#list(c(length(unique(PM$site_ID)), length(unique(CALIFORNIA_PM.daily.EJ$site_ID)), 
#       length(unique(PM$county)), length(unique(CALIFORNIA_PM.daily.EJ$county)), 
#       length(unique(PM$block_group)), length(unique(CALIFORNIA_PM.daily.EJ$block_group)),
#       length(unique(PM$tract)), length(unique(CALIFORNIA_PM.daily.EJ$tract))))

#date range: 20 march - 30 april

#race composition indicators:
#diverse pop (none exceed 50%): 0
#vast majority white >80%: 1
#majority white >50: 2
#vast majority black >80%: 3
#majority black >50%: 4
#vast majority hispanic (any+) >80%: 5
#majority hispanic (any) >50%: 6
#vast majority asian >80%: 7
#majority asian >50% 8

#DP colors
#c(0 = "gray25", 1 = "darkorange3", 2 = "darkorange2", 3 = "brown4", 4 = "brown3", 
#5 = "deepskyblue3", 6 = "deepskyblue", 7 = "springgreen4", 8 = "springgreen2")

#summarize weights
#summary(PM_CTWt$CTWt7Raw)
#summary(PM_CoWt$CoWt7Raw)
#summary(PM_CTWt$CTWt7_01) 
#summary(PM_CoWt$CoWt7_01) 
#summary(PM_CTWt$CTWt7_1c)
#summary(PM_CoWt$CoWt7_1c)

#Estimating Equation variables:
# E = Exposure
# P = PM2.5
# R = race group
# sub i = location: sensor-level (daily PM2.5 sub is) and block-level (race group percent sub ib)
# sub t = time: year (same time period)
# sub b = block level
# sub s = sensor level
# sub d = daily level
# $$E_{it} = P_{is,td} + R_{ib,t}$$


# ---- Prep + Reference ----

setwd("~/Downloads/AQRD/aqrd-final")

PM <- read.csv("data/PM.csv")
PMlong_bins_all <- read.csv("data/PMlong_bins_all.csv")

library(tidyverse)
library(lubridate)
library(ggplot2)
library(vroom)
library(stargazer)
library(webshot2)
library(magick)
library(hrbrthemes)
library(viridis)
library(forcats)
library(magrittr)
library(gt)
library(imagefx)

detach("package:stargazer",unload=T)
# Delete it
remove.packages("stargazer")
# Download the source
download.file("https://cran.r-project.org/src/contrib/stargazer_5.2.3.tar.gz", destfile = "stargazer_5.2.3.tar.gz")
# Unpack
untar("stargazer_5.2.3.tar.gz")
# Read the sourcefile with .inside.bracket fun
stargazer_src <- readLines("stargazer/R/stargazer-internal.R")
# Move the length check 5 lines up so it precedes is.na(.)
stargazer_src[1990] <- stargazer_src[1995]
stargazer_src[1995] <- ""
# Save back
writeLines(stargazer_src, con="stargazer/R/stargazer-internal.R")
# Compile and install the patched package
install.packages("stargazer", repos = NULL, type="source")

library(stargazer)



PM_Obs <- PM |> #reference for number of observations per county
  count(county)
PM_ObsPs <- PM |> #reference for number of observations per sensor
  count(site_ID)
AvObs <- mean(PM_ObsPs$n) 


# ---- County-level PM Recode + Weighting ----

PM_CoWt <- PM |> #create sensor-level weights for data availability
  group_by(site_ID) |>
  mutate(SeWt = length(site_ID) / AvObs) |>
  ungroup()
PM_CoWt <- PM_CoWt |> #average sensor-level weights by county, total county population served, and create county-level weights
  group_by(county) |>
  mutate(CoSeWt = mean(SeWt),
         CoPopSe = sum(unique(race_total_block_group)),
         CoWt7Raw = length(county) / CoPopSe^.7 * CoSeWt) |>
  ungroup() |>
  mutate(CoWt7_01 = CoWt7Raw / max(CoWt7Raw), #standardize weight between >0 and ≤1
         CoWt7_1c = (CoWt7Raw / mean(CoWt7Raw)) - 1, #standardize around 1 at avg then temporarily to 0
         CoWt7_1c = (CoWt7_1c / 2) + 1, #reduce range of weighting then center back around 1
         is20 = case_when(year == 2020 ~ 1, TRUE ~ 0)) #dummy variable, 1 for 2020, 0 for 2019
as.factor(PM_CoWt$is20)

#filter for Mar 20 to April 30, 2019 and 2020 in new DF, recode 2020 weights to be above 1, filter for only PM2.5 parameter
PM_CoWtSel <- PM_CoWt |>
  filter(month == 3 & day > 19 & parameter == "PM2.5" | month == 4 & parameter == "PM2.5") |>
  mutate(CoWt7_01Spl = case_when(year == 2020 ~ 2-CoWt7_01, TRUE ~ CoWt7_01)) |>
  mutate(pctWhite = race_white_block_group / race_total_block_group) 
PM_CoWtSel <- PM_CoWtSel |> #generate county-level yearly averages
  group_by(county, year) |>
  mutate(CoAvYr = mean(sample_mean)) |>
  ungroup()

PM_CWS_WM <- PM_CoWtSel |> #Calculate weighted mean by year
  group_by(year) |>
  summarize(WM = weighted.mean(CoAvYr, CoWt7_1c)) |>
  ungroup()
PM_CWS_WM

PM_CWS_Summ <- PM_CoWtSel |> #create separate data frame for county-level averages and weights per year
  select(county, year, CoWt7_01Spl, CoWt7Raw, CoWt7_01, CoWt7_1c, CoAvYr, is20) |>
  group_by(county, year, CoWt7_01Spl, CoWt7Raw, CoWt7_01, CoWt7_1c, CoAvYr, is20) |>
  summarise() |>
  ungroup()


# ---- Plot County-level Weight RD ----

WtRD <- PM_CWS_Summ |> #plot all pop PM2.5 county averages per year in RD format
  ggplot(aes(x = CoWt7_01Spl, y = CoAvYr, group = year)) +
  geom_point(aes(color = as.factor(year)), alpha = PM_CWS_Summ$CoWt7_01, size = 2.3) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.85) +
  geom_segment(alpha = 0.009, linetype = "3333", color = "red", 
               aes(x = 1, xend = 2, y = filter(PM_CWS_WM, year == 2020)$WM, 
                   yend = filter(PM_CWS_WM, year == 2020)$WM)) +
  geom_segment(alpha = 0.009, linetype = "3333", color = "red", 
               aes(x = 0, xend = 1, y = filter(PM_CWS_WM, year == 2019)$WM, 
                   yend = filter(PM_CWS_WM, year == 2019)$WM, color = "WMline")) +
  #geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  #scale_linetype_manual(name = "Wtd. Mean", values = c("WMline" = "red")) +
  scale_x_continuous(limits = c(0, 2), expand = c(0.002, 0.002)) +
  scale_y_continuous(limits = c(0, 12.5), expand = c(0.0031, 0.0031)) +
  theme_bw() +
  scale_color_manual(values = c("2019" = "#bf401d", "2020" = "#3e39bf")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", margin=margin(0,0,12,0)), 
        legend.title = element_text(hjust = 0.5, face = "bold", size = 12),
        plot.margin = grid::unit(c(5,5,5,5), "mm"),
        legend.text = element_text(size = 11)) +
  guides(color = guide_legend(title = "Year")) +
  labs(x = "Weighting by Data Robustness \n(1.0 Represents Highest Weight and Year Cutoff)",
       y = expression("Mean PM2.5 (ug/m"^-3*") by County"),
       title = "Average PM2.5 Exposure Year Over Year - March-April 2019 vs. 2020")
WtRD
ggsave("figures/RD_CoWt_Scatterplot_19v20.png", width = 8, height = 6.5, units = "in", dpi = 500)


# ---- (Unused) Tract-level PM Weighting & Recode ----

PM_CTWt <- PM |> #calculate tract pops served and filter out tracts with 0 pop
  group_by(tract) |>
  mutate(CTPopSe = sum(unique(race_total_block_group))) |>
  filter(CTPopSe > 20)
PM_CTWt <- PM_CTWt |> #create sensor-level weights for data availability
  group_by(site_ID) |>
  mutate(SeWt = length(site_ID) / AvObs) |>
  ungroup()
PM_CTWt <- PM_CTWt |> #average sensor-level weights by tract, total tract population served, and create tract-level weights
  group_by(tract) |>
  mutate(CTSeWt = mean(SeWt),
         CTWt7Raw = length(tract) / CTPopSe^.9 * CTSeWt) |>
  ungroup() |>
  mutate(CTWt7_01 = (CTWt7Raw / max(CTWt7Raw))^.2, #standardize weight between >0 and ≤1, compress
         CTWt7_1c = ((CTWt7Raw / mean(CTWt7Raw))^.2) - 1, #standardize around 1 at avg then temporarily to 0, compress
         CTWt7_1c = ((CTWt7_1c / 2) + 1), #reduce range of weighting then center back around 1, compress
         CTWt7_1c = (1 - mean(CTWt7_1c)) + CTWt7_1c, #correct mean back to 1
         is20 = case_when(year == 2020 ~ 1, TRUE ~ 0)) #dummy variable, 1 for 2020, 0 for 2019

#filter for Mar 20 to April 30, 2019 and 2020 in new DF, recode 2020 weights to be above 1, filter for only PM2.5 parameter
PM_CTWtSel <- PM_CTWt |>
  filter(month == 3 & day > 19 & parameter == "PM2.5" | month == 4 & parameter == "PM2.5") |>
  mutate(CTWt7_01Spl = case_when(year == 2020 ~ 2-CTWt7_01, TRUE ~ CTWt7_01))
PM_CTWtSel <- PM_CTWtSel |> #generate tract-level yearly averages
  group_by(tract, year) |>
  mutate(CTAvYr = mean(sample_mean)) |>
  ungroup()

PM_CTS_WM <- PM_CTWtSel |> #Calculate weighted mean by year
  group_by(year) |>
  summarize(WM = weighted.mean(CTAvYr, CTWt7_1c)) |>
  ungroup()
PM_CTS_WM

PM_CTS_Mean <- PM_CTWtSel |> #Calculate normal mean by year
  group_by(year) |>
  summarize(mean = mean(CTAvYr)) |>
  ungroup()
PM_CTS_Mean

PM_CTS_Summ <- PM_CTWtSel |> #create separate data frame for tract-level averages and weights per year
  select(tract, year, CTWt7_01Spl, CTWt7Raw, CTWt7_01, CTWt7_1c, CTAvYr, is20) |>
  group_by(tract, year, CTWt7_01Spl, CTWt7Raw, CTWt7_01, CTWt7_1c, CTAvYr, is20) |>
  summarise() |>
  ungroup()


# ---- (Unused) Plot Tract-level Weight RD ----

#plot all pop PM2.5 tract averages per year in RD format
CTWtRD <- PM_CTS_Summ |> 
  ggplot(aes(x = CTWt7_01Spl, y = CTAvYr, group = year)) +
  geom_point(aes(color = as.factor(year)), alpha = PM_CTS_Summ$CTWt7_01, size = 2.3) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.85) +
  geom_segment(alpha = 0.009, linetype = "3333", color = "red", 
               aes(x = 1, xend = 2, y = filter(PM_CTS_WM, year == 2020)$WM, 
                   yend = filter(PM_CTS_WM, year == 2020)$WM)) +
  geom_segment(alpha = 0.009, linetype = "3333", color = "red", 
               aes(x = 0, xend = 1, y = filter(PM_CTS_WM, year == 2019)$WM, 
                   yend = filter(PM_CTS_WM, year == 2019)$WM, color = "WMline")) +
  #geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  #scale_linetype_manual(name = "Wtd. Mean", values = c("WMline" = "red")) +
  scale_x_continuous(limits = c(0, 2), expand = c(0.002, 0.002)) +
  scale_y_continuous(limits = c(0, 26), expand = c(0.0031, 0.0031)) +
  theme_bw() +
  scale_color_manual(values = c("2019" = "#bf401d", "2020" = "#3e39bf")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", margin=margin(0,0,12,0)), 
        legend.title = element_text(hjust = 0.5, face = "bold", size = 12),
        plot.margin = grid::unit(c(5,5,5,5), "mm"),
        legend.text = element_text(size = 11)) +
  guides(color = guide_legend(title = "Year")) +
  labs(x = "Weighting by Data Robustness \n(1.0 Represents Highest Weight and Year Cutoff)",
       y = expression("Mean PM2.5 (ug/m"^-3*") by Tract"),
       title = "Average PM2.5 Exposure Year Over Year - March-April 2019 vs. 2020")
CTWtRD
ggsave("figures/RD_CTWt_Scatterplot_19v20.png", width = 8, height = 6.5, units = "in", dpi = 500)


# ---- Plot County-level Weight as Binary/Discrete ----

#plot all pop PM2.5 county averages per year in discrete format
CoWtBi <- PM_CWS_Summ |> 
  ggplot(aes(x = is20, y = CoAvYr)) +
  geom_point(aes(color = as.factor(year)), alpha = 0.5, size = 2.3) +
  #geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.85) +
  geom_segment(alpha = 0.009, linetype = "3333", color = "red", 
               aes(x = .9, xend = 1.1, y = filter(PM_CWS_WM, year == 2020)$WM, 
                   yend = filter(PM_CWS_WM, year == 2020)$WM)) +
  geom_segment(alpha = 0.009, linetype = "3333", color = "red", 
               aes(x = -.1, xend = .1, y = filter(PM_CWS_WM, year == 2019)$WM, 
                   yend = filter(PM_CWS_WM, year == 2019)$WM, color = "WMline")) +
  #geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  #scale_linetype_manual(name = "Wtd. WM", values = c("WMline" = "red")) +
  scale_x_continuous(limits = c(-.25, 1.25), expand = c(0.002, 0.002)) +
  scale_y_continuous(limits = c(0, 12.5), expand = c(0.0031, 0.0031)) +
  theme_bw() +
  scale_color_manual(values = c("2019" = "#bf401d", "2020" = "#3e39bf")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", margin=margin(0,0,12,0)), 
        legend.title = element_text(hjust = 0.5, face = "bold", size = 12),
        plot.margin = grid::unit(c(5,5,5,5), "mm"),
        legend.text = element_text(size = 11)) +
  guides(color = guide_legend(title = "Year")) +
  labs(x = "Difference by Year",
       y = expression("WM PM2.5 (ug/m"^-3*") by County"),
       title = "Average PM2.5 Exposure Year Over Year \nMarch-April 2019 vs. 2020")
CoWtBi
ggsave("figures/County_Binary_Scatterplot_19v20.png", width = 8, height = 6.5, units = "in", dpi = 500)


# ---- Create Racial Designation Grouping Variable ----

#Filter PM long into date ranges, keep race parameters, generate % white column
PMlong_Sel <- PMlong_bins_all |>
  filter(month == 3 & day > 19 & parameter == "PM2.5" | month == 4 & parameter == "PM2.5") |>
  filter(EJ_name != "ln(Median Income)") |>
  mutate(pctWhite = race_white_block_group / race_total_block_group) |>
  mutate(raceCtg = case_when(pctWhite > .8 ~ 1,
                             pctWhite > .5 ~ 2,
                             (EJ_name == "Share Black (%)" & EJ_value > 80) ~ 3,
                             (EJ_name == "Share Black (%)" & EJ_value > 50) ~ 4,
                             (EJ_name == "Share Hispanic/Latinx (%)" & EJ_value > 80) ~ 5,
                             (EJ_name == "Share Hispanic/Latinx (%)" & EJ_value > 50) ~ 6,
                             (EJ_name == "Share Asian (%)" & EJ_value > 80) ~ 7,
                             (EJ_name == "Share Asian (%)" & EJ_value > 50) ~ 8,
                             TRUE ~ 0),
         raceCtgName = case_when(raceCtg == 0 ~ "Relatively Diverse (None Exceed 50%)",
                                 raceCtg == 1 ~ "Vast Majority White (>80%)",
                                 raceCtg == 2 ~ "Majority White (>50%)",
                                 raceCtg == 3 ~ "Vast Majority Black (>80%)",
                                 raceCtg == 4 ~ "Majority Black (>50%)",
                                 raceCtg == 5 ~ "Vast Majority Hispanic/Latinx (>80%)",
                                 raceCtg == 6 ~ "Majority Hispanic/Latinx (>50%)",
                                 raceCtg == 7 ~ "Vast Majority Asian (>80%)",
                                 raceCtg == 8 ~ "Majority Asian (>50%)"),
         isVast = ifelse(str_detect(raceCtgName, "Vast"), "Relatively Segregated", "Relatively Even Distribution"))

PMlong_Sel$raceCtg <- as.factor(PMlong_Sel$raceCtg)


# ---- Plot Racial Groups in Density format

#multi density chart of exposure by race and year
RLreorder <- levels(factor(PMlong_Sel$raceCtg))

PMlongDP <- PMlong_Sel |>
  #filter(year == 2020) |>
  mutate(raceCtgName = fct_reorder(raceCtgName, as.numeric(raceCtg), .desc = FALSE)) |>
  ggplot(aes(x = sample_mean, group = raceCtgName)) +
  stat_density(aes(color = raceCtgName, linetype = raceCtgName, alpha = raceCtgName, linewidth = raceCtgName), 
               geom = "line", position = "identity") +
  scale_x_continuous(trans = "sqrt", limits = c(0, 50)) +
  scale_y_continuous(limits = c(0, 0.58)) +
  theme_bw() +
  scale_color_manual(values = c("gray25", "darkorange3", "darkorange2", "brown3", 
                                 "royalblue3", "royalblue1", "springgreen4", "springgreen3"),
                     guide = guide_legend(title = "Racial Designation")) +
  scale_alpha_manual(values = c(1, .9, .8, .8, .9, .8, .9, .8),
                     guide = guide_legend(title = "Racial Designation")) +
  scale_linetype_manual(values = c("solid", "solid", "2222", "2222", "solid", "2222", "solid", "2222"),
                        guide = guide_legend(title = "Racial Designation")) +
  scale_linewidth_manual(values = c(1, .9, .7, .7, .9, .7, .9, .7),
                         guide = guide_legend(title = "Racial Designation")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", margin=margin(0,0,12,0)), 
        legend.title = element_text(hjust = 0.22, face = "bold", size = 13),
        plot.margin = grid::unit(c(5,5,5,5), "mm"),
        legend.text = element_text(size = 11.5)) +
  labs(x = expression("PM2.5 (ug/m"^-3*") by Sensor"),
       y = "Frequency",
       title = "Average PM2.5 Exposure by Racial Demographic Category \nMarch-April 2019 and 2020") +
  facet_grid(year ~ isVast, scales = "free_y")
PMlongDP
ggsave("figures/RaceDsgGroup_MargDensPlot_19v20+SegvDiv.png", width = 16, height = 6.5, units = "in", dpi = 400)


# ---- (Unused) Income Calculations & Reorganization ----

#calculate weighted mean in log median income by county
summary(PM_CoWtSel$medincome_block_group)

PM_CWS_Inc <- PM_CoWtSel |>
  filter(medincome_cbg != "NA") |>
  group_by(county, year) |>
  mutate(WMCoInc = weighted.mean(medincome_cbg, CoWt7_1c)) |>
  ungroup()
PM_CWS_Inc

#summarize income pre vs post by county
PM_Inc_Summ <- PM_CWS_Inc |>
  select(county, year, CoWt7_1c, WMCoInc, CoAvYr) |>
  group_by(county, year, CoWt7_1c, WMCoInc, CoAvYr) |>
  summarise() |>
  ungroup()
PM_Inc_Summ

#create income pre vs post scatter with fit line by county
#incplot <- PM_CWS_Inc |>
#  ggplot(aes(x = WMCoInc, y = CoAvYr, group = year)) +
#  geom_point(aes(color = as.factor(year)), size = 2.3) +
#  geom_smooth(method = "lm", formula = y ~ x, se = TRUE)
#incplot



# ---- Regression Table ----

PMlong_SelAll <- PMlong_bins_all |>
  filter(month == 3 & day > 19 & parameter == "PM2.5" | month == 4 & parameter == "PM2.5") |>
  mutate(EJ_value2 = EJ_value)

# PM2.5 March-April 2019 by Race + Income
PMlong_SelAll_income19 <- PMlong_SelAll |>
  filter(year == 2019 & EJ_name == 'ln(Median Income)') |> 
  lm(formula = sample_mean ~ EJ_value)
PM_lnincome <- summary(PMlong_SelAll_income19)

PMlong_SelAll_asian19 <- PMlong_SelAll |>
  filter(year == 2019 & EJ_name == 'Share Asian (%)') |> 
  lm(formula = sample_mean ~ EJ_value)
PM_asian<-summary(PMlong_SelAll_asian19)

PMlong_SelAll_black19 <- PMlong_SelAll |>
  filter(year == 2019 & EJ_name == 'Share Black (%)') |> 
  lm(formula = sample_mean ~ EJ_value)
PM_black<-summary(PMlong_SelAll_black19)

PMlong_SelAll_hisp19 <- PMlong_SelAll |>
  filter(year == 2019 & EJ_name == 'Share Hispanic/Latinx (%)') |> 
  lm(formula = sample_mean ~ EJ_value)
PM_hisp<-summary(PMlong_SelAll_hisp19)

# PM2.5 March-April 2020 by Race + Income
PMlong_SelAll_income20 <- PMlong_SelAll |>
  filter(year == 2020 & EJ_name == 'ln(Median Income)') |> 
  lm(formula = sample_mean ~ EJ_value2)
PM_lnincome <- summary(PMlong_SelAll_income20)

PMlong_SelAll_asian20 <- PMlong_SelAll |>
  filter(year == 2020 & EJ_name == 'Share Asian (%)') |> 
  lm(formula = sample_mean ~ EJ_value2)
PM_asian<-summary(PMlong_SelAll_asian20)

PMlong_SelAll_black20 <- PMlong_SelAll |>
  filter(year == 2020 & EJ_name == 'Share Black (%)') |> 
  lm(formula = sample_mean ~ EJ_value2)
PM_black<-summary(PMlong_SelAll_black20)

PMlong_SelAll_hisp20 <- PMlong_SelAll |>
  filter(year == 2020 & EJ_name == 'Share Hispanic/Latinx (%)') |> 
  lm(formula = sample_mean ~ EJ_value2)
PM_hisp<-summary(PMlong_SelAll_hisp20)

PM_fits <- stargazer(PMlong_SelAll_income19, PMlong_SelAll_asian19, PMlong_SelAll_black19, PMlong_SelAll_hisp19,
                         PMlong_SelAll_income20, PMlong_SelAll_asian20, PMlong_SelAll_black20, PMlong_SelAll_hisp20,
                         digits = 2,
                         title = "Naive OLS of Pre- and Post-COVID-19 Shutdown PM2.5 Exposure by Demographic Characteristics",
                         column.labels = c("Income (log)", "% Asian", "% Black", "% Hispanic/Latinx", 
                                           "Income (log)", "% Asian", "% Black", "% Hispanic/Latinx"),
                         dep.var.labels = c("PM2.5 Exposure - March-April 2019 and 2020"),
                         omit.stat = c("adj.rsq","rsq", "f", "ser"),
                         type = "html",
                         out = "figures/PM_fits.htm",
                         notes.align = "l",
                         notes.label = "Significance Codes:",
                         intercept.bottom = TRUE,
                         covariate.labels = c("2019 Coefficients", "2020 Coefficiencts", "Intercept"))


webshot::webshot("file:///Users/andreeanes/Downloads/AQRD/aqrd-final/figures/PM_fits.htm", 
                 file = "figures/PM_fits.png", zoom = 3, delay = 10)
image_crop(image_read("figures/PM_fits.png"),"2590x980") |>
  image_write("figures/PM_fits.png")
PM_fitsCopy <- image_read("file:///Users/andreeanes/Downloads/AQRD/aqrd-final/figures/PM_fits.png")
BlankWhite <- image_blank(width = 2590, height = 980, color = "white")
PM_fitsWhite <- image_composite(BlankWhite, PM_fitsCopy, offset = "+0+0")
image_write(PM_fitsWhite, path = "figures/PM_fits_white.png")


# ---- Summary Stats Table ----

#Based on PM_CoWtSel

# Daily Average PM2.5 (µg/cubic meter) (sample_mean),
# Average March-April PM2.5 by County and Year (µg/cubic meter) (CoAvYr),
# Block Group Population (race_total_block_group),
# Block Group White Population (race_white_block_group),
# Block Group Black Population (race_black_block_group),
# Block Group Asian Population (race_asian_block_group),
# Block Group Hispanic/Latinx Population (hisp_yes_block_group),
# Block Group Population Percent White (pctWhite),
# Block Group Population Percent Black (perc_black_cbg),
# Block Group Population Percent Asian (perc_asian_cbg),
# Block Group Population Percent Hispanix/Latinx (perc_hisp_cbg),
# Median Income by Census Block (medincome_block_group),
# (log) Median Income by Census Block (medincome_cbg)

#PM_SumVars <- c("sample_mean", "CoAvYr", "race_total_block_group", "race_white_block_group", "
#                race_black_block_group", "race_asian_block_group", "hisp_yes_block_group", "pctWhite", "
#                perc_black_cbg", "perc_asian_cbg", "perc_hisp_cbg", "medincome_block_group", "medincome_cbg")

PMSumDF <- data.frame(
    "Variable" = c("Daily Average PM2.5 (µg/cubic meter)",
                   "Average Mar.-Apr. PM2.5 by County (µg/cubic meter)",
                   "Block Group Population",
                   "Block Group White Population",
                   "Block Group Black Population",
                   "Block Group Asian Population",
                   "Block Group Hispanic/Latinx Population",
                   "Block Group Population Percent White",
                   "Block Group Population Percent Black",
                   "Block Group Population Percent Asian",
                   "Block Group Population Percent Hispanix/Latinx",
                   "Median Income by Census Block",
                   "(log) Median Income by Census Block"),
    "Observations" = c(length(PM_CoWtSel$sample_mean), length(unique(PM_CoWtSel$CoAvYr)), 
                       length(PM_CoWtSel$race_total_block_group), length(PM_CoWtSel$race_white_block_group), 
                       length(PM_CoWtSel$race_black_block_group), length(PM_CoWtSel$race_asian_block_group), 
                       length(PM_CoWtSel$hisp_yes_block_group), length(PM_CoWtSel$pctWhite), 
                       length(PM_CoWtSel$perc_black_cbg), length(PM_CoWtSel$perc_asian_cbg), 
                       length(PM_CoWtSel$perc_hisp_cbg), length(PM_CoWtSel$medincome_block_group), length(PM_CoWtSel$medincome_cbg)),
    "Mean" = c(mean(PM_CoWtSel$sample_mean, na.rm = TRUE), mean(unique(PM_CoWtSel$CoAvYr), na.rm = TRUE), 
               mean(PM_CoWtSel$race_total_block_group, na.rm = TRUE), mean(PM_CoWtSel$race_white_block_group, na.rm = TRUE), 
               mean(PM_CoWtSel$race_black_block_group, na.rm = TRUE), mean(PM_CoWtSel$race_asian_block_group, na.rm = TRUE), 
               mean(PM_CoWtSel$hisp_yes_block_group, na.rm = TRUE), mean(PM_CoWtSel$pctWhite, na.rm = TRUE), 
               mean(PM_CoWtSel$perc_black_cbg, na.rm = TRUE), mean(PM_CoWtSel$perc_asian_cbg, na.rm = TRUE), 
               mean(PM_CoWtSel$perc_hisp_cbg, na.rm = TRUE), mean(PM_CoWtSel$medincome_block_group, na.rm = TRUE), 
               mean(PM_CoWtSel$medincome_cbg, na.rm = TRUE)),
    "Standard Deviation" = c(sd(PM_CoWtSel$sample_mean, na.rm = TRUE), sd(unique(PM_CoWtSel$CoAvYr), na.rm = TRUE), 
                             sd(PM_CoWtSel$race_total_block_group, na.rm = TRUE), sd(PM_CoWtSel$race_white_block_group, na.rm = TRUE), 
                             sd(PM_CoWtSel$race_black_block_group, na.rm = TRUE), sd(PM_CoWtSel$race_asian_block_group, na.rm = TRUE), 
                             sd(PM_CoWtSel$hisp_yes_block_group, na.rm = TRUE), sd(PM_CoWtSel$pctWhite, na.rm = TRUE), 
                             sd(PM_CoWtSel$perc_black_cbg, na.rm = TRUE), sd(PM_CoWtSel$perc_asian_cbg, na.rm = TRUE), 
                             sd(PM_CoWtSel$perc_hisp_cbg, na.rm = TRUE), sd(PM_CoWtSel$medincome_block_group, na.rm = TRUE), 
                             sd(PM_CoWtSel$medincome_cbg, na.rm = TRUE)))

PMSumDF <- format(PMSumDF, scientific = FALSE, digits = 1, big.mark = ",")

PM_SumStats <- PMSumDF |>
  gt() |>
  cols_align(align = "left", columns = Variable) |>
  tab_header(title = md("**Summary Statistics**")) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) |>
  cols_label(Standard.Deviation = "Standard Deviation") |>
  tab_options(table_body.hlines.color = "black",
              column_labels.border.top.color = "black",
              column_labels.border.bottom.color = "black",
              table_body.border.bottom.color = "black",
              table.border.top.color = "black")
PM_SumStats

gtsave(PM_SumStats, filename = "figures/PM_SumStatsTable.png")