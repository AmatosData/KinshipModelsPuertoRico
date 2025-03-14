###################################################
#     Code for Kinship Analysis for Puerto Rico
#     Migration and Kinship in Puerto Rico 2022
#     Prepared by: Amílcar Matos Moreno
#
###################################################
getwd()
setwd("/Users/amilcar/desktop/Max planck Insitute Kinship Models/Data")

library(haven) 
# Import data from American Community Survey to estimate age and sex adjusted migrations rates from Puerto Rico to USA
mydata <- read_dta("ACS_USA_2017_2022.dta")
# Import data for return migration to estimate age and sex adjusted return migration rates to PR. 
retmig_data <- read_dta("IPUMS_PR_2022_returnmig.dta")


######Clean Data frame for migration estimates
#Select variables of interest
mig_data22<- subset(mydata, select=c(perwt,sex,age,migplac1))
retmig_data22<-subset(retmig_data, select=c(perwt,sex,age,migplac1, bpld))

#Select only women migrating from Puerto Rico to USA
mig_data22_PR<-mig_data22[ which(mig_data22$migplac1==110 & mig_data22$sex== 2), ]

# Select females born in Puerto Rico migrating back from any place to PR. 
retmig_data22<- retmig_data22[ which(retmig_data22$migplac1!=110 & retmig_data22$migplac1!=0 & retmig_data22$sex== 2 & retmig_data22$bpld==11000), ]

# Aggregate observations by Age
mig_data22_PR_final <- mig_data22_PR %>% group_by(age) %>%  summarise(count=sum(perwt))
head(mig_data22_PR_final)

retmig_data22_PR_final <- retmig_data22 %>% group_by(age) %>%  summarise(count=sum(perwt))
head(retmig_data22_PR_final)

# Add years "0" 96" "97" "98" "99" "100" "101" to have the same dimensions as other matrixes. 
newdata<-data.frame(c(0,92,96, 97,98,99,100), c(0,0,0,0,0,0,0))
names(newdata)<-c("age", "count")

newretmig_data<-data.frame(c(0,83,84,89,93,95,97,98,99,100), c(0,0,0,0,0,0,0,0,0,0))
names(newretmig_data)<-c("age", "count")

mig_data21_PR_Final<- rbind(mig_data22_PR_final, newdata)

retmig_data22_Final<- rbind(retmig_data22_PR_final, newretmig_data)

library(dplyr)
# Sort data frame by age of people. 
mig_data21_sorted<- print(arrange(mig_data21_PR_Final, age))

retmig_data21_sorted<- print(arrange(retmig_data22_Final, age))


# Upload Total population of females for 2022 in PR.
tot_pop_PR<- read.csv("tot_pop_fem_PR.csv")
tot_pop_PR_2021<- round(tot_pop_PR$X2021 *1000) # Rounding and multiplying decimals into thousands. 

# Add total pop in same data frame as migration to calculate age-specific migration rate
mig_data21_sorted$totpop <- tot_pop_PR_2021

retmig_data21_sorted$totpop <- tot_pop_PR_2021



#Calculate age_speciic migration rate for 2021
mig_rate_PRfem_2021 <- mig_data21_sorted %>% summarise(mig_rate=count/totpop)
#Calculate age_speciic return migration rate for 2021
retmig_rate_PRfem_2021 <- retmig_data21_sorted %>% summarise(mig_rate=count/totpop)

# Add age_specific migration rates to the same data frame. 
mig_data21_sorted$mig_rate<- mig_rate_PRfem_2021

mig_data21_sorted$retmig_rate<- retmig_rate_PRfem_2021
names(mig_data21_sorted)<-c("age", "count", "totpop","mig_rate","retmig_rate")

attach(mig_data21_sorted)
plot(age, mig_rate$mig_rate)
plot(age, retmig_rate$mig_rate)
detach()

save(mig_data21_sorted, file="migdataFEM21PR.rdata")

plot()

# Codigo de ivan para calcular la probabilidad de quedarse en PR de cada persona. 
# mig data
# load("amilcar/migdataFEM21PR.rdata")
# ux <- cbind(mig_data21_sorted$mig_rate)
# convert migr rate to probability
# ux$q <- 1 - exp(-ux$mig_rate)
# probability to stay in PR
# ux$p <- 1 - ux$q
# cohort stay
# ux$lx <- cumprod(ux$p)
# plot
# plot(ux$lx, main = "l(x) stay in PR", ylab = "l(x)", xlab = "Age")
# abline(v = 65)  


####################### Migration data ready #################################
library(DemoKin)
library(dplyr)
library(tidyr)
library(ggplot2)


source("UNWPP_data.R")
#select country, year and sex
data_PR_2021 <- UNWPP_data(country = "Puerto Rico",
                           start_year =  2021,
                           end_year = 2021,
                           sex = "Female")

save(data_PR_2021, file="data_PR_2021.rdata")

country_fert <- data_PR_2021 %>%
  select(age, year, fx) %>%
  pivot_wider(names_from = year, values_from = fx) %>%
  select(-age) %>%
  as.matrix()
#Survival
country_surv <- data_PR_2021 %>%
  select(age, year, px) %>%
  pivot_wider(names_from = year, values_from = px) %>%
  select(-age) %>%
  as.matrix()



country_mig<-subset(mig_data21_sorted, select=c(mig_rate,retmig_rate))

# Rearrange dataframe to have the stages in order of ocurrence  1, 2, 3.
country_mig<-relocate(country_mig, mig_rate, retmig_rate)

country_surv_3<-as.data.frame(country_surv)
country_surv_3$px1 <- country_surv_3$"2021"
names(country_surv_3)<-c("px1","px2")

save(country_mig, file="Mig_returnmig_stages.rdata")


# Translate data to a list format and use those as inputs in DemoKin function.
PR_2021_surv <- data_PR_2021 %>%
  select(px1 = px) 

####################################################################
### Build the multistage model
###################################################################
source("kin_multi_stage.R")
# create lists
lists_PR2021 <- make_mulstistate_parity_matrices(f_parity = country_mig, p_parity = country_surv_3)

PR_2021_kin_multi <-
  kin_multi_stage(
    U = lists_PR2021[["U"]],
    f = lists_PR2021[["F."]],
    D = lists_PR2021[["D"]],
    H = lists_PR2021[["H"]])

########      ######      ######      ######    ######
# Distribution of Daughters and Migration in 2021 for Puerto Rico. 
PR_2021_kin_multi %>%
  filter(kin %in% "d") %>%
  mutate(migration = as.integer(stage_kin)) %>%
  group_by(age_focal, kin, migration) %>%
  summarise(count= sum(living)) %>%
  DemoKin::rename_kin() %>% print(n=220) %>%
  ggplot() +
  geom_bar(aes(x=age_focal, y = count, fill=migration), stat = "identity") +
  labs(y = "Daughters") +
  theme_bw() +
  facet_wrap(~kin, nrow = 2)

# Distribution of Living kin and migration by kin type
PR_2021_kin_multi %>%
  filter(kin %in% c("d", "m", "ys", "os")) %>%
  mutate(migration = as.integer(stage_kin)) %>%
  group_by(age_focal, kin, migration) %>%
  summarise(count= sum(living)) %>%
  DemoKin::rename_kin() %>%
  ggplot() +
  geom_bar(aes(x=age_focal, y = count, fill=migration), stat = "identity") +
  labs(y = "Daughters") +
  theme_bw() +
  facet_wrap(~kin, nrow = 2) + ggtitle("Expected Migration by Kin type")

# Distribution of Younger Sister and Migration in 2021 for Puerto Rico. 
PR_2021_kin_multi %>%
  filter(kin %in% c("ys","os","nos", "nys")) %>%
  mutate(migration = stage_kin) %>%
  group_by(age_focal, kin, migration) %>%
  summarise(count= sum(living)) %>%
  DemoKin::rename_kin() %>%
  ggplot() +
  geom_bar(aes(x=age_focal, y = count, fill=migration), stat = "identity") +
  labs(y = "Sisters") +
  theme_bw() +
  facet_wrap(~kin, nrow = 2)


#Mean age of daughters by migration status for a focal 65 years old. 
PR_2021_kin_multi %>% filter(age_focal==65, kin=="d") %>% 
  group_by(as.factor(stage_kin)) %>%
  summarise(mean_age = sum(living*age_kin)/sum(living))



PR_2021_kin_multi %>% 
  plot_diagram(kin_total, rounding=3)







#######################################     ############################################

# Lets see the Keytz diagram for Puerto Rico in 2021. 
pr_2021 <- 
  kin(p = country_surv, f = country_fert, time_invariant = TRUE)

pr_2021$kin_summary %>% 
  filter(age_focal == 65) %>% 
  select(kin, count = count_living) %>% 
  plot_diagram(rounding = 2)

# Plot Kinship structure
PR_2021_kin_multi$kin_summary %>%
  filter(age_focal == 65) %>%
  select(kin, count = count_living) %>%
  plot_diagram(rounding = 2)




###################################################
# Two sex model without migration for Puerto Rico
##################################################

source("UNWPP_data.R")
#select country, year and sex
data_PR_2021_female <- UNWPP_data(country = "Puerto Rico",
                           start_year =  1950,
                           end_year = 2021,
                           sex = "Female")

source("UNWPP_data.R")
#select country, year and sex
data_PR_2021_male <- UNWPP_data(country = "Puerto Rico",
                           start_year =  1950,
                           end_year = 2021,
                           sex = "Male")


# Reshaping data to matrix
PR_asfr_females <- data_PR_2021_female %>%
  select(age, year, fx) %>%
  pivot_wider(names_from = year, values_from = fx) %>%
  select(-age) %>%
  as.matrix()
PR_surv_females <- data_PR_2021_female %>%
  select(age, year, px) %>%
  pivot_wider(names_from = year, values_from = px) %>%
  select(-age) %>%
  as.matrix()
PR_surv_males <- data_PR_2021_male %>%
  select(age, year, px) %>%
  pivot_wider(names_from = year, values_from = px) %>%
  select(-age) %>%
  as.matrix()

##### Male fertility assumption (male mean fatherhood age (29 yd; female mean childbearing age 25yd)
PR_asfr_male<- PR_asfr_females
Male_test<- matrix(data = 0, nrow = 4, ncol = 72)
PR_asfr_male<-rbind(Male_test,PR_asfr_females)
PR_asfr_male<- head(PR_asfr_male, -4)





# E1.1 In a time invariant model: how many living children can 
# a woman expect to have at age 65 in 2000, 2010, and in 2021,
# and what are their mean ages?

PR_kin_2sex <- kin2sex(
  pf = PR_surv_females, pm = PR_surv_males, ff = PR_asfr_females,
  fm = PR_asfr_male, time_invariant = TRUE, sex_focal = "f")


kin_out_1950 <- kin2sex(
  pf = PR_surv_females[,"1950"],pm = PR_surv_males[,"1950"],
  ff = PR_asfr_females[,"1950"],fm = PR_asfr_male[,"1950"],
  time_invariant = TRUE,
  sex_focal = "f")

kin_out_2000 <- kin2sex(
  pf = PR_surv_females[,"2000"],pm = PR_surv_males[,"2000"],
  ff = PR_asfr_females[,"2000"],fm = PR_asfr_male[,"2000"],
  time_invariant = TRUE,
  sex_focal = "f")

kin_out_2010 <- kin2sex(
  pf = PR_surv_females[,"2010"],pm = PR_surv_males[,"2010"],
  ff = PR_asfr_females[,"2010"],fm = PR_asfr_male[,"2010"],
  time_invariant = TRUE)

kin_out_2021 <- kin2sex(
  pf = PR_surv_females[,"2021"], pm = PR_surv_males[,"2021"], 
  ff = PR_asfr_females[,"2021"],
  fm = PR_asfr_male[,"2021"],
  time_invariant = TRUE,
  sex_focal = "f")


bind_rows(
  kin_out_2000$kin_summary %>%
    filter(age_focal==65, kin == "d") %>%
    select(count_living, mean_age, kin) %>%
    mutate(year = 2000),
  kin_out_2010$kin_summary %>%
    filter(age_focal==65, kin == "d") %>%
    select( count_living, mean_age, kin) %>%
    mutate(year = 2010),
  kin_out_2021$kin_summary %>%
    filter(age_focal==65, kin == "d") %>%
    select(count_living, mean_age, kin) %>%
    mutate(year = 2021)
)

# Conclusion: The burden of caregiving is much higher in these times compared to Puerto Rico
# in the 1950s. This is not a bad things, is just a change is longevity, however, 
# is in important to account for the increase number of ascendent relatives in Puerto Rico 2021. 

# Look at the number of living kin in the time invariant model
# lets group aunts and siblings and visualize the number of living kin
PR_kin_2sex <- kin2sex(
  pf = PR_surv_females[,"2021"], pm = PR_surv_males[,"2021"], ff = PR_asfr_females[,"2021"],
  fm = PR_asfr_male[,"2021"], time_invariant = TRUE, sex_focal = "f")

kin_out <- PR_kin_2sex$kin_summary %>%
  mutate(kin = case_when(kin %in% c("os", "ys") ~ "s",
                         kin %in% c("ya", "oa") ~ "a",
                         T ~ kin)) %>%
  filter(kin %in% c("d", "m", "gm", "ggm", "s", "a"))

kin_out %>%
  summarise(count=sum(count_living), .by = c(kin, age_focal, sex_kin)) %>%
  ggplot(aes(age_focal, count, fill=sex_kin))+
  geom_area()+
  theme_bw() +
  ggtitle("Number of total living kin in Puerto Rico 2021 by Focal Age") +
  facet_wrap(~kin)

###### Conclusion: Children and Siblings still the greatest assets for family support in 
##### Puerto Rico. Within these kin types, female presence is greater than male presence.

############
# Plot the time variant and time invariant model. 
############
PR_kin_2sex <- kin2sex(
  pf = PR_surv_females[,"1950"], pm = PR_surv_males[,"1950"], ff = PR_asfr_females[,"1950"],
  fm = PR_asfr_male[,"1950"], time_invariant = TRUE, sex_focal = "f")

PR_kin_2sex_variant <- kin2sex(
  pf = PR_surv_females, pm = PR_surv_males, ff = PR_asfr_females,
  fm = PR_asfr_male, time_invariant = FALSE, sex_focal = "f")


PR_kin_2sex_variant$kin_summary %>%
  filter(cohort == 1950) %>% mutate(type = "variant") %>%
  bind_rows(PR_kin_2sex$kin_summary %>% mutate(type = "invariant")) %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s",
                         kin %in% c("ya", "oa") ~ "a",
                         T ~ kin)) %>%
  filter(kin %in% c("d", "m", "gm", "ggm", "s", "a")) %>%
  group_by(type, kin, age_focal, sex_kin) %>%
  summarise(count=sum(count_living)) %>%
  ggplot(aes(age_focal, count, linetype=type))+
  geom_line()+ theme_bw() +
  ggtitle("Time invariant and variant models for Puerto Rico 1950-2021")+
  facet_grid(cols = vars(kin), rows=vars(sex_kin), scales = "free")

##############
# Keytz Diagram for Puerto Rico 2000, 2010, 2021
#############

PR_kin_2sex_variant_65 <- kin2sex(
  pf = PR_surv_females, pm = PR_surv_males, ff = PR_asfr_females,
  fm = PR_asfr_male, time_invariant = FALSE, sex_focal = "f")

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(age_focal == 65, year=="2000") %>% 
  select(kin, count = count_living) %>% 
  plot_diagram(rounding = 2)

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(age_focal == 65, year=="2010") %>% 
  select(kin, count = count_living) %>% 
  plot_diagram(rounding = 2)

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(age_focal == 65, year=="2021") %>% 
  select(kin, count = count_living) %>% 
  plot_diagram(rounding = 2)


 ############
####
####  Data for table 1: Number of living close kin by decade and age. 
####
#############

### Average number of living kin for an older adults (65-100) in Puerto Rico. 
PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2010", kin %in% c("d", "m", "ys","os")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2010_65plus.csv")

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2010", kin %in% c("d")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2010_65plus_children.csv")

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2010", kin %in% c("os","ys")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2010_65plus_siblings.csv")

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2010", kin %in% c("m")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2010_65plus_parents.csv")

# 2000
PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2000", kin %in% c("d", "m", "ys","os")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2000_65plus.csv")

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2000", kin %in% c("d")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2000_65plus_children.csv")

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2000", kin %in% c("os","ys")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2000_65plus_siblings.csv")

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2000", kin %in% c("m")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2000_65plus_parents.csv")

# 2021
PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2021", kin %in% c("d", "m", "ys","os")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2021_65plus.csv")

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2021", kin %in% c("d")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2021_65plus_children.csv")

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2021", kin %in% c("os","ys")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2021_65plus_siblings.csv")

PR_kin_2sex_variant_65$kin_summary %>% 
  filter(year=="2021", kin %in% c("m")) %>% 
  summarise(count=sum(count_living), .by = c(age_focal, year)) %>%
  print(n=150) %>% write.csv("kin_estimate_2021_65plus_parents.csv")



## Number of living kin for a person 55-85 years old in Puerto Rico 2000, 2010, 2021
PR_kin_2sex_variant_65$kin_summary %>% 
  filter(age_focal %in% c(55,65,75,85), year %in% c("2000", "2010","2021"), kin %in% c("d", "m", "ys","os")) %>% 
  summarise(count=sum(count_living), .by = c(kin, age_focal, year)) %>% print(n=150)


## Number of living close kin for a person 55 to 85 in Puerto Rico 2000, 2010, 2021. 
PR_kin_2sex_variant_65$kin_summary %>% 
  filter(age_focal %in% c(55,65,75,85), year %in% c("2000","2010","2021")) %>% 
  summarise(total_count=sum(count_living[kin=="d"], na.rm=T) + sum(count_living[kin=="m"], na.rm=T) + sum(count_living[kin=="os"], na.rm=T) + sum(count_living[kin=="ys"], na.rm=T), .by= c(age_focal, year)) %>% print(n=150)


## Mean age of kin for a focal 55 years old.
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(age_focal==55, year %in% c("2000","2010","2021"), kin %in% c("d","m","s")) %>%
  select(count_living, mean_age, kin, year) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year)) %>% print(n=100)


## Mean age of kin for a focal 65 years old. 
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(age_focal==65, year %in% c("2000","2010","2021"), kin %in% c("d","m","s")) %>%
  select(count_living, mean_age, kin, year) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year)) %>% print(n=100)

## Mean age of kin for a focal 75 years old.
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(age_focal==75, year %in% c("2000","2010","2021"), kin %in% c("d","m","s")) %>%
  select(count_living, mean_age, kin, year) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year)) %>% print(n=100)

## Mean age of kin for a focal 85 years old.
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(age_focal==85, year %in% c("2000","2010","2021"), kin %in% c("d","m","s")) %>%
  select(count_living, mean_age, kin, year) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year)) %>% print(n=100)

## Mean age of kin for a focal 65plus years old (Children 2000). 
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(year=="2000", kin %in% c("d")) %>%
  select(count_living, mean_age, kin, year, age_focal) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year, age_focal)) %>% print(n=150)

## Mean age of kin for a focal 65 plus years old (Siblings 2000). 
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(year=="2000", kin %in% c("s")) %>%
  select(count_living, mean_age, kin, year, age_focal) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year, age_focal)) %>% print(n=150)

## Mean age of kin for a focal 65 plus years old (Parents 2000). 
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(year=="2000", kin %in% c("m")) %>%
  select(count_living, mean_age, kin, year, age_focal) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year, age_focal)) %>% print(n=150)

## Mean age of kin for a focal 65plus years old (Children 2010). 
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(year=="2010", kin %in% c("d")) %>%
  select(count_living, mean_age, kin, year, age_focal) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year, age_focal)) %>% print(n=150)

## Mean age of kin for a focal 65 plus years old (Siblings 2010). 
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(year=="2010", kin %in% c("s")) %>%
  select(count_living, mean_age, kin, year, age_focal) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year, age_focal)) %>% print(n=150)

## Mean age of kin for a focal 65 plus years old (Parents 2010). 
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(year=="2010", kin %in% c("m")) %>%
  select(count_living, mean_age, kin, year, age_focal) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year, age_focal)) %>% print(n=150)


## Mean age of kin for a focal 65plus years old (Children 2021). 
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(year=="2021", kin %in% c("d")) %>%
  select(count_living, mean_age, kin, year, age_focal) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year, age_focal)) %>% print(n=150)

## Mean age of kin for a focal 65 plus years old (Siblings 2021). 
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(year=="2021", kin %in% c("s")) %>%
  select(count_living, mean_age, kin, year, age_focal) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year, age_focal)) %>% print(n=150)

## Mean age of kin for a focal 65 plus years old (Parents 2021). 
PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "s", T ~ kin)) %>%
  filter(year=="2021", kin %in% c("m")) %>%
  select(count_living, mean_age, kin, year, age_focal) %>%
  summarise(mean_age=sum(mean_age/length(mean_age), na.rm=T), .by = c(kin, year, age_focal)) %>% print(n=150)


############
####
###
####
#############





##### Kin Ratio
kin_ratio_PR <- PR_kin_2sex_variant_65$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "si", T ~ kin)) %>%
  filter(kin %in% c("d", "m", "gm", "si"))

## Expected Number of Living close kin for 2021
kin_ratio_PR %>%
  filter(year=="2021") %>%
  summarise(count=sum(count_living), .by = c(kin, age_focal, sex_kin, year)) %>%
  ggplot(aes(age_focal, count, fill=sex_kin))+ ggtitle("Expected # of Living Kin") +
  geom_area()+
  theme_bw() +
  facet_wrap(~kin)

kin_ratio_PR %>%
  filter(year=="2021") %>% 
  summarise(count_living = sum(count_living), .by = c(age_focal, kin, sex_kin)) %>% 
  select(age_focal, kin, sex_kin, count_living) %>% 
  pivot_wider(names_from = sex_kin, values_from = count_living) %>% 
  mutate(sex_ratio = m/f) %>% print(n=510)
  


###############################################
## Figure 1 code 
# Expected number of total kin by kintype
kin_ratio_PR %>%
  filter(year %in% c("2000","2021")) %>%
  summarise(count=sum(count_living), .by = c(kin, age_focal, year)) %>%
  ggplot(aes(age_focal, count, fill=kin))+
  geom_area()+
  geom_area(colour = "black")+
  xlab("Age of Focal")+
  ylab("Expected Number of Close Living Kin")+
  scale_fill_discrete(name = "Kin Type", labels = c("Children", "Grandparents", "Parents", "Siblings"))+
  theme_bw()+
  theme(legend.position = "bottom") +
  facet_wrap(~year, nrow=1)

z <- (kin_ratio_PR %>%
  filter(year %in% c("2021")) %>%
  summarise(count=sum(count_living), .by = c(kin, age_focal, year)) %>%
  ggplot(aes(age_focal, count, fill=kin))+ ggtitle("Expected Number of Close Living Kin in Puerto Rico in 2021") +
  geom_area()+
  geom_area(colour = "black")+
  ylab("Expected Number of Close Living Kin")+
  scale_fill_discrete(name = "Kin Type", labels = c("Children", "Grandparents", "Parents", "Siblings"))+
  theme_bw()+
  theme(legend.position = "top", axis.title.x = element_blank(),
  axis.text.x = element_blank() ) )


# For figure 2
w <-(kin_ratio_PR %>%
       filter(year=="2021") %>% 
       summarise(count_living = sum(count_living), .by = c(age_focal, kin, sex_kin)) %>% 
       select(age_focal, kin, sex_kin, count_living) %>% 
       pivot_wider(names_from = sex_kin, values_from = count_living) %>% 
       mutate(sex_ratio = m/f) %>%
       ggplot(aes(age_focal, sex_ratio, col=kin)) +
       geom_line()+
       ggtitle("Sex ratio by Close Living Kin for 2021") +
       xlab("Age of Focal Person")+
       ylab("Sex Ratio of Close Kin (m/f)")+
       scale_fill_discrete(name = "Kin Type")+
       theme_bw() +
  theme(legend.position = "none") )

figure_2<- ggarrange(z, w, nrow=2, ncol=1, labels = c("A","B"))



################################################
# Total number of expected living close kin
kin_ratio_PR %>% 
  filter(year %in% c("2000","2010","2021")) %>%
  summarise(total_count=sum(count_living[kin=="d"], na.rm=T) + sum(count_living[kin=="m"], na.rm=T) + sum(count_living[kin=="si"], na.rm=T), .by= c(age_focal, year)) %>% print(n=150)

#Total count of children for every focal age 
kin_ratio_PR %>% 
  filter(year %in% c("2000","2010","2021")) %>%
  summarise(total_count=sum(count_living[kin=="d"], na.rm=T), .by= c(age_focal, sex_kin, year)) %>% print(n=150)

#Total count of siblings for every focal 
kin_ratio_PR %>% 
  filter(year %in% c("2000","2010","2021")) %>%
  summarise(total_count=sum(count_living[kin=="si"], na.rm=T), .by= c(age_focal, sex_kin, year)) %>% print(n=150)







# E1.2 Compare ‘kin sex ratios’ of grandparents (gm), parents (m), daughters (d) and siblings (m and n) in 
# a time-variant framework for the cohort 1956, at each of Focal.

PR_kin_2sex_1956 <- kin2sex(
  pf = PR_surv_females, pm = PR_surv_males, ff = PR_asfr_females,
  fm = PR_asfr_male, time_invariant = FALSE, sex_focal = "f", output_cohort = 1956)

PR_kin_2sex_1956$kin_summary

kin_ratio <- PR_kin_2sex_1956$kin_summary %>%
  mutate(kin = case_when(kin %in% c("ys", "os") ~ "si", T ~ kin)) %>%
  filter(kin %in% c("d", "m", "gm", "si"))

kin_ratio %>%
  summarise(count=sum(count_living), .by = c(kin, age_focal, sex_kin)) %>%
  ggplot(aes(age_focal, count, fill=sex_kin))+ ggtitle("Expected # of Living Kin for People born in Puerto Rico's 1956") +
  geom_area()+
  theme_bw() +
  facet_wrap(~kin)

# Total number of expected living close kin
kin_ratio %>% summarise(total_count=sum(count_living[kin=="d"], na.rm=T) + sum(count_living[kin=="m"], na.rm=T) + sum(count_living[kin=="d"], na.rm=T), .by= c(age_focal)) %>% print(n=150)


#Total count of children for every focal age 
kin_ratio %>% summarise(total_count=sum(count_living[kin=="d"], na.rm=T), .by= c(age_focal, sex_kin)) %>% print(n=150)

#Total count of siblings for every focal 
kin_ratio %>% summarise(total_count=sum(count_living[kin=="si"], na.rm=T), .by= c(age_focal, sex_kin)) %>% print(n=150)



# Sex ratio of kin: Will have a (number) of kin by (this number) of other kin. Grandfathers for every grandmother, etc. 
kin_ratio %>%
  group_by(kin, age_focal) %>%
  summarise(sex_ratio = sum(count_living[sex_kin=="m"], na.rm=T)/sum(count_living[sex_kin=="f"], na.rm=T)) %>%
  ggplot(aes(age_focal, sex_ratio))+ ggtitle("Sex ratio by Kin Type") +
  geom_line(title())+
  theme_bw() +
  facet_wrap(~kin, scales = "free")

#Mean age
kin_ratio_PR %>% filter(age_focal==65) %>% 
  group_by(kin) %>%
  summarise(mean_age = sum(count_living*mean_age)/sum(count_living))


# ages
kin_ratio_PR %>%
  filter(age_focal == 65, kin != "focal") %>%
  ggplot(aes(mean_age, count_living)) +
  geom_line() +
  facet_wrap(~kin, scales = "free_y")





