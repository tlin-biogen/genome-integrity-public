library(tidyverse)
library(kableExtra)
#############################################################################
# plotting function & color blindness friendly color palette
cols<-c("#767676" ,"#6388b4","#55ad89" , "#ffae34", "#8cc2ca" , "#ef6f6a","#c3bc3f" ,"#bb7693" )

#plotting function
aav_GI_plot <- function(dataset, x, y) {
  ggplot(dataset, aes({{x}}, {{y}},color=Sample_num)) + 
    geom_point() +  
    geom_line(aes(group=Sample_num),stat = "summary", fun = mean, linewidth=1) +
    geom_hline(aes(yintercept=expected,color=Sample_num),linetype='dashed',linewidth=1) +
    scale_colour_manual(values=cols,name="% Integrity sample",
                        labels=c("84%","63%","42%","21%","10%","5%","0%"))+
    scale_x_log10()+
    theme_classic(#base_size = 14
    )+
    ylim(-3.5,101)+
    xlab("Concentration (copies/uL)")+
    ylab("Full genome (%)")
}
#############################################################################

#############################################################################
# clean and sort data, compute genome integrity with BioRad linkage (avg) model
#############################################################################
# read in raw ddPCR files from a folder of AAV integrity experiments
setwd("~/R/genome_integrity/aav_files")
raw_file_list <- list.files(path = "~/R/genome_integrity/aav_files",pattern="^XZ")
aav_raw_compiled<-raw_file_list %>% read_csv(id = "file")

#select relevant cols of data
aav_raw_compiled <- select(aav_raw_compiled, file, Sample, Target, Concentration, Linkage)

# keep beginning of file name as Experiment ID
aav_raw_compiled <- aav_raw_compiled %>% separate_wider_delim(file, "_raw", names = c("Experiment", NA))

# drop dilutions from exp -131 that don't match the other experiments
aav_raw_compiled <- aav_raw_compiled %>% filter(!grepl("_(01|02|03|04|05|10|11|12)$", Sample))

# make Sample names match across all experiments
aav_raw_compiled$Sample <- gsub("^S0", "S", aav_raw_compiled$Sample)

#samples that starts with S only (drop neg and pos controls)
aav_raw_compiled <- aav_raw_compiled[str_detect(aav_raw_compiled$Sample, "S"), ]
#exclude not S1-S7
aav_raw_compiled <- aav_raw_compiled %>% filter(!grepl("^S(8|9|10|11|12)", Sample))

aav_raw_compiled<-separate(aav_raw_compiled,Sample,into=c("Sample_num","Replicate"),sep = "_")

# update target names to match manuscript
aav_raw_compiled <- aav_raw_compiled %>%
  mutate(Target = str_replace_all(Target, c("1" = "prom", "2" = "pA")))

# filter into separate dfs by target (since droplet info is in duplicate)
aav_prom<- aav_raw_compiled %>% filter(Target=="prom") 
aav_pA<- aav_raw_compiled %>% filter(Target=="pA") 

# expected integrity values for each exp
aav_prom$expected <- rep(as.numeric(c("84","63","42","21","10","5","0")), each=4,6)

# calculate %integrity with: BioRad %Linkage calculation equation 1 from Pranter 2023 S1 Text
aav_prom$BR1 <- aav_prom$Linkage/((aav_prom$Concentration+aav_pA$Concentration)/2)*100

# order data by experiment & sample number to accurately label theoretical copies for experiments
aav_prom <- aav_prom %>% arrange(Experiment,Sample_num)
aav_prom$copies<-rep(as.numeric(c("498","249","125","62")),42)

#sep expected to create a label
exp_lab = paste(aav_prom$expected, "%", sep = "")

BR1_aav<-aav_GI_plot(aav_prom,copies, BR1) + theme(legend.position = "none") +   
  geom_text(aes(650,expected,label = exp_lab, vjust = -0.5),show.legend = FALSE)

############################################################################################################
# clean and sort analysed data from Poisson-multinomial Shiny app, plot calculated integrity
############################################################################################################

# read Shiny app *.csv files
app_file_list <- list.files(path = "~/R/genome_integrity/aav_files",pattern="^0")
aav_app<-app_file_list %>% read_csv(id = "file")

#select relevant cols of data
aav_app <- select(aav_app,file, SAMPLE, TARGET,ESTIMATED_PERCENTAGE)
# keep beginning of file name as Experiment ID
aav_app<- aav_app %>% separate_wider_delim(file, "_appdata", names = c("Experiment", NA))
# adjust Experiment name to match the raw data files
aav_app <- aav_app %>%mutate(Experiment = paste0("XZ_2022", Experiment))

# drop dilutions from exp -131 that don't match the other experiments
aav_app <- aav_app %>% filter(!grepl("_(01|02|03|04|05|10|11|12)$", SAMPLE))

# make Sample names match across all experiments
aav_app$SAMPLE <- gsub("^S0", "S", aav_app$SAMPLE)

#samples that starts with S only (drop neg and pos controls)
aav_app <- aav_app[str_detect(aav_app$SAMPLE, "S"), ]
#exclude not S1-S7
aav_app <- aav_app %>% filter(!grepl("^S(8|9|10|11|12)", SAMPLE))

aav_app<-separate(aav_app,SAMPLE,into=c("Sample_num","Replicate"),sep = "_")

# update target names to match Exp_131
aav_app<- aav_app %>% 
  mutate(TARGET = str_replace_all(TARGET, c("1"="prom","2"="pA")))

#filtered arbitrarily by prom target (since droplet info is in duplicate)
aav_app <- aav_app %>% filter(TARGET=="prom") 

#create columns with theoretical expected % integrity + copy number
aav_app$expected<-rep(as.numeric(c("84","63","42","21","10","5","0")),each=4,6)
aav_app$copies<-rep(as.numeric(c("498","249","125","62")),42)

exp_aav = paste(aav_app$expected, "%", sep = "")
shiny_aav <- aav_GI_plot(aav_app,copies, ESTIMATED_PERCENTAGE) + theme(legend.position = "none") +
  geom_text(aes(650,expected,label = exp_aav, vjust = -0.5),show.legend = FALSE)

library(patchwork)
BR1_aav+shiny_aav + plot_annotation(tag_levels = "A")
###############################################################################
#linearity linkage vs shiny
###############################################################################

# avg %integrity and %recovery of biorad models
aav_avg_raw <- aav_prom %>% 
  filter(!expected== c("0")) %>% 
  group_by(Sample_num,
           Experiment,
           expected) %>% 
  summarise(avg_BR1= mean(BR1),
            stdevBR1=sd(BR1),
            rsdBR1=(stdevBR1/avg_BR1*100),#calculate % RSD of BR1 model
            rec_BR1=mean(BR1/expected*100))#calculate % recovery of BR1 model

aav_avg_app <- aav_app %>% 
  filter(!expected== c("0")) %>% 
  group_by(Sample_num,
           Experiment,
           expected) %>% 
  summarise(avg_shiny= mean(ESTIMATED_PERCENTAGE),
            stdev_shiny=sd(ESTIMATED_PERCENTAGE),
            rsd_shiny=(stdev_shiny/avg_shiny*100),#calculate % RSD of poisson-multinomial model
            rec_shiny=mean(ESTIMATED_PERCENTAGE/expected*100))#calculate % recovery of poisson-multinomial model

# Table 4
aav_lin<-merge(aav_avg_raw, aav_avg_app) 
aav_lin %>% kbl(digits = 1) %>% kable_styling(latex_options = "striped") #data for manuscript table 4
aav_lin<-pivot_longer(aav_lin,cols=c(avg_BR1,avg_shiny), names_to=c("model"))

#### Start of GLS (generalized least squared) model ####

# split data into two, and then apply GLS to each
dat1 <- aav_lin %>% filter(model == "avg_BR1")
dat2 <- aav_lin %>% filter(model == "avg_shiny")

library(nlme)
library(broom.mixed)

# correlation=corCompSymm(form = ~ 1 | Sample_num) indicating
# Sample_num is the cluster variable, and we are correcting for
# within-cluster correlation among replicates

# Apply gls to the averaged BR1 data
gfit1 <- nlme::gls(value ~ expected, data=dat1,
                   correlation=corCompSymm(form = ~ 1 | Sample_num))

# view model summary information
broom.mixed::tidy(gfit1)

# residuals of gfit1
est.gls.1 <- predict(gfit1, newdata = dat1)
dat1$est.gls <- est.gls.1
dat1$resid.gls <- with(dat1, value - est.gls)

# Apply gls to the averaged Shiny data
gfit2 <- nlme::gls(value ~ expected, data=dat2,
                   correlation=corCompSymm(form = ~ 1 | Experiment))
est.gls.2 <- predict(gfit2, newdata = dat2)
dat2$est.gls <- est.gls.2
dat2$resid.gls <- with(dat2, value-est.gls)

broom.mixed::tidy(gfit2)

# combine GLS fit of dat1 and dat2 to make plot
dat <- dat1 %>% bind_rows(dat2)

# because this plot look similar to Fig6, with same color scheme.
dat %>% ggplot(aes(x = est.gls, y = resid.gls, color = model)) + geom_point() + 
  geom_hline(yintercept = 0,linetype="dashed", color="darkgrey") +
  theme_classic(base_size = 12)+
  xlab("Fitted values")+
  ylab("Residuals")+ 
  #ggtitle("Fitted vs Residuals, Generalized Least Squared (GLS)")+   
  scale_color_manual(labels = c("Linkage (avg)", "Poisson-multinomial"),
                     values = c("#6388b4" ,"#c3bc3f", "#bb7693")) 


# calculate pseudo-R2 for GLS
cor(dat1$value,dat1$est.gls)^2
cor(dat2$value,dat2$est.gls)^2

# linearity plot (manuscript fig 6)
dat %>% 
  ggplot(aes(expected, value,color=model)) + 
  geom_point() +
  geom_line(aes(y = est.gls), linewidth = 1,)+
  theme_classic(base_size = 12)+
  annotate("text",
           x = 5, y = 90, 
           label = "italic(y) == 8.31 + 1.05 * italic(x) * ',' ~~ italic(R) [pseudo] ^2 ~ '=' ~ 0.989",
           parse = TRUE,color = "#6388b4", hjust = 0, vjust = 0) +
  annotate("text",
           x = 5, y = 90, 
           label = "italic(y) == -1.16 + 1.01 * italic(x) * ',' ~~ italic(R) [pseudo] ^2 ~ '=' ~ 0.996",
           parse = TRUE,color = "#c3bc3f", hjust = 0, vjust = 1.1) +
  xlab("Expected integrity (%)")+
  ylab("Calculated integrity (%)")+
  scale_color_manual(labels = c("Linkage (avg)", "Poisson-multinomial"), 
                     values = c("#6388b4" ,"#c3bc3f"))

