library(tidyverse)
library(kableExtra)
#############################################################################
# plotting function & color blindness friendly color palette
cols<-c("#767676", "#6388b4", "#55ad89", "#ffae34", "#8cc2ca", "#ef6f6a", "#c3bc3f", "#bb7693" )
#plotting function
GI_plot <- function(dataset, x, y) {
  ggplot(dataset, aes({{x}}, {{y}},color=Sample_num)) + 
    geom_point() +  
    geom_line(aes(group=Sample_num),stat = "summary", fun = mean, linewidth=1) +
    geom_hline(aes(yintercept=expected,color=Sample_num),linetype='dashed',linewidth=1) +
    scale_colour_manual(values=cols,breaks=c("S8","S7","S6","S5","S4","S3","S2","S1"),
                        name="% Integrity sample",
                        labels=c("100%","82%","60%","43%","29%","18%","8%","0%"))+
    scale_x_log10()+
    theme_classic(#base_size = 14
    )+
    ylim(-3.5,101)+
    xlab(expression(paste("Concentration (copies/",mu,"L)")))+
    ylab("Calculated integrity (%)")
}

#############################################################################
# # clean and sort data, compute genome integrity with different models
#############################################################################

# read in raw ddPCR files from plasmid_files folder
setwd("~/R/genome_integrity/plasmid_files")
raw_file_list <- list.files(path = "~/R/genome_integrity/plasmid_files",pattern="^raw")
raw_compiled<-raw_file_list %>% read_csv(id = "file")

# select relevant cols, add experiment ID,sep sample & replicate
raw_compiled <- select(raw_compiled,file, Sample, Target,Concentration,`Ch1+Ch2+`,`Ch1-Ch2-`,AcceptedDroplets,Linkage)

# keep beginning of file name as Experiment ID
raw_compiled <- raw_compiled %>% 
  separate_wider_delim(file, "data_", names = c(NA, "Experiment")) %>% 
  separate_wider_delim(Experiment, ".csv", names = c("Experiment",NA))

# update target names to match manuscript
raw_compiled <- raw_compiled %>%
  mutate(Target = str_replace_all(Target, c("1" = "pA", "2" = "CMV")))

# separate Sample into sample # and replicate #
raw_compiled<-separate(raw_compiled,Sample,into=c("Sample_num","Replicate"),sep = "_")

# filter into separate dfs by target (since droplet info is in duplicate)
CMV<- raw_compiled %>% filter(Target=="CMV") 
pA<- raw_compiled %>% filter(Target=="pA") 

# calculate %integrity with formula 1: double positive/total positive droplets
CMV$percentmodel <- CMV$`Ch1+Ch2+`/(CMV$AcceptedDroplets-CMV$`Ch1-Ch2-`)*100

# calculate %integrity with: BioRad %Linkage calculation equation 1 from Pranter 2023 S1 Text
CMV$BR1 <- CMV$Linkage/((CMV$Concentration+pA$Concentration)/2)*100

# calculate %integrity with: BioRad %Linkage calculation equation 3 from Pranter 2023 S1 Text
# first solve BioRad eq 2 from S1 text
CMV$linkcomp<-CMV$Linkage+(abs(CMV$Concentration-pA$Concentration))
# second find max CMV vs pA concentration and solve eq 3 
max_vals<-raw_compiled%>%
  group_by(Sample_num,Replicate,Experiment) %>%
  summarize(max_val = max(Concentration, na.rm=TRUE))

CMV <- left_join(CMV, max_vals, by = join_by(Experiment == Experiment,Replicate == Replicate,Sample_num == Sample_num, ))

CMV$BR3 <- CMV$linkcomp/CMV$max_val*100

# create columns with theoretical expected % integrity + copy number
CMV$expected <-rep(as.numeric(c("0","8","18","29","43","60","82","100")),each=12,2)

CMV$copies <- rep(as.numeric(c("5000","4000","3000","2000","1000","500","250","125",
                               "63","31","16","8")),16) 

# Table 1: Percent recovery of integrity calculated by percentage of double positive droplets/total positive droplets
table1 <- CMV %>% 
  group_by(copies, expected) %>% 
  summarise(rec_percentmodel=mean(percentmodel/expected*100))
table1 <-table1 %>% pivot_wider(names_from = expected, values_from = rec_percentmodel)
table1$`0`[table1$`0`  %in% c("Inf", "NaN")] <- "NA"
table1 %>% kbl(digits = 1) %>% kable_styling(latex_options = "striped")

# Table 2:Percent recovery of integrity calculated by linkage models
table2 <- CMV %>% 
  group_by(Sample_num) %>% 
  summarise(avg_BR1= mean(BR1),
            avg_BR3= mean(BR3),
            rsd_BR1=sd(BR1)/mean(BR1)*100,
            rsd_BR3=sd(BR3)/mean(BR3)*100,
            rec_BR1=mean(BR1/expected*100),
            rec_BR3=mean(BR3/expected*100))
table2 %>% kbl(digits = 1) %>% kable_styling(latex_options = "striped")

# plot model results

# labels for geom_hline lines
exp_val = paste(CMV$expected, "%", sep = "")

#plot results for each model
percentmodel<-GI_plot(CMV,copies, percentmodel)  + theme(legend.position = "none") +expand_limits(x = 8500)+  
  geom_text(aes(8000,expected,label = exp_val, vjust = -0.5),show.legend = FALSE)
percentmodel

BR1<-GI_plot(CMV,copies, BR1) + theme(legend.position = "none") +expand_limits(x = 8500)+ #ggtitle("BioRad %Linked, eq.1") +   
  geom_text(aes(8000,expected,label = exp_val, vjust = -0.5),show.legend = FALSE)

BR3<-GI_plot(CMV,copies, BR3) + theme(legend.position = "none") +expand_limits(x = 8500)+#ggtitle("BioRad %Linked, eq.3") +   
  geom_text(aes(8000,expected,label = exp_val, vjust = -0.5),show.legend = FALSE)


############################################################################################################
# clean and sort analysed data from Poisson-multinomial Shiny app, plot calculated integrity
############################################################################################################

app_file_list <- list.files(path = "~/R/genome_integrity/plasmid_files",pattern="^2021")
app_compiled<-app_file_list %>% read_csv(id = "file")

# select relevant cols, add experiment ID,sep sample & replicate
app_compiled <- select(app_compiled,file,CONCENTRATION, SAMPLE, TARGET,ESTIMATED_PERCENTAGE)

app_compiled <- app_compiled %>% 
  separate_wider_delim(file, "_app", names = c("Experiment",NA)) 

# update target names to match manuscript
app_compiled <- app_compiled %>%
  mutate(TARGET = str_replace_all(TARGET, c("1" = "pA", "2" = "CMV")))

#filtered arbitrarily by CMV target (since droplet info is in duplicate)
app_compiled <- app_compiled %>% filter(TARGET=="CMV") 

#create columns with theoretical expected % integrity + copy number
app_compiled$expected <-rep(as.numeric(c("0","8","18","29","43","60","82","100")),each=12,2)
app_compiled$copies <- rep(as.numeric(c("5000","4000","3000","2000","1000","500","250","125",
                                        "63","31","16","8")),16) 
app_compiled$Experiment <- rep(c("Exp_0805","Exp_0903"),each=96) 

app_compiled<-separate(app_compiled,SAMPLE,into=c("Sample_num","Replicate"),sep = "_")

# avg %integrity & percent recovery of both exps
table3 <- app_compiled %>% 
  group_by(Sample_num) %>% 
  summarise(avg_integ= mean(ESTIMATED_PERCENTAGE),
            stdev_integ=sd(ESTIMATED_PERCENTAGE),
            rsd_integ=sd(ESTIMATED_PERCENTAGE)/mean(ESTIMATED_PERCENTAGE)*100,
            recovery=mean(ESTIMATED_PERCENTAGE/expected*100))

table3 %>% kbl(digits = 1) %>% kable_styling(latex_options = "striped")

#plot shiny app results

#labels for geom_hline lines
expapp = paste(app_compiled$expected, "%", sep = "")

# results for shiny app
shiny<-GI_plot(app_compiled,copies, ESTIMATED_PERCENTAGE) + theme(legend.position = "none")+ 
  #ggtitle("Poisson multinomal %Full Genome")+   
  geom_text(aes(8000,expected,label = expapp, vjust = -0.5),show.legend = FALSE)+ expand_limits(x = 8500)

###############################################################################
#linearity linkage vs shiny
###############################################################################

#compare linearity of models
# summarize linkage models
linkage_linearity <- CMV %>% 
  group_by(Sample_num,
           Experiment,
           expected) %>% 
  summarise(avg_BR1= mean(BR1),
            avg_BR3= mean(BR3))%>%
  pivot_longer(cols=c(avg_BR1,avg_BR3), names_to=c("model"))

# summarize shiny app model
app_linearity <- app_compiled %>% 
  group_by(Sample_num,
           Experiment,
           expected) %>% 
  summarise(avg_shiny= mean(ESTIMATED_PERCENTAGE))%>%
  pivot_longer(cols=avg_shiny,names_to = "model")

linkagevsapp <- rbind(linkage_linearity,app_linearity)
linkagevsapp <-filter(linkagevsapp, !expected=="0")

#### Start of GLS (generalized least squared) model ####
library(broom.mixed)
library(nlme)
# split data into 3 models, and then apply GLS to each
dat1 <- linkagevsapp %>% filter(model == "avg_BR1")
dat2 <- linkagevsapp %>% filter(model == "avg_BR3")
dat3 <- linkagevsapp %>% filter(model == "avg_shiny")

# correlation=corCompSymm(form = ~ 1 | Sample_num) indicating
# Sample_num is the cluster variable, and we are correcting for
# within-cluster correlation among replicates

gfit1 <- nlme::gls(value ~ expected, data=dat1,
                   correlation=corCompSymm(form = ~ 1 | Sample_num))

gfit2 <- nlme::gls(value ~ expected, data=dat2,
                   correlation=corCompSymm(form = ~ 1 | Sample_num))

gfit3 <- nlme::gls(value ~ expected, data=dat3,
                   correlation=corCompSymm(form = ~ 1 | Sample_num))

# view model summary information
broom.mixed::tidy(gfit1)
broom.mixed::tidy(gfit2)
broom.mixed::tidy(gfit3)

# The "predict" function is actually nlme::predict.gls
est.gls.1 <- predict(gfit1, newdata = dat1)
dat1$est.gls <- est.gls.1
dat1$resid.gls <- with(dat1, value - est.gls)

est.gls.2 <- predict(gfit2, newdata = dat2)
dat2$est.gls <- est.gls.2
dat2$resid.gls <- with(dat2, value - est.gls)

est.gls.3 <- predict(gfit3, newdata = dat3)
dat3$est.gls <- est.gls.3
dat3$resid.gls <- with(dat3, value - est.gls)

# calculate pseudo-R2 for GLS
cor(dat1$value,dat1$est.gls)^2
cor(dat2$value,dat2$est.gls)^2
cor(dat3$value,dat3$est.gls)^2


# combine GLS fit of dat1 and dat2 to make plot
dat <- dat1 %>% bind_rows(dat2) %>% bind_rows(dat3)

# because this plot look similar to p2, with same color scheme.
dat %>% ggplot(aes(x = est.gls, y = resid.gls, color = model)) + geom_point() + 
  geom_hline(yintercept = 0,linetype="dashed", color="darkgrey") +
  theme_classic(base_size = 12)+
  xlab("Fitted values")+
  ylab("Residuals")+  
  scale_color_manual(labels = c("Linkage (avg)","Linkage (comp)", "Poisson-multinomial"), 
                     values = c("#6388b4" ,"#c3bc3f","#bb7693"))



# linearity plot
lin<-dat %>% 
  ggplot(aes(expected, value,color=model)) + 
  geom_point() +
  geom_line(aes(y = est.gls), linewidth = 1,)+
  theme_classic(base_size = 12)+
  annotate("text",
           x = 5, y = 95,
           label = "italic(y) == 14.3 +0.899 * italic(x) * ',' ~~ italic(R) [pseudo] ^2 ~ '=' ~ '0.972'",
           parse = TRUE,color = "#6388b4", hjust = 0, vjust = 0) +
  annotate("text",
           x = 5, y = 95,
           label = "italic(y) == 17.3 + 0.877 * italic(x) * ',' ~~ italic(R) [pseudo] ^2 ~ '=' ~ '0.970'",
           parse = TRUE,color = "#c3bc3f", hjust = 0, vjust = 1.1) +
  annotate("text",
           x = 5, y = 95,
           label = "italic(y) == -0.156 + 0.968 * italic(x) * ',' ~~ italic(R) [pseudo] ^2 ~ '=' ~ '1.0'",
           parse = TRUE,color = "#bb7693", hjust = 0, vjust = 2.1) +
  xlab("Expected integrity (%)")+
  ylab("Calculated integrity (%)")+
  scale_color_manual(labels = c("Linkage (avg)","Linkage (comp)", "Poisson-multinomial"), 
                     values = c("#6388b4" ,"#c3bc3f","#bb7693"))+ theme(legend.position = c(0.8, 0.3))


#plot linkage and shiny model results together 
#"&" adds the element to all subplots
library(patchwork)
(BR1+BR3)/(shiny+lin) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold'))
