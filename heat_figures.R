library(tidyverse)
#############################################################################
# process raw data files
# compute AAV %genome integrity with different models
#############################################################################

# read in raw ddPCR files from heat_files folder
setwd("~/R/genome_integrity/heat_files")
heat_raw_file_list <- list.files(path = "~/R/genome_integrity/heat_files",pattern="^raw")
heat_compiled<-heat_raw_file_list %>% read_csv(id = "file")

#select relevant cols of data
heat_raw <- select(heat_compiled,Sample, Target,Concentration,`Ch1+Ch2+`,`Ch1-Ch2-`,AcceptedDroplets,Linkage)
heat_raw$experiment <- rep(c("02JUNE23","08JUN23"),each=100)

#samples that dont start with N (drop NC)
heat_raw <- heat_raw [!str_detect(heat_raw$Sample, "N"), ]

# order the data
heat_raw <- heat_raw %>% arrange(Sample)

heat_raw$treat <- rep(c("95C_0","95C_01","95C_05","95C_10","95C_20","95C_30"),each=24)

heat_raw$expected <- rep(as.numeric(c("48"),c(144)))

# filter into separate dfs by target (since droplet info is in duplicate)
heat_T1<- heat_raw %>% filter(Target=="1") 
heat_T2<- heat_raw %>% filter(Target=="2") 

# calculate %integrity with: BioRad %Linkage calculation equation 1 from Pranter 2023 S1 Text
heat_T1$BR1 <- heat_T1$Linkage/((heat_T1$Concentration+heat_T2$Concentration)/2)*100

# calculate %integrity with: BioRad %Linkage calculation equation 3 from Pranter 2023 S1 Text
#first solve BioRad eq 2 from S1 text
heat_T1$linkcomp<-heat_T1$Linkage+(abs(heat_T1$Concentration-heat_T2$Concentration))
# second find max cmv vs hGH_pA concentration and solve eq 3 
heat_max_vals<-heat_raw%>%
  group_by(Sample) %>%
  summarize(max_val = max(Concentration, na.rm=TRUE))

heat_T1 <- merge(heat_T1, heat_max_vals)
heat_T1$BR3 <- heat_T1$linkcomp/heat_T1$max_val*100

# theoretical copies for treateriment -131
heat_T1$copies<- rep(as.numeric(c("960","480","240","120","60","30")),each=2,6)

heat_T1 <- heat_T1 %>% 
  separate(Sample,into=c("Sample_num","Replicate"),sep = "_") 

# avg %integrity and %recovery of biorad models
heat_avg_compiled <- heat_T1 %>% 
  group_by(Sample_num,treat,expected) %>% 
  summarise(
    avg_BR1= mean(BR1),
    avg_BR3= mean(BR3),
    rec_BR1=mean(BR1/expected*100),
    rec_BR3=mean(BR3/expected*100))


heat_T1_longer<-heat_T1 %>% 
  pivot_longer(
    cols = c("BR1","BR3"), 
    names_to = "model", 
    values_to = "model_result"
  )


############################################################################################################
# process web app data (poisson multinomial results)
############################################################################################################
heat_app_file_list <- list.files(path = "~/R/genome_integrity/heat_files",pattern="^app")
heat_app_compiled<-heat_app_file_list %>% read_csv(id = "file")

heat_app <- select(heat_app_compiled,SAMPLE, TARGET,ESTIMATED_PERCENTAGE)
heat_app$experiment <- rep(c("02JUNE23","08JUN23"),each=100)

#samples that dont start with N (drop NC)
heat_app <- heat_app [!str_detect(heat_app$SAMPLE, "N"), ]

# order the data
heat_app <- heat_app %>% arrange(SAMPLE)

heat_app$treat<- rep(c("95C_0","95C_01","95C_05","95C_10","95C_20","95C_30"),each=24)

heat_app$expected <- rep(as.numeric(c("48"),c(144)))

heat_app_T1<- heat_app %>% filter(TARGET=="1") 
heat_app_T2<- heat_app %>% filter(TARGET=="2") 

heat_app_T1$copies<- rep(as.numeric(c("960","480","240","120","60","30")),each=2,6)

heat_app_T1 <- heat_app_T1 %>% 
  separate(SAMPLE,into=c("Sample_num","Replicate"),sep = "_") 

heat_app_compiled <- heat_app_T1 %>% 
  group_by(Sample_num,treat,expected) %>% 
  summarise(
    avg_integ= mean(ESTIMATED_PERCENTAGE),
    rec_integ=mean(ESTIMATED_PERCENTAGE/expected*100))

heat_app_T1$model <- rep(c("shiny"),c(72))

heat_app_T1 <- heat_app_T1 %>% 
  rename("model_result" = "ESTIMATED_PERCENTAGE")

compiled<-bind_rows(heat_app_T1, heat_T1_longer)
compiled <- compiled %>% arrange(model)
# New facet label names
treat_labels <- c("0 min.","1 min.","5 min.","10 min.","20 min.","30 min.")
names(treat_labels) <- c("95C_0","95C_01","95C_05","95C_10","95C_20","95C_30")


cols<-c("#8cc2ca" ,"#c3bc3f" ,"#55ad89")
r<-ggplot(compiled, aes(model, model_result,color=model)) + 
  stat_summary(geom = "bar", fun = mean, aes(fill = model)) +
  facet_grid(~ treat)+
  geom_hline(aes(yintercept=expected),linetype='dashed', color="grey35",linewidth=1.1) +
  scale_fill_manual(values=cols,name="Model",
                    labels=c("Linkage (avg)","Linkage (comp)",
                             "Poisson multinomial"))+
  guides(color = "none", size = "none")+
  theme_classic(base_size = 12)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +  
  ylab("Calculated integrity (%)")+
  scale_colour_manual(values=cols)
r+facet_grid(~ treat,labeller = labeller(treat = treat_labels))
