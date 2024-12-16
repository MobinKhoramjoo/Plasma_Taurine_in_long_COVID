{
library(readxl)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(reshape2)
library(gtsummary)
library(survival)
library(survminer)
library(pammtools)
library(adjustedCurves)
library(tidycmprsk)
library(ggsurvfit)
  
}

## Import supplemental data from https://doi.org/10.1016/j.xcrm.2023.101254
A<- read_excel("supplemental Data 2.xlsx", sheet = 1)
metadata <- read_excel("supplemental Data 2.xlsx", sheet = 2)

# data cleaning:
data <- data.frame(sapply(A, function(x) as.numeric(as.character(x)))) #NAs to empty cells
data[,1:4] <- A[,1:4] #import the annotation columns again
rm(A) 
colnames(data)[1] <- colnames(metadata)[1] #same ID for patients

sum(is.na(data) == TRUE) #number of NAs in the data
sum(data == 0, na.rm = TRUE) # number of zero values in data
data <- replace(data, data == 0, NA) #making the zero values to NAs
## transpose the data: 
t.data <- as.data.frame(t(data))

## missing value vector:
p <-c()
for (i in 1:nrow(t.data)) {
  p[i] <- sum(is.na(t.data[i,]))/262
} 
#input the vector to the data
t.data <- t.data %>% 
  mutate(percent_of_missing_Vlues= p) %>% 
  select(percent_of_missing_Vlues, everything())

#filter the data based on the generated column (percent of missing values)
filtered_t.data <- t.data %>% 
  filter(percent_of_missing_Vlues < 0.5)
##re-transpose the data:
cleaned_data <- as.data.frame(t(filtered_t.data))
cleaned_data <- cleaned_data[-1,]
row.names(cleaned_data) <- row.names(data)

## imputation of remaining missing values
for (i in 5:782){
  for(j in 1:262){
    if (is.na(cleaned_data[j,i]) == "TRUE"){
      cleaned_data[j,i] = min(cleaned_data[,i], na.rm = TRUE)
    }
  }
}

##check NAs and zero values to be removed
sum(is.na(cleaned_data))
sum(cleaned_data[,5:782]== 0, na.rm = TRUE)
write.csv(cleaned_data, "cleaned & imputted data (long, control).csv")

#normalize the data:
cleaned_data[,5:782] <- data.frame(sapply(cleaned_data[,5:782], function(x) as.numeric(x))) #make the columns numeric
normal_data <- cbind(cleaned_data[,1:4], as.data.frame(log10(cleaned_data[,5:782])+10))
normal_data  

long_ctr <- cleaned_data %>% 
  filter(Type != "Acute" ) ## excluding the acute samples

count_LC_C <- long_ctr[,5:782] %>% ### making the count matrix 
  t() %>% 
  as.data.frame()
colnames(count_LC_C) <- long_ctr$`Study ID`
count_LC_C <- as.matrix(count_LC_C)

coldata <- long_ctr[,2] %>% ### making the meta data
  as.data.frame()
row.names(coldata)<- long_ctr$`Study ID`
colnames(coldata) <- 'Type'
coldata$Type <- factor(coldata$Type)


sapply(as.data.frame(count_LC_C), function(x) which(is.na(x))) #locate the NA value in matrix
count_LC_C[343,11] <- mean(count_LC_C[343,],na.rm = TRUE) #replace it with mean value

####  Differentially expression analysis: 
y <- DGEList(counts=count_LC_C, group=coldata$Type)
y <-estimateCommonDisp(y)
y<-calcNormFactors(y)

FC <- exactTest(y)
topTags(FC)
FC$table ## table of DE analysis 


########################### making the volcano plot from MetaboAnalyst data: 
FC_LC_CTRL <- read_excel("volcano.xlsx")
colnames(FC_LC_CTRL)[1] <- "mols"
FC_LC_CTRL$`-LOG10(p)'`

LCvsCTRL <- FC_LC_CTRL %>% 
  select("mols", "log2(FC)", "p.ajusted")

range(LCvsCTRL$`log2(FC)`)
log2(2)

cat_LC_CTRL <- c() #create a vector for categories
for (i in 1:length(row.names(LCvsCTRL))) { #assign corresponding names to vector
  if(LCvsCTRL$`log2(FC)`[i] > log2(2) & LCvsCTRL$p.ajusted[i] < 0.05){
    cat_LC_CTRL[i] <- "Upregulated"
  }
  else if(LCvsCTRL$`log2(FC)`[i] < -log2(2) & LCvsCTRL$p.ajusted[i] < 0.05){
    cat_LC_CTRL[i] <- "Downregulated"
  }
  else{
    cat_LC_CTRL[i] <- "Non-significant"
  }
}

LCvsCTRL <- LCvsCTRL %>% # assigning the vector to a column 
  mutate(category = cat_LC_CTRL)

showing_label_LC_CTRL <- c() #create a vector for the molecules that should be shown 
for (i in 1:length(row.names(LCvsCTRL))) { # assigning corresponding lables
  if(LCvsCTRL$category[i] == "Upregulated" | LCvsCTRL$category[i] == "Downregulated"){
    showing_label_LC_CTRL[i] <- LCvsCTRL$mols[i]
  }
  else {
    showing_label_LC_CTRL[i] <- NA
  }
}

LCvsCTRL <- LCvsCTRL %>% #assigning the created vector to a column
  mutate(present_label = showing_label_LC_CTRL)

LCvsCTRL<- LCvsCTRL %>% 
  mutate(FDR= -log10(p.ajusted))
sum(LCvsCTRL$category == "Downregulated")
2^LCvsCTRL$`log2(FC)`[3]


ggplot(data = LCvsCTRL, 
       aes(x= `log2(FC)`,
           y= `FDR`, 
           colour= `category`)) + 
  geom_text_repel(aes(label= `present_label`),size = 4.5, force = 3, max.overlaps = 1.8) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-2, 8)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title="Convalescence vs control")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="bottom", 
        legend.title = element_blank())+
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 14, face= "bold"), 
        legend.position = "top")
####### boxplot for control acute long:
Data_a.c.t <- data %>% 
  select(Type, Taurine)
Data_a.c.t$Type <- as.character(Data_a.c.t$Type)
mycomparison <- list(c("Control","Acute"), c("Acute","Convalescence"), c("Control","Convalescence"))


ggplot(Data_a.c.t, aes(x = factor(Type, levels = c("Control", "Acute", "Convalescence")), y = Taurine, color = Type)) + 
  geom_boxplot(outlier.shape = NA, color= "black", size= 0.75) + 
  ggtitle("") +
  ylab("Plasma Taurine (μM)") +
  xlab("") +
  ylim(0,810)+
  scale_color_manual(values = c("#EC6B56","#377B2B", "#FFC154")) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14), axis.text = element_text(size = 14, face= "bold"), legend.position = "none") +
  geom_jitter(size = 3, aes(col = Type), alpha= 0.5) +
  stat_compare_means(comparisons = mycomparison,
                     label = "p.format",
                     label.y = NULL,
                     na.rm = TRUE,
                     method = "wilcox.test",
                     bracket.size = 1,
                     size= 7)


############## dot & line plot:

par(mfcol = c(1,2))  

line_increased <- read_excel('line.xlsx', sheet = 'Increased')

increase <- ggpaired(line_increased,
                     id= 'Study ID',
                     cond1 = "Acute", 
                     cond2 = "Convalescence",
                     fill = "condition", 
                     palette = c("#EC6B56", "#FFC154"), 
                     line.color = "darkgray", 
                     point.size = 1.5,
                     width= 0.3,)+
  ggtitle("") +
  ylab("Plasma Taurine (μM)") +
  xlab("") +
  ylim(0,700)+
  theme_bw() +
  theme(axis.title.y = element_text(size = 14), 
        axis.text = element_text(size = 14, face= "bold"), 
        legend.position = "none")+
  stat_compare_means(label = "p.format",
                     label.x = 1.45,
                     label.y = 650,
                     na.rm = TRUE,
                     method = "wilcox.test",
                     bracket.size = 1,
                     size= 7,
                     paired = TRUE)

line_decreased <- read_excel('line.xlsx', sheet = 'Decreased')

decrease <- ggpaired(line_decreased,
                     id= 'Study ID',
                     cond1 = "Acute", 
                     cond2 = "Convalescence",
                     fill = "condition", 
                     palette = c("#EC6B56", "#FFC154"), 
                     line.color = "darkgray", 
                     point.size = 1.5,
                     width= 0.3,)+
  ggtitle("") +
  ylab("Plasma Taurine (μM)") +
  xlab("") +
  ylim(0,700)+
  theme_bw() +
  theme(axis.title.y = element_text(size = 14), 
        axis.text = element_text(size = 14, face= "bold"), 
        legend.position = "none")+
  stat_compare_means(label = "p.format",
                     label.x = 1.45,
                     label.y = 650,
                     na.rm = TRUE,
                     method = "wilcox.test",
                     bracket.size = 1,
                     size= 7,
                     paired = TRUE)

combined_plot <- plot_grid(increase, decrease, ncol = 2, align = 'v')
combined_plot
#p value calculation

melted_decrease<- melt(line_decreased) #decreased 
wilcox.test(value~variable,data = melted_decrease, paired = TRUE)
median(line_decreased$Acute)
median(line_decreased$Convalescence)

melted_increase<- melt(line_increased) #increased 
wilcox.test(value~variable,data = melted_increase, paired = TRUE)
median(line_increased$Acute)
median(line_increased$Convalescence)

######################## Taurine and symptoms
## Import supplemental data from https://doi.org/10.1016/j.xcrm.2023.101254
symptoms <- read_excel("Supplemental Data 2.xlsx", sheet = 4)

ggplot(data = symptoms, 
       aes(x= factor(Symptom,levels = unique(symptoms$Symptom)), y=Taurine)) + 
  geom_boxplot(outlier.shape = NA, 
               aes(fill=Status),
               width = 0.55) +
  ylim(0,400)+
  ylab("Plasma Taurine (μM)") +
  xlab("") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "right")  

sym <- read_excel("Supplemental Data 2.xlsx", sheet = 3)


wilcox.test(Taurine~Fatigue, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~Weakness, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~`Chills/Temperature lability`, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~`Night Sweats`, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~`Runny nose`, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~`Muscle ache`, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~SoB, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~`Chronic cough`, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~Palpitation, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~Tachycardia, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~`Cognitive impairment`, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~Insomnia, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~`Change in smell`, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~`Change in taste`, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~Headache, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~`Mood disturbance`, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~Nausea, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~`Abdominal Pain`, data = sym, paired = FALSE)$p.value
wilcox.test(Taurine~Diarrhea, data = sym, paired = FALSE)$p.value


######### dot and line plot for with and event free groups: 

with_event <- read_excel('line.xlsx', sheet = 4)

with_event_1 <- ggpaired(with_event,
                         id= 'Study ID',
                         cond1 = "Acute", 
                         cond2 = "Convalescence",
                         fill = "condition", 
                         palette = c("#EC6B56", "#FFC154"), 
                         line.color = "darkgray", 
                         point.size = 1.5,
                         width= 0.5,)+
  ggtitle("With-event") +
  ylab("Plasma Taurine (μM)") +
  xlab("") +
  ylim(0,420)+
  theme_bw() +
  theme(axis.title.y = element_text(size = 14), 
        axis.text = element_text(size = 14, face= "bold"), 
        legend.position = "none")+
  stat_compare_means(label = "p.format",
                     label.x = 1.45,
                     label.y = 410,
                     na.rm = TRUE,
                     method = "wilcox.test",
                     bracket.size = 1,
                     size= 7,
                     paired = TRUE)
with_event_1

no_event <- read_excel('line.xlsx', sheet = 5)

no_event_1 <- ggpaired(no_event,
                       id= 'Study ID',
                       cond1 = "Acute", 
                       cond2 = "Convalescence",
                       fill = "condition", 
                       palette = c("#EC6B56", "#FFC154"), 
                       line.color = "darkgray", 
                       point.size = 1.5,
                       width= 0.5,)+
  ggtitle("Event-free") +
  ylab("Plasma Taurine (μM)") +
  xlab("") +
  ylim(0,420)+
  theme_bw() +
  theme(axis.title.y = element_text(size = 14), 
        axis.text = element_text(size = 14, face= "bold"), 
        legend.position = "none")+
  stat_compare_means(label = "p.format",
                     label.x = 1.45,
                     label.y = 410,
                     na.rm = TRUE,
                     method = "wilcox.test",
                     bracket.size = 1,
                     size= 7,
                     paired = TRUE)

no_event_1


(median(no_event$Convalescence)/median(no_event$Acute))
median(with_event$Convalescence)/median(with_event$Acute)


combined_plot_event <- plot_grid(no_event_1, with_event_1, ncol = 2, align = 'v')
combined_plot_event

############# QoL and the symptoms of decreased and increased groups:
a <- read_excel('line.xlsx', sheet = 1)

a <- a %>% 
  filter(Phase == 'Convalescence') %>% 
  select(`Study ID`, Taurine, Status) %>% 
  mutate(`EQ-VAS` = as.numeric(metadata$`EQ-VAS`), `Total SF-12`= as.numeric(metadata$`Total SF-12`))

ggplot(a, aes(x= Status, y= `Total SF-12`))+ geom_boxplot()+
  stat_compare_means(label = "p.signif",
                     label.x = 1.45,
                     label.y = 110,
                     na.rm = TRUE,
                     method = "wilcox.test",
                     bracket.size = 1,
                     size= 7,
                     paired = FALSE) ### both QoL scores are not significant

colnames(b)


b <- read_excel('Supplemental Data 2.xlsx', sheet = 3)
b <- b %>% 
  mutate(Status = a$Status)

b %>% 
  select(Status, Fatigue, Weakness, `Chills/Temperature lability`, `Night Sweats`, `Runny nose`, 
         `Muscle ache`, SoB, `Chronic cough`, Palpitation, Tachycardia, `Cognitive impairment`, 
         Insomnia, `Change in smell`, `Change in taste`, `Headache`, `Mood disturbance`, Nausea,
         `Abdominal Pain`,Diarrhea, ) %>%
  tbl_summary(
    by= Status) %>% 
  add_n() %>% 
  add_overall() %>% 
  add_p()


ggplot(b, aes(x= Status, y= `N of symptoms`))+ geom_boxplot()+
  stat_compare_means(label = "p.signif",
                     label.x = 1.45,
                     label.y = 11,
                     na.rm = TRUE,
                     method = "wilcox.test",
                     bracket.size = 1,
                     size= 7,
                     paired = FALSE) ### there is no significant difference between the # of symptoms in decreased or increased taurine.


############# stack bar plot:
stack <- read_excel('line.xlsx', sheet= 1)
stack <- stack[118:234, ]

contable <-as.data.frame(table(stack$Status[1:117], stack$Event[1:117]))
contable <- contable %>% 
  mutate(percent = c(45,74.2, 55, 25.8))

ggplot(contable, aes(x = Var1, y = percent, fill = Var2)) + 
  geom_bar(position = "stack", stat = "identity") +
  ylab("Percent") +
  scale_fill_manual(values = c("#0000D1", "#D60000")) +  # Change fill colors to blue and red
  theme_classic() +
  theme(
    legend.position = "top",     # Place the legend on top
    plot.margin = margin(1, 10, 2, 3, "cm")
  )

################ boxplot for event in convalescence:

mycomparison_event<- list(c("Event-Free","With Event"))
ggplot(stack, aes(x = factor(Event, levels = c("Event-Free", "With Event")), y = Taurine, color = Event)) + 
  geom_boxplot(outlier.shape = NA, color= "black", size= 0.75) + 
  ggtitle("") +
  ylab("Plasma Taurine (μM)") +
  xlab("") +
  ylim(0,600)+
  scale_color_manual(values = c("#0000D1", "#D60000")) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14), axis.text = element_text(size = 14, face= "bold"), legend.position = "none") +
  geom_jitter(size = 3, aes(col = Event), alpha= 0.5) +
  stat_compare_means(comparisons = mycomparison_event,
                     label = "p.format",
                     label.y = 500,
                     na.rm = TRUE,
                     method = "wilcox.test",
                     bracket.size = 1,
                     size= 7)
?stat_compare_means
################## boxplots for severity groups:
stack<- stack %>% 
  mutate(severity = long_ctr$severity[29:145])
mycomparison_severity<- list( c("recovered","severe"))
ggplot(stack, aes(x = factor(severity, levels = c("recovered", "mild", "severe")), y = Taurine, color = severity)) + 
  geom_boxplot(outlier.shape = NA, color= "black", size= 0.75) + 
  ggtitle("") +
  ylab("Plasma Taurine (μM)") +
  xlab("") +
  ylim(0,670)+
  scale_color_manual(values = c("#b873d1","#69B0F8","#F2757B")) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14), axis.text = element_text(size = 14, face= "bold"), legend.position = "none") +
  geom_jitter(size = 3, aes(col = severity), alpha= 0.8) +
  stat_compare_means(comparisons = mycomparison_severity,
                     label = "p.signif",
                     label.y = NULL,
                     na.rm = TRUE,
                     method = "wilcox.test",
                     bracket.size = 1,
                     size= 7,
                     hide.ns= TRUE)
########### survival analysis:
tau_data <- read_excel('tau_outcome.xlsx')

tau_data <- tau_data %>% 
  mutate(Event = as.factor(Event))
##### performing cumulative hazard plot (univariate)
A<-survfit(Surv(tau_data$difftime.m, 
                as.numeric(as.character(tau_data$Event)))~tau_data$tau, 
           data = tau_data)
summary(A)
plot(A)

tau.KMcurve<- ggsurvplot(A, fun = "cumhaz", 
                         pval = TRUE,
                         legend.labs=c("Decreased Taurine", "Increased Taurine"),
                         palette = c("red","blue"),
                         size= 1.5,
                         cumevents = FALSE,
                         break.time.by = 6,
                         risk.table = 'nrisk_cumevents')

tau.KMcurve$plot<-tau.KMcurve$plot + expand_limits(y=1)
tau.KMcurve

########### performing crr and cox regression for Taurine (multivariate)

adjust_data <- tau_data[,c(-1,-3,-4)]
adjust_data <- adjust_data %>% 
  mutate(tau = ifelse(tau == 'increased',1,0))
colnames(adjust_data)[3] <- "Taurine"

adjust_data <- data.frame(adjust_data)
colnames(adjust_data)[7] <- 'Supplemental O2'
colnames(adjust_data)[12] <- 'WHO Scale'
class(adjust_data$`WHO Scale`) <- "numeric"
adjust_data$Taurine<-as.factor(adjust_data$Taurine)

cox <- coxph(Surv(difftime.m, ifelse(Event == 1,1,0)) ~ Taurine
             + Age
             + Sex
             + Diabetes
             + CKD
             + `Supplemental O2`
             + Dexamethasone
             + Antibiotic
             + Tocilizumab
             + Remdesvir
             + Vaccination
             + `WHO Scale`,
             data = adjust_data)
summary(cox)
cox %>% 
  tbl_regression(exp = TRUE) %>% 
  bold_p(t = 0.05)

ggforest(cox, data = adjust_data,
         fontsize = 1,
         cpositions = c(0.01, NA, 0.25),
         refLabel = NA)

##### scatter plot for taurine and other biomarkers of ineterest. 

ggplot(normal_data, aes(x= Taurine, y= Lipopolysaccharide.binding.protein)) + 
  geom_point(size = 3, alpha = 0.75)+
  geom_smooth(method=lm)+
  stat_cor(method = "pearson", size = 8, label.x = 2.2) +
  theme_bw()+
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold")) +
  xlab("log10(Taurine)") + ylab("log10(LBP)")
