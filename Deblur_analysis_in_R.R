# ABALONE WS Deblur Analysis in R

## Figure 1: ----

# Part 1: plot qPCR normalized copies over time in fecal samples ----

library(readr)
library(ggplot2)
abalone_dna_metadata_filtered <- read_delim("~/d_abalone_caxc/qiime2/abalone_dna_metadata_filtered.tsv", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)
fecal<-subset(abalone_dna_metadata_filtered, sample_type == "fecal")
fecal$collection_date<-as.Date(fecal$collection_date, "%m/%d/%Y")

ggplot(data = fecal, aes(collection_date, copies_normalized_pseudo, color = Treatment))+
  geom_point()+
  geom_smooth()+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  xlab("Date")+
  scale_y_log10()+
  scale_color_manual(values = c("darkcyan", "darkorange2"), labels = c("Control", "Exposed"))+
  theme_bw()

ggplot(data = fecal, aes(Time, copies_normalized_pseudo, color = Treatment))+
  geom_point()+
  geom_smooth()+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  xlab("Date")+
  scale_y_log10()+
  scale_color_manual(values = c("darkcyan", "darkorange2"), labels = c("Control", "Exposed"))+
  theme_bw()+
  scale_x_continuous(breaks=seq(0,11,1))+
  labs(x = "Month")

fecal<-subset(fecal, Exposed == "Positive")
fecal<-subset(fecal, mort_sacrifice == "sacrifice")
low<-subset(fecal, tube_id %in% c('31','32','43','48'))
ggplot(data = low, aes(Time, copies_normalized_pseudo, color = tube_id))+
  geom_point()+
  geom_line()+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  xlab("Date")+
  scale_y_log10()+
  theme_bw()+
  scale_x_continuous(breaks=seq(0,11,1))+
  labs(x = "Month")

# Part 2: microbiome volatility density plots ----

# Create vector of distances between timepoints for either condition
## timepoints include months: May - July, July - Oct, Oct - Dec, Dec - Jan, Jan - March, March - April
# animals 1-24
control<- c(0.403692901180735, 0.395824022730404, 0.522156617549389, 0.403144944926932, 0.557958772016018, 0.561206372394416,
           0.342060474493932, 0.291271208224926, 0.573780547538106, 0.642619250943927, 0.536111677289212, 0.558790629687434,
           0.30755922727333, 0.299319886805304, 0.369840621985648, 0.62043018213459, 0.434567179349663, 0.465501957825808,
           0.415484911086796, 0.313853782164732, 0.295074813158825, 0.536915182081554, 0.385870194466998, 0.469861055054011,
           0.436472721103403, 0.437248095360562, 0.341263957142199, 0.453654871438955, 0.518744449938829,
           0.460377245867727, 0.407562717460364, 0.552819069941584, 0.507012615253728, 0.323003870762575, 0.499718343987108,
           0.219288861690679, 0.515336593760697, 0.612251076281258, 0.611604002478737, 0.66166185588502,
           0.481873070526528, 0.521954112245838, 0.551776803376221, 0.391985571771912, 0.41120343381212, 0.47746719983543,
           0.335138034814036, 0.354966574282246, 0.58775585493793, 0.385751809476335, 0.402810072378008,
           0.508435228233891, 0.391765379541209, 0.664316203263437, 0.466037159117092, 0.397346811989761,
           0.374943751873819, 0.281525381947777, 0.331562899208252, 0.476662193005074, 0.505111628282757, 0.543543924080146,
           0.307963363011015, 0.195031958979167, 0.31229354982845, 0.520097308131422, 0.586219875680424, 0.514346204562978,
           0.519130292441425, 0.212037068839671, 0.399038736619873, 0.529724951039812, 0.324088732797038, 0.588671934726917,
           0.26728723904588, 0.321798323452508, 0.315208121327393, 0.767316367100379, 0.58397905942958, 0.376043317556674,
           0.464362931938917, 0.251860483452508, 0.520551511118791, 0.328889126513238, 0.500241521579254,
           0.627867910857291, 0.224279299682666, 0.356128244672314, 0.607777550678169, 0.406914899997422, 0.560711276210481,
           0.533367972111831, 0.520730215336581, 0.490731012758057, 0.707548127921525, 0.605206797167886, 0.257364807910352)
  
#animals 25-48
exposed<-c(0.336297611826307, 0.297417029398962, 0.287905532038042, 0.55641842183168, 0.434599706491461, 0.533519777475217,
           0.446608300385498, 0.50758529526185, 0.477110573096593, 0.312051534693181, 0.607359783930601,
           0.245709747443765, 0.305992821875149, 0.566525696566597, 0.314033463851792, 0.489241462233394,
           0.276853205792788, 0.357177229127925, 0.447842477424673, 0.50964976534647, 0.28517988301359, 0.422542354008808,
           0.507523269792513, 0.495008100773835, 0.32251083251387, 0.547670861461372, 0.319795374506177, 0.333059920046721,
           0.406858226783548, 0.23452884469437, 0.617359322273927, 0.406029857071407, 0.274454719682811,
           0.327987780241737, 0.599062973095681, 0.52891657754545, 0.468482976276331, 0.459243270274311,
           0.348657205601067, 0.408711445143216, 0.277584197536296, 0.437896470801139,
           0.247782390664317, 0.322683580696375, 0.537250116385777, 0.570978498312689, 0.532341699014778,
           0.320129163575081, 0.460330549507781, 0.353659854815383, 0.465996323680743, 0.404304526321923, 0.478660267878825,
           0.274571524173659, 0.244710461473932, 0.254643203645604, 0.440885892849073, 0.599496598762292, 0.575158141343031,
           0.353855250183178, 0.247102049304644, 0.321059012080125, 0.48174041890816, 0.478837851176654, 0.398234091707475,
           0.441650465204872, 0.461200478879983, 0.219835169071612, 0.709056604939942, 0.642645225028707, 0.62915737076822,
           0.429134648897105, 0.502777085240355, 0.583855961456432, 0.461167610095901, 0.61487218139559, 0.577941083730973,
           0.322590873017617, 0.306833992898455, 0.366242458458616, 0.503278680884526, 0.651137279845216, 0.65449292529986,
           0.36776174192756, 0.460634989182071, 0.617539958953279, 0.379431080583626, 0.571086274899624, 0.593356667797448,
           0.320689965999387, 0.216547007319426, 0.46153506153269, 0.608119794380021, 0.556408784205222)

df<- data.frame(control, "Control")
df2<- data.frame(exposed, "Exposed")
colnames(df)<- c("distance", "treatment")
colnames(df2)<- c("distance", "treatment")
df3<-rbind(df, df2)

hist(df3$distance)
library(ggplot2)
ggplot(data = df3, aes(x=distance, group = treatment, fill = treatment)) +
  geom_density(alpha=0.4)+
  theme_bw()

t.test(distance ~ treatment, data = df3)
# mean of control = 0.4506151
# mean of exposed = 0.4354719
# p = 0.3969


compare low to high pathogen load

low<- c(0.327987780241737, 0.599062973095681, 0.52891657754545, 0.468482976276331, 0.459243270274311,
        0.36776174192756, 0.460634989182071, 0.617539958953279, 0.379431080583626, 0.571086274899624, 0.593356667797448,
        0.446608300385498, 0.50758529526185, 0.477110573096593, 0.312051534693181, 0.607359783930601,
        0.320689965999387, 0.216547007319426, 0.46153506153269, 0.608119794380021, 0.556408784205222,
        0.322590873017617, 0.306833992898455, 0.366242458458616, 0.503278680884526, 0.651137279845216, 0.65449292529986,
        0.276853205792788, 0.357177229127925, 0.447842477424673, 0.50964976534647, 0.28517988301359, 0.422542354008808,
        0.245709747443765, 0.305992821875149, 0.566525696566597, 0.314033463851792, 0.489241462233394,
        0.429134648897105, 0.502777085240355, 0.583855961456432, 0.461167610095901, 0.61487218139559, 0.577941083730973,
        0.274571524173659, 0.244710461473932, 0.254643203645604, 0.440885892849073, 0.599496598762292, 0.575158141343031
)
high <- c(0.247782390664317, 0.322683580696375, 0.537250116385777, 0.570978498312689, 0.532341699014778,
          0.348657205601067, 0.408711445143216, 0.277584197536296, 0.437896470801139,
          0.353855250183178, 0.247102049304644, 0.321059012080125, 0.48174041890816, 0.478837851176654, 0.398234091707475,
          0.507523269792513, 0.495008100773835, 0.32251083251387, 0.547670861461372, 0.319795374506177, 0.333059920046721,
          0.406858226783548, 0.23452884469437, 0.617359322273927, 0.406029857071407, 0.274454719682811,
          0.336297611826307, 0.297417029398962, 0.287905532038042, 0.55641842183168, 0.434599706491461, 0.533519777475217,
          0.320129163575081, 0.460330549507781, 0.353659854815383, 0.465996323680743, 0.404304526321923, 0.478660267878825,
          0.441650465204872, 0.461200478879983, 0.219835169071612, 0.709056604939942, 0.642645225028707, 0.62915737076822
)
df<- data.frame(low, "Low")
df2<- data.frame(high, "High")
df4<- data.frame(control, "Control")
colnames(df)<- c("distance", "load")
colnames(df2)<- c("distance", "load")
colnames(df4)<- c("distance", "load")
df3<-rbind(df, df2, df4)

hist(df3$distance)
library(ggplot2)
ggplot(data = df3, aes(x=distance, group = load, fill = load)) +
  geom_density(alpha=0.6)+
  theme_bw()

# Part 3: alpha diversity comparison between exposed and unexposed abalone at end ----
setwd("~/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS")
library(qiime2R)
physeq<-qza_to_phyloseq(features="FECAL/table-deblur-fecal-survived.qza", 
                        tree="insertion-tree-deblur.qza", 
                        taxonomy="silva-taxonomy-deblur.qza", 
                        metadata= "abalone_dna_metadata_filtered.tsv")
library(phyloseq)
library(ggplot2)
dFecal11<-subset_samples(physeq, Time == "11")
p <- plot_richness(dFecal11, color = "Treatment", x = "Treatment", measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"))
p <- p + geom_boxplot(aes(fill = Treatment), alpha=0.2)
plot(p)

# Panel D: Diff abundance analysis over time and between conditions ----

# Over Time (physeq phyloseq object has all fecal samples from both treatments from only survived abalone)

# Maaslin2
library(Maaslin2)
library(mia)
library(qiime2R)
library(phyloseq)
setwd("~/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS")
physeq<-qza_to_phyloseq(features="FECAL/table-deblur-fecal-survived.qza", 
                        tree="insertion-tree-deblur.qza", 
                        taxonomy="silva-taxonomy-deblur.qza", 
                        metadata= "abalone_dna_metadata_filtered.tsv")
# convert phyloseq to TSE
fecal_tse <- makeTreeSummarizedExperimentFromPhyloseq(physeq)
# set seed to avoid random variance in some tools
set.seed(13253)
# maaslin expects features as columns and samples as rows 
# for both the asv/otu table as well as meta data 
asv <- t(assay(fecal_tse))
meta_data <- data.frame(colData(fecal_tse))

Maaslin2(input_data = asv, 
         input_metadata = meta_data, 
         analysis_method = "LM", #default - when working with count data that is transformed by TSS, CPLM has not shown performance exceeding LM, and I am filtering out rare taxa to deal with some sparcity in the data
         normalization = "TSS", # default - this divides the count per OTU by the total counts in each sample as a normalization, effectively creating relative abundances
         correction = "BH", #default - correction for multiple hypothesis testing
         min_prevalence = 0.85, #10-90% is reasonable, ~200 features are present across 85% of my samples
         min_abundance = 3, # my input data is n counts, so this is minimum number of counts per each sequence
         transform = "LOG", # default and my preferred transformation for interpretation
         output = "outputTIME", # output directory 
         fixed_effects = c("Time"), # I am just modeling all samples (either treatment) over time because no significant differences were found btn treatments across all timepoints
         reference = c("Time, 0"), # 0 is the reference because this is pre-exposure to header tank abalone
         random_effects = c("tube_id")) # the tube id is the individual animal identifier to account for repeated sampling of same animals

# OUTPUT GOES DIRECTLY TO FOLDER

# Between conditions at T=11 (dFecal11 phyloseq object has just fecal samples from T11)

dFecal11<-subset_samples(physeq, Time == "11")
fecal11_tse <- makeTreeSummarizedExperimentFromPhyloseq(dFecal11)
asv <- t(assay(fecal11_tse))
meta_data <- data.frame(colData(fecal11_tse))

Maaslin2(input_data = asv, 
         input_metadata = meta_data, 
         analysis_method = "LM",
         normalization = "TSS",
         correction = "BH",
         min_prevalence = 0.85,
         min_abundance = 3, 
         transform = "LOG",
         output = "outputTREATMENT", 
         fixed_effects = c("Treatment"), # here I am comparing the difference between exposed or control at only the final timepoint
         reference = c("Treatment, Negative")) # control abalone are the reference

#NOTE: The heatmap will not be generated if you only use one fixed effect while using the model

# Run Model on BOth Time and Exposure

asv <- t(assay(fecal_tse))
meta_data <- data.frame(colData(fecal_tse))

Maaslin2(input_data = asv, 
         input_metadata = meta_data, 
         analysis_method = "LM", #default - when working with count data that is transformed by TSS, CPLM has not shown performance exceeding LM, and I am filtering out rare taxa to deal with some sparcity in the data
         normalization = "TSS", # default - this divides the count per OTU by the total counts in each sample as a normalization, effectively creating relative abundances
         correction = "BH", #default - correction for multiple hypothesis testing
         min_prevalence = 0.85, #10-90% is reasonable, ~200 features are present across 85% of my samples
         min_abundance = 3, # my input data is n counts, so this is minimum number of counts per each sequence
         transform = "LOG", # default and my preferred transformation for interpretation
         output = "outputTimeTreatment", # output directory 
         fixed_effects = c("Time", "Treatment"), # use 2 fixed effects
         reference = c("Time, 0", "Treatment, Negative"), # set reference for both effects
         random_effects = c("tube_id"))

# ANCOM-BC (argument for ancom-bc: https://www.nature.com/articles/s41522-020-00160-w)

# tutorial https://microbiome.github.io/OMA/differential-abundance.html#differential-abundance-analysis
library(ANCOMBC)
# Between Exposed and Control

dFecal11<-subset_samples(physeq, Time == "11")
fecal11_tse <- makeTreeSummarizedExperimentFromPhyloseq(dFecal11)

abcTREATMENT <- ancombc2(
  data = fecal11_tse,
  tax_level="Genus",
  fix_formula = "Treatment",
  p_adj_method = "fdr", 
  prv_cut = 0.5,
  lib_cut = 100, 
  group = "Treatment", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 100, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100), # use max_iter >= 100 on real data 
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)
res_Treatment = abcTREATMENT$res

# Over Time # currently NOT WORKING
fecal_tse <- makeTreeSummarizedExperimentFromPhyloseq(physeq)

abcTIME <- ancombc2(
  data = fecal_tse,
  tax_level="Genus",
  fix_formula = "Time + Treatment",
  rand_formula = "(Time | tube_id)",
  p_adj_method = "fdr", 
  prv_cut = 0.5,
  lib_cut = 100, 
  group = "Treatment", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 100, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100), # use max_iter >= 100 on real data 
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)
library(tidyverse)
res_TIME = abcTIME$res %>%
  mutate_if(is.numeric, function(x) round(x, 2))
res_TIME %>%
  datatable(caption = "ANCOM-BC2 Primary Results")

## Figure 2 ----

# Panel A: qPCR copies between tissue types and feces ----

# read in metadata with copy numbers
library(readr)
abalone_dna_metadata_filtered <- read_delim("~/d_abalone_caxc/qiime2/abalone_dna_metadata_filtered.tsv", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)
#  subset just to T=11 final samples
T11<-subset(abalone_dna_metadata_filtered, month == "Apr")

# subset just to Exposed abalone, since none of control were positive 
T11_pos<-subset(T11, Exposed == "Positive")

# plot
ggplot(T11_pos, aes(x=factor(sample_type_2, level=c('PE', 'DG', 'DI', 'GC', 'F')), copies_normalized_pseudo, fill=sample_type_2))+
  geom_boxplot(color = "black", lwd = 0.5)+
  scale_y_continuous(n.breaks = 7, trans = "log10")+
  scale_fill_brewer(palette = "Oranges")+
  labs(y = "Normalized CaXc copies/uL", x = "Sample Type")+
  theme_bw()

# run statistical test
hist(T11_pos$log10_copies_normalized)
shapiro.test(T11_pos$log10_copies_normalized) # not normally distributed
# run Kruskal Wallis
kruskal.test(log10_copies_normalized ~ sample_type_2, T11_pos) # p = 6.107e-8
library(dunn.test)
dunn.test(T11_pos$log10_copies_normalized, T11_pos$sample_type_2, method = "BH")
# DI < DG
# F = DG
# F > DI
# GC = DG
# GC > DI
# GC = F
# PE > DG
# PE >> DI
# PE > F
# PE > GC
# DI < DG = F = GC < PE

# Panel B: Alpha Diversity ----
setwd("~/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS")
library(qiime2R)
physeq<-qza_to_phyloseq(features="T11/table-deblur-T11-control.qza", 
                        tree="insertion-tree-deblur.qza", 
                        taxonomy="silva-taxonomy-deblur.qza", 
                        metadata= "abalone_dna_metadata_filtered.tsv")
library(microbiome)
pielou_control<-evenness(physeq, index = "pielou", zeroes = TRUE, detection = 0)
T11_control_pielou <- cbind(sample_data(physeq), pielou_control)
library(ggplot2)
ggplot(T11_control_pielou, aes(x =factor(sample_type_2, level=c('PE', 'DG', 'DI', 'GC', 'F')), y = pielou, fill = sample_type_2))+
  geom_boxplot(color = "black", lwd = 0.5)+
  scale_fill_brewer(palette = "Greens")+
  labs(y = "Microbiome Evenness", x = "Sample Type")+
  theme_bw()

physeq_exposed<-qza_to_phyloseq(features="T11/table-deblur-T11-exposed.qza", 
                        tree="insertion-tree-deblur.qza", 
                        taxonomy="silva-taxonomy-deblur.qza", 
                        metadata= "abalone_dna_metadata_filtered.tsv")
pielou_exposed<-evenness(physeq_exposed, index = "pielou", zeroes = TRUE, detection = 0)
T11_exposed_pielou <- cbind(sample_data(physeq_exposed), pielou_exposed)
ggplot(T11_exposed_pielou, aes(x =factor(sample_type_2, level=c('PE', 'DG', 'DI', 'GC', 'F')), y = pielou, fill = sample_type_2))+
  geom_boxplot(color = "black", lwd = 0.5)+
  scale_fill_brewer(palette = "Oranges")+
  labs(y = "Microbiome Evenness", x = "Sample Type")+
  theme_bw()

# Final Panel ----

# Maaslin2 to compare tissue to fecal samples

setwd("~/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS")
library(qiime2R)
physeq<-qza_to_phyloseq(features="T11/table-deblur-T11.qza", 
                        tree="insertion-tree-deblur.qza", 
                        taxonomy="silva-taxonomy-deblur.qza", 
                        metadata= "abalone_dna_metadata_filtered.tsv")
library(phyloseq)
tissuefecal_physeq<-subset_samples(physeq, sample_type == c("fecal","tissue"))
tissuefecal_physeq_f<-subset_samples(physeq, sampleID == c("fecal","tissue"))
View(sample_data(tissuefecal_physeq))
library(mia)
tse <- makeTreeSummarizedExperimentFromPhyloseq(tissuefecal_physeq)
asv <- t(assay(tse))
meta_data <- data.frame(colData(tse))
library(Maaslin2)
Maaslin2(input_data = asv, 
         input_metadata = meta_data, 
         analysis_method = "LM",
         normalization = "TSS",
         correction = "BH",
         min_prevalence = 0.85,
         min_abundance = 3, 
         transform = "LOG",
         output = "outputSampleType", 
         fixed_effects = c("Treatment", "sample_type"), # include both treatment and sample type to get heatmap
         reference = c("Treatment, Negative", "sample_type, fecal")) # control always reference, using fecal as reference to see what is more present in tissue

# need to rename OTUIDs as taxonomic groups - merge OTUIDs with silva taxonomy by left join
library(readr)
OTUs <- read_csv("/Users/emilykunselman/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS/outputSampleType/Maaslin2_OTUIDs_sampleType_Treatment.csv")
silva_tax <- read_csv("/Users/emilykunselman/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS/exported-silva-taxonomy-deblur/taxonomy_edited.csv")
#not to unlibrary mia
detach("package:mia", unload = TRUE)
library(dplyr)
maaslin2_OTUs<-left_join(OTUs, silva_tax)
write_csv(maaslin2_OTUs, "/Users/emilykunselman/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS/outputSampleType/maaslin2_OTUs.csv")

# Figure 3 ----

# Taxa Bar Plot ----

library(tidyverse)
library(qiime2R)
metadata<-readr::read_tsv("/Users/emilykunselman/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS/abalone_dna_metadata_filtered.tsv")
metadata<-as.data.frame(metadata)
metadata<- filter(metadata, Time == "11")
metadata<- filter(metadata, sample_type != "fecal")
table<-read_qza("/Users/emilykunselman/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS/T11/table-deblur-T11-tissue.qza")$data
taxonomy<-read_qza("/Users/emilykunselman/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS/silva-taxonomy-deblur.qza")$data %>% parse_taxonomy()
taxasums<-summarize_taxa(table, taxonomy)$Genus
taxa_barplot(taxasums, metadata, "sample_type")
ggsave("taxa_barplot.pdf", height = 4, width = 15, device = "pdf")


# what about merging all samples from the same sample type into one average?
setwd("~/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS")
physeq<-qza_to_phyloseq(features="T11/table-deblur-T11-tissue.qza", 
                        tree="insertion-tree-deblur.qza", 
                        taxonomy="silva-taxonomy-deblur.qza", 
                        metadata= "abalone_dna_metadata_filtered.tsv")
library(phyloseq)
str(sample_data(physeq))
# Want to MERGE SAMPLE TYPES
# first, convert all sample data to factors
df <- as.data.frame(lapply(sample_data(physeq),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq)
sample_data(physeq) <- sample_data(df)
# then merge samples
merged_replicates = merge_samples(physeq,"sample_type_2")  # abundances are summed
# Want to agglomerate taxa to remove uncultured and NAs
physeq1<-merged_replicates %>% tax_fix(unknowns = c("uncultured"))
physeq2<-physeq1 %>% tax_fix() %>% tax_agg(rank = "Genus")
# Want to make taxa bar plot 
library(microViz)
sample_names(physeq2) # use these to order taxa bar plot
comp_barplot(physeq2,
             tax_level = "Genus", n_taxa = 10,
             bar_outline_colour = NA,
             sample_order = c("PE", "DG", "DI", "GC"),
             bar_width = 1,
             taxon_renamer = toupper
)
comp_barplot(physeq2,
             tax_level = "Genus", n_taxa = 10,
             bar_outline_colour = "black",
             sample_order = c("PE", "DG", "DI", "GC"),
             bar_width = 0.9,
             taxon_renamer = toupper
)

# do this again but for control vs exposed
# merge samples
merged_replicates2 = merge_samples(physeq,"Treatment") # abundances are summed
# Want to agglomerate taxa to remove uncultured and NAs
physeq3<-merged_replicates2 %>% tax_fix(unknowns = c("uncultured"))
physeq4<-physeq3 %>% tax_fix() %>% tax_agg(rank = "Genus")
# Want to make taxa bar plot 
library(microViz)
sample_names(physeq4) # use these to order taxa bar plot
#plot
comp_barplot(physeq4,
             tax_level = "Genus", n_taxa = 10,
             bar_outline_colour = "black",
             sample_order = "asis",
             bar_width = 0.9,
             taxon_renamer = toupper
)

# A - Plot Taxa relative abundances per sample type between control and exposed
#Post esophagus
setwd("~/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS")
library(microViz)
physeqPE<-qza_to_phyloseq(features="T11/table-deblur-T11-PE.qza", 
                        tree="insertion-tree-deblur.qza", 
                        taxonomy="silva-taxonomy-deblur.qza", 
                        metadata= "abalone_dna_metadata_filtered.tsv")
# merge samples
merged_replicatesPE = merge_samples(physeqPE,"Treatment")
# Want to agglomerate taxa to remove uncultured and NAs
physeqPE<-merged_replicatesPE %>% tax_fix(unknowns = c("uncultured"))
physeqPE<-physeqPE %>% tax_fix() %>% tax_agg(rank = "Genus")
# Want to make taxa bar plot 
sample_names(physeqPE) # use these to order taxa bar plot
#plot
comp_barplot(physeqPE,
             tax_level = "Genus", n_taxa = 10,
             bar_outline_colour = "black",
             sample_order = "asis",
             bar_width = 0.9,
             taxon_renamer = toupper
)

#Digestive Gland
setwd("~/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS")
library(microViz)
physeqDG<-qza_to_phyloseq(features="T11/table-deblur-T11-DG.qza", 
                          tree="insertion-tree-deblur.qza", 
                          taxonomy="silva-taxonomy-deblur.qza", 
                          metadata= "abalone_dna_metadata_filtered.tsv")
# merge samples
merged_replicatesDG = merge_samples(physeqDG,"Treatment")
# Want to agglomerate taxa to remove uncultured and NAs
physeqDG<-merged_replicatesDG %>% tax_fix(unknowns = c("uncultured"))
physeqDG<-physeqDG %>% tax_fix() %>% tax_agg(rank = "Genus")
# Want to make taxa bar plot 
sample_names(physeqDG) # use these to order taxa bar plot
#plot
comp_barplot(physeqDG,
             tax_level = "Genus", n_taxa = 10,
             bar_outline_colour = "black",
             sample_order = "asis",
             bar_width = 0.9,
             taxon_renamer = toupper
)

#Distal Instestine
setwd("~/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS")
library(microViz)
physeqDI<-qza_to_phyloseq(features="T11/table-deblur-T11-DI.qza", 
                          tree="insertion-tree-deblur.qza", 
                          taxonomy="silva-taxonomy-deblur.qza", 
                          metadata= "abalone_dna_metadata_filtered.tsv")
# merge samples
merged_replicatesDI = merge_samples(physeqDI,"Treatment")
# Want to agglomerate taxa to remove uncultured and NAs
physeqDI<-merged_replicatesDI %>% tax_fix(unknowns = c("uncultured"))
physeqDI<-physeqDI %>% tax_fix() %>% tax_agg(rank = "Genus")
# Want to make taxa bar plot 
sample_names(physeqDI) # use these to order taxa bar plot
#plot
comp_barplot(physeqDI,
             tax_level = "Genus", n_taxa = 10,
             bar_outline_colour = "black",
             sample_order = "asis",
             bar_width = 0.9,
             taxon_renamer = toupper
)

#Gut Contents
setwd("~/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS")
library(microViz)
physeqGC<-qza_to_phyloseq(features="T11/table-deblur-T11-GC.qza", 
                          tree="insertion-tree-deblur.qza", 
                          taxonomy="silva-taxonomy-deblur.qza", 
                          metadata= "abalone_dna_metadata_filtered.tsv")
# merge samples
merged_replicatesGC = merge_samples(physeqGC,"Treatment")
# Want to agglomerate taxa to remove uncultured and NAs
physeqGC<-merged_replicatesGC %>% tax_fix(unknowns = c("uncultured"))
physeqGC<-physeqGC %>% tax_fix() %>% tax_agg(rank = "Genus")
# Want to make taxa bar plot 
sample_names(physeqGC) # use these to order taxa bar plot
#plot
comp_barplot(physeqGC,
             tax_level = "Genus", n_taxa = 10,
             bar_outline_colour = "black",
             sample_order = "asis",
             bar_width = 0.9,
             taxon_renamer = toupper
)

# Panel ?: Differential Abundance between Treatment and Control ----

library(qiime2R)
setwd("~/d_abalone_caxc/qiime2/DEBLUR_ANALYSIS")
physeqPE<-qza_to_phyloseq(features="T11/table-deblur-T11-PE.qza", 
                        tree="insertion-tree-deblur.qza", 
                        taxonomy="silva-taxonomy-deblur.qza", 
                        metadata= "abalone_dna_metadata_filtered.tsv")
physeqDI<-qza_to_phyloseq(features="T11/table-deblur-T11-DI.qza", 
                          tree="insertion-tree-deblur.qza", 
                          taxonomy="silva-taxonomy-deblur.qza", 
                          metadata= "abalone_dna_metadata_filtered.tsv")
# Maaslin 2 to determine DA
library(mia)
tse <- makeTreeSummarizedExperimentFromPhyloseq(physeqPE)
asv <- t(assay(tse))
meta_data <- data.frame(colData(tse))
library(Maaslin2)
# decrease prevalence filter to 0.1 because looking at rarer OTUs also
Maaslin2(input_data = asv, 
         input_metadata = meta_data, 
         analysis_method = "LM",
         normalization = "TSS",
         correction = "BH",
         min_prevalence = 0.1,
         min_abundance = 3, 
         transform = "LOG",
         output = "outputPE", 
         fixed_effects = c("Treatment"),
         reference = c("Treatment, Negative"))

tse <- makeTreeSummarizedExperimentFromPhyloseq(physeqDI)
asv <- t(assay(tse))
meta_data <- data.frame(colData(tse))
# decrease prevalence filter to 0.1 because looking at rarer OTUs also
Maaslin2(input_data = asv, 
         input_metadata = meta_data, 
         analysis_method = "LM",
         normalization = "TSS",
         correction = "BH",
         min_prevalence = 0.1,
         min_abundance = 3, 
         transform = "LOG",
         output = "outputDI", 
         fixed_effects = c("Treatment"),
         reference = c("Treatment, Negative"))


# ANCOM-BC to determine DA
library(ANCOMBC)
library(knitr)
library(tidyverse)
library(ggplot2)
tse <- makeTreeSummarizedExperimentFromPhyloseq(physeqPE)
asv <- t(assay(tse))
meta_data <- data.frame(colData(tse))
out <- ancombc2(
  data = tse,
  tax_level="Genus",
  fix_formula = "Treatment", 
  p_adj_method = "BH", 
  prv_cut = 0.1,
  lib_cut = 1000, 
  group = "Treatment", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 1000, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 1000), # use max_iter >= 100 on real data 
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)


res_PE <- out$res #get results

res_PE_TRUE<-filter(res_PE, diff_TreatmentPositive == TRUE) #filter to significant results
df_fig = res_PE_TRUE %>% 
  dplyr::arrange(desc(lfc_TreatmentPositive)) %>%
  dplyr::mutate(direct = ifelse(lfc_TreatmentPositive> 0, "Exposed", "Control"))

df_fig$direct = factor(df_fig$direct, 
                       levels = c("Exposed", "Control"))
#plot
ggplot(df_fig, aes(x = reorder(taxon, -lfc_TreatmentPositive), y = lfc_TreatmentPositive, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_TreatmentPositive - se_TreatmentPositive, ymax = lfc_TreatmentPositive + se_TreatmentPositive), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes between Exposed and Control in PE") + 
  scale_fill_manual(values = c("darkorange2", "darkcyan")) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("ancombcPE.pdf", height = 6, width = 3, device = "pdf")

tse <- makeTreeSummarizedExperimentFromPhyloseq(physeqDI)
asv <- t(assay(tse))
meta_data <- data.frame(colData(tse))
out2 <- ancombc2(
  data = tse,
  tax_level="Genus",
  fix_formula = "Treatment", 
  p_adj_method = "BH", 
  prv_cut = 0.1,
  lib_cut = 1000, 
  group = "Treatment", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 1000, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 1000), # use max_iter >= 100 on real data 
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)

res_DI <- out2$res #get results

res_DI_TRUE<-filter(res_DI, diff_TreatmentPositive == TRUE) #filter to significant results
df_fig = res_DI_TRUE %>% 
  dplyr::arrange(desc(lfc_TreatmentPositive)) %>%
  dplyr::mutate(direct = ifelse(lfc_TreatmentPositive> 0, "Exposed", "Control"))

df_fig$direct = factor(df_fig$direct, 
                       levels = c("Exposed", "Control"))
#plot
ggplot(df_fig, aes(x = reorder(taxon, -lfc_TreatmentPositive), y = lfc_TreatmentPositive, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_TreatmentPositive - se_TreatmentPositive, ymax = lfc_TreatmentPositive + se_TreatmentPositive), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes between Exposed and Control in DI") + 
  scale_fill_manual(values = c("darkorange2", "darkcyan")) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("ancombcDI.pdf", height = 6, width = 3, device = "pdf")
