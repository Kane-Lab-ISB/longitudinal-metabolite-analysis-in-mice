library(ggplot2)
library(dplyr)
library(ggfortify)
library(ggpubr)
library(pheatmap)
library(tidyr)
library(readxl)
library(tidyverse)
library(stringr)
library(caret)
library(limma)
library(edgeR)
library(ComplexUpset)
library(GGally)
library(viridis)
library(UpSetR)
library(ComplexHeatmap)
library(splines)
library(glmnet)
library(lme4)
library(nlme)
library(lmerTest)
library(FELLA)
library(WGCNA)
library(rmcorr)
library(reshape2)
library(limmaDE2)
library(ggnewscale)
library(org.Mm.eg.db)
library(KEGGREST)
library(magrittr)
library(e1071)
library(randomForest)

# 1. import metabolomics, meta-data, and some useful database
source("metabolite_functions.R")
load("discovery_cohort_data.RData") # discovery cohort data, including metabolites and metadata
load("validation_dat.RData")
load("metabolite_names.RData") # metabolite names match
load("fella.RData") # fella enrichement analysis database

## 1a. PCA plot (suppl_fig.1a)
ggplot2::theme_set(theme_bw() + 
                     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           plot.title = element_text(hjust = 0.5), text = element_text(size = 16), 
                           panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                           axis.text.y = element_text(color = "black"),
                           axis.text.x = element_text(color = "black")))
# exclude metabolites with less variants
dis_meta_dat = dis_dat[, -c(129:942)]
dis_quant_dat = dis_dat[, c(129:942)]
rownames(dis_quant_dat) = dis_meta_dat$id
dis_quant_dat_keep = dis_quant_dat[,-nearZeroVar(dis_quant_dat)]
dis_metab_pca = stats::prcomp(dis_quant_dat_keep)
ggplot(data = dis_metab_pca$x,
       aes(x = PC1,
           y = PC2,
           color = dis_meta_dat$Time.label,
           shape = dis_meta_dat$Sex)) +
  geom_point(size = 3) + 
  scale_color_manual(values = c("#F1E51DFF", "#BADE28FF", "#5BC863FF", "#21A685FF", "#2F6B8EFF")) +
  labs(x = "PC1", y = "PC2", 
       color = "Time points", 
       shape = "Sex")

## 1b. PC associations with factors (suppl_fig.1b)
p_mat_cat = apply(dis_meta_dat[, c("Mouse.ID", "Time.label", "Sex", "cage")], 2, function(factor){
  factor = as.factor(factor)
  p = apply(dis_metab_pca$x[, 1:10], 2, function(x) {
    -log10(summary(aov(lm(x ~ factor)))[[1]][1, "Pr(>F)"])
  })
  p
})
p_mat_cts = apply(dis_meta_dat[, c("Age.at.assesment..days.", "Age.at.Death..days.", "Time.to.death..days.")], 2, function(factor){
  factor = as.numeric(factor)
  p = apply(dis_metab_pca$x[, 1:10], 2, function(x) {
    -log10(summary(aov(lm(x ~ factor)))[[1]][1, "Pr(>F)"])
  })
  p
})
p_mat = cbind(p_mat_cat, p_mat_cts) %>% t() %>% as.matrix()
pheatmap(p_mat,
         border_color = "white",
         color = colorRampPalette(c("grey90", "#FDE725FF", "#2F6B8EFF"))(50),
         cluster_cols = F,
         angle_col = "45",
         show_rownames = TRUE)

# 2. differentially abundant metabolites (DAM) analysis with time splines 
dis_meta_dat_excl = dis_meta_dat %>%
  filter(Time.label != "T5")
time_point_excl = as.numeric(factor(dis_meta_dat_excl$Time.label))
X = ns(time_point_excl, df = 3)
sex = factor(dis_meta_dat_excl$Sex)
design_ns = model.matrix(~ 0 + sex + sex:X)
colnames(design_ns) = make.names(colnames(design_ns))
dis_quant_dat_excl = dis_quant_dat_keep %>%
  filter(rownames(.) %in% dis_meta_dat_excl$id)
counts_excl = dis_quant_dat_excl %>% t() %>% as.data.frame() 
corfit = duplicateCorrelation(counts_excl, design_ns, block = dis_meta_dat_excl$Mouse.ID)
fit_ns = lmFit(counts_excl, design_ns, block = dis_meta_dat_excl$Mouse.ID, correlation = corfit$consensus)
fit2_ns = eBayes(fit_ns, trend = TRUE, robust = TRUE)
cont_mix = makeContrasts((sexfemale.X3 + sexmale.X3)/2 - (sexfemale.X2 + sexmale.X2)/2,
                         (sexfemale.X2 + sexmale.X2)/2 - (sexfemale.X1 + sexmale.X1)/2,
                         levels = colnames(coef(fit2_ns)))
cont_female = makeContrasts(sexfemale.X2 - sexfemale.X1,
                            sexfemale.X3 - sexfemale.X2,
                            levels = colnames(coef(fit2_ns)))
cont_male = makeContrasts(sexmale.X2 - sexmale.X1,
                          sexmale.X3 - sexmale.X2,
                          levels = colnames(coef(fit2_ns)))
cont_sex = makeContrasts(sexfemale.X1 - sexmale.X1,
                         sexfemale.X2 - sexmale.X2,
                         sexfemale.X3 - sexmale.X3,
                         levels = colnames(coef(fit2_ns)))

cont_lst = list(cont_mix, cont_female, cont_male, cont_sex)
tt_global_lst = lapply(cont_lst, function(x) {
  fit_cont = contrasts.fit(fit2_ns, x)
  fit_cont = eBayes(fit_cont)
  tt_cont = topTable(fit_cont, n = Inf, adjust.method = "BH") %>%
    filter(P.Value < 0.05/781)
  return(tt_cont)
})

## 2a. mix female and male expression heatmap (suppl_fig.2a)
mix_sex_order = dis_dat %>%
  arrange(., Sex, `Time.label`) 
mix_sex_expression = mix_sex_order %>%
  dplyr::select(unlist(rownames(tt_global_lst[[1]]))) %>%
  t() %>% as.matrix()
ComplexHeatmap::pheatmap(mix_sex_expression,
                         color = colorRampPalette(c("#FDE725FF", "grey95", "#46307EFF"))(20),
                         scale = "row",
                         cluster_cols = F,
                         annotation_col = data.frame(`Time_point` = mix_sex_order$Time.label,
                                                     Sex = mix_sex_order$Sex),
                         annotation_colors = list(`Time_point` = c(BL = "#F1E51DFF", T2 = "#BADE28FF", T3 = "#5BC863FF", T4 = "#21A685FF", T5 = "#2F6B8EFF"),
                                                  Sex = c(female = "#D1E11CFF", male = "#443983FF")),
                         show_rownames = FALSE)

## 2b. co-expression network analysis on 527 metabolites (suppl_fig.2b)
mix_metabolite_expr = dis_dat %>% 
  dplyr::select(unlist(rownames(tt_global_lst[[1]])))
cor <- WGCNA::cor # bug with WGCNA, need to change to cor from WGCNA package
net = blockwiseModules(mix_metabolite_expr, power = 9, maxBlockSize = 800, ## pickSoftThreshold 9
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.3,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "mix_dem_aging",
                       verbose = 3)
cor <- stats::cor # convert cor back to stat package

mergedColors = labels2colors(net$colors, colorSeq = c("#404688FF", "#D1E11CFF", "#481B6DFF"))
datME = moduleEigengenes(mix_metabolite_expr, mergedColors)$eigengenes %>%
  as.data.frame() 
dis_dat_ME = cbind(dis_dat, datME) %>%
  as.data.frame()
plotDendroAndColors(dendro = net$dendrograms[[1]], 
                    colors = mergedColors,
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    groupLabels = TRUE)

## 2c. module eigenvalue associations
eigen_assoc_res = lapply(943:945, function(x) {
  lmer(`Age.at.assesment..days.` ~ dis_dat_ME[, x] + as.factor(Sex)  + (1 |`Mouse.ID`), 
       data = dis_dat_ME) %>% summary()
})

## 2d. hub metabolite metabolic data enrichment (suppl_fig.2c and d) 
metab_set = cbind(color = mergedColors,
                  gene = colnames(mix_metabolite_expr)) %>%
  as.data.frame()
subset1 = metab_set %>% filter(color == "#404688FF") %>% mutate(metab = gene)
subset2 = metab_set %>% filter(color == "#D1E11CFF") %>% mutate(metab = gene)
age_subset_lst = list(subset1$metab, subset2$metab)
age_subset_kegg_res = lapply(age_subset_lst, function(x) {
  kegg_res = get_kegg_res(x)
  return(kegg_res)
})
age_subset_kegg_combine = Reduce(rbind, list(tail(age_subset_kegg_res[[1]], 20) %>% mutate(level = "subset1"),
                                             tail(age_subset_kegg_res[[2]], 20) %>% mutate(level = "subset2")))
age_subset_kegg_combine$kegg_path = factor(age_subset_kegg_combine$kegg_path, levels = age_subset_kegg_combine$kegg_path)
ggplot(data = age_subset_kegg_combine,
       aes(x = as.numeric(CompoundHits),
           y = kegg_path,
           color = -log10(p.value),
           size = ratio)) +
  geom_point() +
  scale_color_viridis(direction = -1) +
  facet_wrap(. ~ as.factor(level), scales = "free") + 
  labs(x = "Number of metabolite hits",
       y = "",
       size = "Metabolite ratio")

## 2e. metabolite abundance for hub metabolites (fig.2a)
age_metabolites = dis_dat %>%
  dplyr::select(unlist(subset1$metab), unlist(subset2$metab), 
                `Mouse.ID`, id, `Age.at.assesment..days.`, Sex) %>%
  melt(., id = c("Mouse.ID", "id", "Age.at.assesment..days.", "Sex")) %>%
  as.data.frame() %>%
  dplyr::group_by(Sex, variable) %>%
  dplyr::mutate(scaled = scale(value)) %>%
  ungroup() %>%
  as.data.frame() 

age_metabolites_mean = age_metabolites %>%
  dplyr::group_by(Sex, variable, `Age.at.assesment..days.`) %>%
  summarise_at("scaled", mean) %>%
  as.data.frame() %>%
  mutate(subset = ifelse(variable %in% subset1$metab, "subset1", "subset2"))

ggplot(data = age_metabolites_mean,
       aes(x = as.numeric(`Age.at.assesment..days.`),
           y = as.numeric(scaled),
           color = as.factor(subset),
           linetype = as.factor(Sex),
           group = paste(Sex, variable))) +
  geom_line(alpha = 0.1) +
  stat_summary(aes(color = as.factor(subset),
                   linetype = as.factor(Sex),
                   group = paste(Sex, subset)),
               fun = mean, geom = "line", size = 2) +
  scale_color_manual(values = c("#404688FF", "#D1E11CFF")) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  labs(x = "Age at asessment (days)",
       y = "Metabolite abundance",
       color = "Metabolite subset",
       linetype = "Sex")

# 2f. hub metabolite selection (suppl_fig.2e)
age_at_assess = dis_dat$`Age.at.assesment..days.`
hub_metabs = lapply(list(subset1, subset2), function(x) {
  metabolite_list = x$gene
  geneModuleMembership = as.data.frame(cor(mix_metabolite_expr, datME, use = "p"))
  cluster_MM = geneModuleMembership[rownames(geneModuleMembership) %in% metabolite_list, ] %>%
    mutate(metab = rownames(.)) 
  geneTraitSignificance = as.data.frame(cor(mix_metabolite_expr, age_at_assess, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(mix_metabolite_expr)))
  GS_table = geneTraitSignificance %>%
    mutate(metab = rownames(.))
  cluster_GS = GS_table[GS_table$metab %in% metabolite_list, ]
  cluster_MMGS = merge(cluster_MM, cluster_GS, by = "metab")
  return(cluster_MMGS)
})
subset1_dat = hub_metabs[[1]] %>% 
  mutate(MM = `ME#404688FF`,
         type = rep("subset1", nrow(.))) %>%
  dplyr::select(metab, MM, V1, type)
subset2_dat = hub_metabs[[2]] %>% 
  mutate(MM = `ME#D1E11CFF`,
         type = rep("subset2", nrow(.))) %>%
  dplyr::select(metab, MM, V1, type)
subset_dat = rbind(subset1_dat, subset2_dat)
hub_dat = subset_dat %>%
  filter(abs(MM) > 0.8 & abs(V1) > 0.2)

hub_subset1 = hub_dat %>% filter(type == "subset1")
hub_subset2 = hub_dat %>% filter(type == "subset2")

ggplot(data = subset_dat,
       aes(x = abs(MM),
           y = abs(V1),
           color = as.factor(type))) +
  geom_point(size = 0.5) + 
  geom_point(data = hub_dat,
             aes(x = abs(MM),
                 y = abs(V1),
                 color = as.factor(type)),
             shape = 2,
             size = 1) +
  xlim(0, 1) + 
  ylim(0, 1) +
  geom_vline(xintercept = 0.8, lty = 3) +
  geom_hline(yintercept = 0.2, lty = 3) +
  scale_color_manual(values = c("#404688FF", "#D1E11CFF")) +
  labs(x = "Absolute ModuleMembership",
       y = "Absolute Significance",
       color = "Metabolite subset")

# 2g. hub metabolite enrichment (fig.2b)
hub_subsets = list(hub_subset1$metab, hub_subset2$metab)
hub_kegg_combine = lapply(hub_subsets, function(x) {
  get_kegg_res(x)
}) %>%
  Reduce(rbind, .)
hub_kegg_combine$kegg_path = factor(hub_kegg_combine$kegg_path, levels = hub_kegg_combine$kegg_path)
#pdf(file = "fig2b.pdf", height = 3, width = 7)
ggplot(data = hub_kegg_combine,
       aes(x = as.numeric(CompoundHits),
           y = kegg_path,
           color = -log10(p.value),
           size = ratio)) +
  geom_point() +
  scale_color_viridis(direction = -1) +
  xlim(3, 8) +
  labs(x = "Number of metabolite hits",
       y = "",
       size = "Metabolite ratio")

# 3. DAMs for multiple comparison groups

## 3a. upset plot including 4 comparison plots (fig.3a)
DEM_lst = list(`Mix` = rownames(tt_global_lst[[1]]),
               `Females` = rownames(tt_global_lst[[2]]),
               `Males` = rownames(tt_global_lst[[3]]),
               `Sex differences` = rownames(tt_global_lst[[4]]))
m = make_comb_mat(DEM_lst)
ComplexHeatmap::UpSet(m,
                      comb_col = "#20A387FF",
                      set_order = c("Sex differences", "Females", 
                                    "Males", "Mix"),
                      top_annotation = upset_top_annotation(m, add_numbers = TRUE))

## 3b. common metabolites expression (suppl_fig.3)
overlap = Reduce(intersect, lapply(tt_global_lst, function(x) rownames(x)))
common_metabilites_order = dis_dat %>%
  arrange(., `Time.label`, Sex) 

common_metabilites_expression = common_metabilites_order %>%
  dplyr::select(unlist(overlap)) %>%
  t() %>% as.matrix()

ComplexHeatmap::pheatmap(common_metabilites_expression,
                         color = colorRampPalette(c("#FDE725FF", "grey95", "#46307EFF"))(20),
                         scale = "row",
                         cluster_cols = F,
                         annotation_col = data.frame(Sex = common_metabilites_order$Sex,
                                                     `Time_point` = common_metabilites_order$Time.label),
                         annotation_colors = list(`Time_point` = c(BL = "#F1E51DFF", T2 = "#BADE28FF", T3 = "#5BC863FF", T4 = "#21A685FF", T5 = "#2F6B8EFF"),
                                                  Sex = c(female = "#D1E11CFF", male = "#443983FF")),
                         show_rownames = FALSE)

## 3c. kegg pathway enrichment for hub metabolites (fig.3b)
overlap_kegg = metab_name_match[metab_name_match$PLOT_NAME %in% overlap, ]
overlap_kegg_res = get_kegg_res(overlap_kegg$PLOT_NAME)
overlap_kegg_res$kegg_path = factor(overlap_kegg_res$kegg_path, levels = overlap_kegg_res$kegg_path)
ggplot(data = overlap_kegg_res,
       aes(x = as.numeric(CompoundHits),
           y = kegg_path,
           color = -log10(p.value),
           size = ratio)) +
  geom_point() +
  scale_x_continuous(breaks = c(4, 6, 8, 10), limits = c(3, 11)) + 
  scale_color_viridis(direction = -1) +
  theme(axis.text.y = element_text(size = 10)) + 
  labs(x = "Number of metabolite hits",
       y = "",
       size = "Metabolite ratio")

## 3d. female specific DAMs that change with age (fig.3c)
female_excl_over = setdiff(rownames(tt_global_lst[[2]]), overlap)
female_diff_inter = intersect(female_excl_over, rownames(tt_global_lst[[4]]))
female_diff_res = get_kegg_res(female_diff_inter)
female_diff_res = tail(female_diff_res, 15) %>%
  arrange(desc(p.value))
female_diff_res$kegg_path = factor(female_diff_res$kegg_path, levels = female_diff_res$kegg_path)
ggplot(data = female_diff_res,
       aes(x = as.numeric(CompoundHits),
           y = kegg_path,
           color = -log10(p.value),
           size = ratio)) +
  geom_point() + 
  scale_color_viridis(direction = -1) +
  theme(axis.text.y = element_text(size = 12, face = "bold")) + 
  labs(x = "Number of metabolite hits",
       y = "",
       size = "Metabolite ratio")
```

## 3e. abundance of corticosterone (suppl_fig.4)
corticosterone_plot = dis_dat %>%
  dplyr::select(id, Sex, `Age.at.assesment..days.`, corticosterone) %>%
  melt(., id = c("id", "Sex", "Age.at.assesment..days.")) %>% as.data.frame() %>%
  mutate(value = as.numeric(value))
ggplot(data = corticosterone_plot,
       aes(x = `Age.at.assesment..days.`,
           y = value,
           color = as.factor(Sex))) +
  stat_summary(fun = "median", 
               fun.min = function(z) {quantile(z, 0.25)},
               fun.max = function(z) {quantile(z, 0.75)}) +
  stat_summary(fun = "median", aes(color = as.factor(Sex)), geom = "line", size = 1.5) + 
  scale_color_manual(values = c("#D1E11CFF", "#443983FF"))

# 4. FI/devFI calculation and association 

## 4a. FI and devFI boxplot (suppl_fig.5) 
dis_dat_fi = dis_dat %>%
  mutate(fi = FI.Score) %>%
  filter(id != "60.4_T3")
time_label = c("BL", "T2", "T3", "T4")
time_fi_lst = lapply(time_label, function(x) {
  dis_meta_time = dis_dat_fi[dis_dat_fi$Time.label == x, ]
  female_time = dis_meta_time[dis_meta_time$Sex == "female", ]
  male_time = dis_meta_time[dis_meta_time$Sex == "male", ]
  female_median_fi = median(female_time$fi)
  male_median_fi = median(male_time$fi)
  female_time_fi = female_time %>%
    mutate(devfi = fi - female_median_fi)
  male_time_fi = male_time %>%
    mutate(devfi = fi - male_median_fi)
  time_fi = rbind(female_time_fi, male_time_fi)
  return(time_fi)
})

time_fi_all_exclt5 = Reduce(rbind, time_fi_lst)
time_fi_all_t5 = dis_dat_fi[dis_dat_fi$Time.label == "T5", ] %>%
  mutate(devfi = fi - median(.$fi)) 
time_fi_all = rbind(time_fi_all_exclt5, time_fi_all_t5)
time_fi_all_plot = time_fi_all %>%
  dplyr::select(id, Sex, `Time.label`, fi, devfi) %>%
  melt(., id = c("id", "Sex", "Time.label")) %>%
  mutate(fi_type = c(rep("FI", 160), 
                     rep("devFI", 160))) %>%
  as.data.frame()
time_fi_all_plot$fi_type = factor(time_fi_all_plot$fi_type, levels = c("devFI", "FI"))
ggplot(data = time_fi_all_plot,
       aes(x = paste(Sex, Time.label),
           y = value)) +
  geom_boxplot(aes(group = paste(fi_type, Sex, Time.label)), outlier.shape = NA, color = "grey70") +
  geom_jitter(aes(color = as.factor(fi_type)),
              position = position_jitterdodge(),
              alpha = 0.75) +
  scale_color_manual(values = c("#2E6D8EFF", "#470E61FF")) +
  labs(x = "", y = "FI/devFI score", color = "")

## 4b. FI and devFI feature selection 
set.seed(2024)
seeds = vector(mode = "list", length = 101)
for (i in 1:100) seeds = sample.int(5000, 100)
seeds[[101]] = sample.int(5000, 1)

male_id = (time_fi_all %>% filter(Sex == "male"))$id
female_id = (time_fi_all %>% filter(Sex == "female"))$id
fi_dat = time_fi_all %>%
  dplyr::select(unlist(colnames(dis_quant_dat_excl)), fi) %>%
  `row.names<-`(time_fi_all$id) 
devfi_dat = time_fi_all %>%
  dplyr::select(devfi, unlist(colnames(dis_quant_dat_excl))) %>%
  `row.names<-`(time_fi_all$id) %>%
  mutate(fi = devfi) %>%
  dplyr::select(-devfi)
fi_features_parameters = continuous_features(fi_dat)
devfi_features_parameters = continuous_features(devfi_dat)

## 4c. frequency of selection (suppl_fig.6a)
frequency_lst = lapply(list(fi_features_parameters,
                            devfi_features_parameters), function(x) {
                              feature_freq = Reduce(cbind, lapply(1:100, function(y) {x[[y]][[1]]}))
                              colnames(feature_freq) = paste0("run", 1:100)
                              freq_sum = rowSums(feature_freq)
                              return(freq_sum)
                            })
freq_lst_tbl = cbind(meatbolites = colnames(dis_quant_dat_excl),
                     Reduce(cbind, frequency_lst)) %>%
  as.data.frame() %>%
  `colnames<-`(c("metabolites", "FI", "devFI"))
freq_lst_melt = melt(freq_lst_tbl, id = "metabolites") %>% 
  as.data.frame() %>%
  mutate(value = as.numeric(value))
vinter = data.frame(xint = c(quantile(as.numeric(freq_lst_tbl$FI), prob = 0.8),
                             quantile(as.numeric(freq_lst_tbl$devFI), prob = 0.8)),
                    variable = c("FI", "devFI"))

ggplot(data = freq_lst_melt,
       aes(x = as.numeric(value),
           color = as.factor(variable))) +
  geom_density(size = 1.25) +
  scale_color_manual(values = c("#470E61FF", "#2E6D8EFF")) +
  xlim(0, 100) + 
  geom_vline(data = vinter, 
             aes(xintercept = xint, 
                 color = as.factor(variable)), 
             linetype = "dotted",
             size = 2) + 
  labs(x = "Frequency over 100 runs",
       y = "Density",
       color = "") 

## 4d. whole cohort features selection (fig.4b)
freq_lst_tbl = freq_lst_tbl %>% mutate(FI = as.numeric(FI), devFI = as.numeric(devFI))
FI_features = freq_lst_tbl[freq_lst_tbl$FI > quantile(freq_lst_tbl$FI, prob = 0.80), ]$metabolites 
devFI_features = freq_lst_tbl[freq_lst_tbl$devFI > quantile(freq_lst_tbl$devFI, prob = 0.80), ]$metabolites 

age_features = c(hub_subset1$metab, hub_subset2$metab) # age features
fi_age_features = intersect(age_features, FI_features) # fi_age features
fi_var_features = intersect(devFI_features, FI_features) # deltafi features
union_features = union(fi_age_features, fi_var_features)

clinical_features_lst = list(`Age_features` = age_features,
                             `FI_features` = FI_features,
                             `devFI_features` = devFI_features)

clinical_features_m = make_comb_mat(clinical_features_lst)
ComplexHeatmap::UpSet(clinical_features_m,
                      comb_col = "#20A387FF",
                      set_order = c("FI_features", "Age_features", "devFI_features"),
                      top_annotation = upset_top_annotation(clinical_features_m, add_numbers = TRUE))

## 4e. devFI unique features
devfi_uni_features = setdiff(devFI_features, union(FI_features, age_features))
devfi_uni_res = get_kegg_res(devfi_uni_features)

## 4f. FI and devFI features enrichement analysis (suppl_fig.6c and d)
FI_feature_res = get_kegg_res(FI_features)
devFI_feature_res = get_kegg_res(devFI_features)
FI_features_comb_kegg_dat = rbind(FI_feature_res %>% arrange(p.value) %>% dplyr::slice(1:15) %>% mutate(fi = rep("FI_fea", nrow(.))) , 
                                  devFI_feature_res %>% arrange(p.value) %>% dplyr::slice(1:15) %>% mutate(fi = rep("devFI_fea", nrow(.)))) %>%
  as.data.frame() %>% arrange(desc(p.value))
FI_features_comb_kegg_dat = within(FI_features_comb_kegg_dat, kegg_path <- ave(as.character(FI_features_comb_kegg_dat$kegg_path), FUN = make.unique))
FI_features_comb_kegg_dat$kegg_path = factor(FI_features_comb_kegg_dat$kegg_path, levels = FI_features_comb_kegg_dat$kegg_path)
ggplot(data = FI_features_comb_kegg_dat,
       aes(x = as.numeric(CompoundHits),
           y = kegg_path,
           color = -log10(p.value),
           size = ratio)) +
  geom_point() +
  scale_color_viridis(direction = -1) +
  labs(x = "Number of metabolite hits",
       y = "",
       color = "-log10(adjusted p-value)",
       size = "Metabolite ratio") +
  facet_wrap(.~ as.factor(fi), nrow = 2, ncol = 2, scales = "free")

## 4g. union feature enrichement analysis (fig.4c)
union_features_kegg = get_kegg_res(union_features) %>% arrange(desc(p.value))
union_features_kegg$kegg_path = factor(union_features_kegg$kegg_path, levels = union_features_kegg$kegg_path)
ggplot(data = union_features_kegg,
       aes(x = as.numeric(CompoundHits),
           y = kegg_path,
           color = -log10(p.value),
           size = ratio)) +
  geom_point() +
  scale_color_viridis(direction = -1) +
  labs(x = "Number of metabolite hits",
       y = "",
       color = "-log10(adjusted p-value)",
       size = "Metabolite ratio") 

# 5. Univariate associations of metabolites with FI/devFI

## 5a. association with current FI/devFI (w/o interaction) (fig.5b)
whole_features_assoc = time_fi_all %>%
  dplyr::select(unlist(union_features), `Age.at.assesment..days.`,
                Sex, `Mouse.ID`, fi, devfi, `Time.label`, id) %>%
  mutate(age = `Age.at.assesment..days.`/1000)

fic_lmer = sapply(1:104, function(x) {
  lmer_res = lmer(fi ~ whole_features_assoc[, x] + as.numeric(age) + as.factor(Sex) + (1 |`Mouse.ID`), data = whole_features_assoc) 
  return(lmer_res)
})

devfic_lmer = sapply(1:104, function(x) {
  lmer_res = lmer(devfi ~ whole_features_assoc[, x] + as.factor(Sex) + (1 |`Mouse.ID`), data = whole_features_assoc)
  return(lmer_res)
})

fic_lmer_coeff = get_coeff(union_features, fic_lmer)
devfic_lmer_coeff = get_coeff(union_features, devfic_lmer)

fic_bound_res = base::suppressWarnings(sapply(1:104, function(x) {
  coeff_tbl1 = summary(fic_lmer[[x]])$coefficients
  coeff_tbl2 = summary(devfic_lmer[[x]])$coefficients
  bound1 = confint(fic_lmer[[x]])
  bound2 = confint(devfic_lmer[[x]])
  return(c(coeff_tbl1[2, 1], bound1[4, 1], bound1[4, 2], coeff_tbl1[2, 5],
           coeff_tbl2[2, 1], bound2[4, 1], bound2[4, 2], coeff_tbl2[2, 5]))
}))

current_overlap = intersect(fic_lmer_coeff$metab, devfic_lmer_coeff$metab)[-9]
current_overlap_dat = fic_bound_res %>% t() %>% as.data.frame() %>%
  mutate(metab = unlist(union_features)) %>%
  `colnames<-`(c("coeff1", "lb1", "ub1", "p1", 
                 "coeff2", "lb2", "ub2", "p2", "metab")) %>%
  mutate(adjp1 = p.adjust(p1, method = "BH"),
         adjp2 = p.adjust(p2, method = "BH")) %>% 
  filter(metab %in% current_overlap) 
current_overlap_dat_convert = mapply(c, current_overlap_dat[, c("coeff1", "lb1", "ub1", "p1", "metab", "adjp1")] %>% 
                                       mutate(type = "fi") %>% 
                                       arrange(as.numeric(coeff1)),
                                     current_overlap_dat[, c("coeff2", "lb2", "ub2", "p2", "metab", "adjp2")] %>% mutate(type = "devfi")) %>%
  as.data.frame() %>%
  mutate(adjp1 = as.numeric(adjp1)) %>%
  mutate(sig = case_when(
    adjp1 < 0.001 ~ "***",
    adjp1 < 0.01 ~ "**",
    adjp1 < 0.05 ~ "*",
    adjp1 < 0.1 ~ "°",
    TRUE ~ ""))

current_overlap_dat_convert$metab = factor(current_overlap_dat_convert$metab, levels = current_overlap_dat_convert$metab[1:8])
current_overlap_dat_convert$type = factor(current_overlap_dat_convert$type, levels = c("fi", "devfi"))

ggplot(data = current_overlap_dat_convert) +
  geom_hline(aes(yintercept = metab), color = "grey95", size = 3) +
  ggstance::geom_pointrangeh(
    aes(x = as.numeric(coeff1),
        y = metab,
        xmin = as.numeric(lb1),
        xmax = as.numeric(ub1),
        color = as.factor(type)),
    size = 0.35, 
    position = position_dodge(width = 0.7)) +
  scale_color_manual(name = "",
                     values = c("#470E61FF", "#2E6D8EFF")) +
  geom_vline(xintercept = 0, lty = 3) +
  labs(x = "Coefficient (95% CI)",
       y = "") +
  xlim(-0.1, 0.1) +
  geom_text(aes(y = metab,
                x = Inf,
                label = sig,
                color = as.factor(type)),
            position = position_dodge(width = 0.7),
            hjust = 1.5, vjust = 0.9, size = 5,
            show.legend = FALSE)

## 5b. association with current FI/devFI, w/ interaction term (suppl_fig.8a)
fic_inter = sapply(1:104, function(x) {
  lmer_res = lmer(fi ~ whole_features_assoc[, x]*as.numeric(age) + as.factor(Sex) + (1 |`Mouse.ID`), data = whole_features_assoc) 
  coeff_tbl = summary(lmer_res)$coefficients
  return(c(coeff_tbl[2, 1], coeff_tbl[2, 2], coeff_tbl[2, 5],
           coeff_tbl[5, 1], coeff_tbl[5, 2], coeff_tbl[5, 5]))
})

fic_inter_coeff = fic_inter %>% t() %>% as.data.frame() %>%
  mutate(metab = unlist(union_features)) %>%
  `colnames<-`(c("metab_coeff", "metab_se", "metab_p", 
                 "metabtime_coeff", "metabtime_se", "metabtime_p", "metab")) %>%
  mutate(metabtime_adjp = p.adjust(metabtime_p, method = "BH")) %>%
  filter(metabtime_adjp < 0.05)

current_fi_p1 = ggplot(data = whole_features_assoc,
                       aes(x = `leucine`,
                           y = fi,
                           color = paste(as.factor(Sex), as.factor(age)))) +
  geom_point(size = 0.8) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  scale_color_manual(values = viridis(100)[seq(10, 100, 10)])

current_fi_p2 = ggplot(data = whole_features_assoc,
                       aes(x = `N-acetylthreonine`,
                           y = fi,
                           color = paste(as.factor(Sex), as.factor(age)))) +
  geom_point(size = 0.8) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  scale_color_manual(values = viridis(100)[seq(10, 100, 10)]) +
  labs(y = "") 

current_fi_p3 = ggplot(data = whole_features_assoc,
                       aes(x = `X-25422`,
                           y = fi,
                           color = paste(as.factor(Sex), as.factor(age)))) +
  geom_point(size = 0.8) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  scale_color_manual(values = viridis(100)[seq(10, 100, 10)]) +
  labs(y = "") 

ggarrange(current_fi_p1, current_fi_p2, current_fi_p3,
          nrow = 1, ncol = 3, align = "hv",
          common.legend = TRUE)

## 5c. association with future FI of baseline metabolites (suppl_fig.8b)
tp_group = data.frame(pre = c(rep("BL", 4), rep("T2", 3), rep("T3", 2), "T4"),
                      post = c("T2", "T3", "T4", "T5", "T3", "T4", "T5", "T4", "T5", "T5"))
whole_prepost_dat = prepost_data_merge(tp_group, whole_features_assoc, union_features)
bl_whole_prepost_dat = whole_prepost_dat %>% filter(Time.label.x == "BL")

future_fi_bl_lmer = sapply(2:105, function(x) {
  lmer_res = lmer(fi.y ~ bl_whole_prepost_dat[, x] + as.numeric(age.x) + as.numeric(age_change) + as.factor(Sex.x) + (1 |`Mouse.ID`), data = bl_whole_prepost_dat) 
  return(lmer_res)
})

future_devfi_bl_lmer = sapply(2:105, function(x) {
  lmer_res = lmer(devfi.y ~ bl_whole_prepost_dat[, x] + as.numeric(age_change) + as.factor(Sex.x) + (1 |`Mouse.ID`), data = bl_whole_prepost_dat) 
  return(lmer_res)
})

future_fi_bl_lmer_coeff = get_coeff(union_features, future_fi_bl_lmer)
future_devfi_bl_lmer_coeff = get_coeff(union_features, future_devfi_bl_lmer)

akg_bl_plot_dat = bl_whole_prepost_dat %>% 
  dplyr::select(Mouse.ID, level, fi.y, devfi.y, `alpha-ketoglutarate.x`) %>%
  melt(., id = c("Mouse.ID", "level", "alpha-ketoglutarate.x")) %>% 
  as.data.frame()  

ggplot(data = akg_bl_plot_dat,
       aes(x = `alpha-ketoglutarate.x`,
           y = value,
           color = as.factor(variable))) +
  geom_point(size = 0.8) + 
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  scale_color_manual(values = c("#470E61FF", "#2E6D8EFF"))

## 5d. association with future FI/devFI change (w/o interaction) (fig.5c)
tp_group = data.frame(pre = c(rep("BL", 4), rep("T2", 3), rep("T3", 2), "T4"),
                      post = c("T2", "T3", "T4", "T5", "T3", "T4", "T5", "T4", "T5", "T5"))
whole_prepost_dat = prepost_data_merge(tp_group, whole_features_assoc, union_features)

delta_fi_lmer = sapply(2:105, function(x) {
  lmer_res = lmer(fi_change ~ whole_prepost_dat[, x] + as.numeric(age.x) + as.numeric(age_change) + as.factor(Sex.x) + (1 |`Mouse.ID`), data = whole_prepost_dat) 
  return(lmer_res)
})
delta_fi_lmer_coeff = get_coeff(union_features, delta_fi_lmer) %>% arrange(desc(coeff))

delta_fi_bound = base::suppressWarnings(sapply(1:104, function(x) {
  coeff_tbl = summary(delta_fi_lmer[[x]])$coefficients
  bound = confint(delta_fi_lmer[[x]])
  return(c(coeff_tbl[2, 1], bound[4, 1], bound[4, 2], coeff_tbl[2, 5]))
})) %>% t() %>% as.data.frame() %>%
  mutate(metab = unlist(union_features)) %>%
  `colnames<-`(c("coeff", "lb", "ub", "p", "metab")) %>%
  mutate(adjp = p.adjust(p, method = "BH")) %>% 
  filter(adjp < 0.05) %>% 
  arrange(coeff) %>%
  mutate(sig = case_when(
    adjp < 0.001 ~ "***",
    adjp < 0.01 ~ "**",
    adjp < 0.05 ~ "*",
    adjp < 0.1 ~ "°",
    TRUE ~ ""))

delta_fi_bound$metab = factor(delta_fi_bound$metab, levels = delta_fi_bound$metab)
ggplot(data = delta_fi_bound) +
  geom_hline(yintercept = delta_fi_bound$metab[seq(2, 27, 2)], 
             color = "grey95", size = 3) + 
  ggstance::geom_pointrangeh(
    aes(x = as.numeric(coeff),
        y = metab,
        xmin = as.numeric(lb),
        xmax = as.numeric(ub)),
    size = 0.35, 
    position = position_dodge(width = 0.7),
    color = "#470E61FF") +
  geom_vline(xintercept = 0, lty = 3) +
  labs(x = "Coefficient (95% CI)",
       y = "") +
  xlim(-0.1, 0.15) +
  geom_text(aes(y = metab,
                x = Inf,
                label = sig),
            position = position_dodge(width = 0.7),
            hjust = 1.5, vjust = 0.9, size = 5,
            show.legend = FALSE) +
  theme(axis.text.y = element_text(size = 11, face = "bold"))

## 5e. association with FI/devFI change w/ interaction term
delta_fi_inter = sapply(2:105, function(x) {
  lmer_res = lmer(fi_change ~ whole_prepost_dat[, x]*as.numeric(age.x) + whole_prepost_dat[, x]* as.numeric(age_change) + as.factor(Sex.x) + (1 |`Mouse.ID`), data = whole_prepost_dat) 
  return(lmer_res)
})

delta_fi_inter_coeff = sapply(1:104, function(x) {
  coeff_tbl = summary(delta_fi_inter[[x]])$coefficients
  return(c(coeff_tbl[6, 1], coeff_tbl[6, 2], coeff_tbl[6, 5],
           coeff_tbl[7, 1], coeff_tbl[7, 2], coeff_tbl[7, 5]))
}) %>% t() %>% 
  as.data.frame() %>%
  mutate(metab = unlist(union_features)) %>%
  `colnames<-`(c("coeff1", "se1", "p1", "coeff2", "se2", "p2", "metab")) %>%
  mutate(adjp1 = p.adjust(p1, method = "BH"),
         adjp2 = p.adjust(p2, method = "BH")) %>% 
  filter(adjp1 < 0.05 | adjp2 < 0.05) %>%
  arrange(adjp1)
# no interaction term significant found for delta_fi_lmer_coeff$metab

## 5f. examples of metabolites associated with frailty change (suppl_fig.8c)
future_fi_p1 = ggplot(data = whole_prepost_dat,
                      aes(x = `creatinine.x`,
                          y = `fi.y`,
                          color = paste(as.factor(Sex.y), as.factor(age.y)))) +
  geom_point(size = 0.8) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  scale_color_manual(values = viridis(100)[seq(10, 100, 10)])

future_fi_p2 = ggplot(data = whole_prepost_dat,
                      aes(x = `creatinine.x`,
                          y = `fi_change`,
                          color = paste(as.factor(Sex.y), as.factor(age.y)))) +
  geom_point(size = 0.8) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  scale_color_manual(values = viridis(100)[seq(10, 100, 10)])

future_fi_p3 = ggplot(data = whole_prepost_dat,
                      aes(x = `phenyllactate (PLA).x`,
                          y = `fi.y`,
                          color = paste(as.factor(Sex.y), as.factor(age.y)))) +
  geom_point(size = 0.8) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  scale_color_manual(values = viridis(100)[seq(10, 100, 10)])

future_fi_p4 = ggplot(data = whole_prepost_dat,
                      aes(x = `phenyllactate (PLA).x`,
                          y = `fi_change`,
                          color = paste(as.factor(Sex.y), as.factor(age.y)))) +
  geom_point(size = 0.8) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  scale_color_manual(values = viridis(100)[seq(10, 100, 10)])

ggarrange(future_fi_p1, future_fi_p2, future_fi_p3, future_fi_p4, 
          nrow = 2, ncol = 2, align = "hv", common.legend = TRUE)

## 5g. combine biomarkers (fig.5d)
FI_assoc_lst = list(fic_lmer_coeff$metab[-c(26, 47)], 
                    delta_fi_lmer_coeff$metab, 
                    devfic_lmer_coeff$metab)
FI_features_union = Reduce(union, FI_assoc_lst)
all_FI = Reduce(c, FI_assoc_lst)
all_FI_df = data.frame(names = names(table(all_FI)),
                       count = unname(table(all_FI))) %>%
  arrange(desc(count.Freq))
test_metab = all_FI_df[all_FI_df$count.Freq >=2, ]$names
presence_tbl = lapply(test_metab, function(x) {
  lst = sapply(1:3, function(y){
    if (x %in% FI_assoc_lst[[y]]) {
      c = 1} else {c = -1}
  })
}) %>% Reduce(rbind, .) %>% as.data.frame() %>%
  `colnames<-`(c("FIc", "deltaFI", "devFIc")) %>%
  mutate(metab = test_metab) %>% 
  melt(.)
order = all_FI_df[all_FI_df$count.Freq >=2, ] %>% arrange(count.Freq)
presence_tbl$metab = factor(presence_tbl$metab, levels = order$names)

ggplot(data = presence_tbl,
       aes(x = variable,
           y = metab,
           color = as.factor(value))) +
  geom_hline(yintercept = unique(order$names[seq(1, 23, 2)]), color = "grey95", size = 3) + 
  geom_point(size = 5) +
  scale_color_manual(values = c("grey80", "#20A387FF")) +
  theme(legend.position="none",
        axis.text.y = element_text(size = 10, face = "bold"))

# 6. sex specific FI/deltaFI feature selection

## 6a. performance comparisons in whole cohort, females and males (suppl_fig.9a)
performance_dat = sapply(2:4, function(y){
  performance_lst = c(sapply(1:100, function(x) {fi_features_parameters[[x]][[y]]}),
                      sapply(1:100, function(x) {devfi_features_parameters[[x]][[y]]}))
})

performance_dat_melt = cbind(level = c(rep("FI", 100), rep("devFI", 100)),
                             performance_dat) %>% as.data.frame() %>%
  `colnames<-`(c("level", "Whole_cohort", "Females", "Males")) %>%
  melt(., id = "level") %>% as.data.frame()
performance_dat_melt$level = factor(performance_dat_melt$level, levels = c("FI", "devFI"))
performance_dat_melt$group = paste(performance_dat_melt$level, performance_dat_melt$variable)
performance_dat_melt$group = factor(performance_dat_melt$group, levels = c("FI Whole_cohort", "FI Females", "FI Males",
                                                                           "devFI Whole_cohort", "devFI Females", "devFI Males"))

ggplot(data = performance_dat_melt,
       aes(x = as.factor(level),
           y = as.numeric(value))) +
  geom_boxplot(aes(group = group),
               outlier.shape = NA,
               color = "grey70") +
  geom_jitter(aes(color = as.factor(variable)),
              position = position_jitterdodge(),
              alpha = 0.75,
              size = 0.75) +
  scale_color_manual(values = c("#20A387FF", "#D1E11CFF", "#443983FF")) +
  ylim(0, 1) + 
  labs(x = "", 
       y = "rsquared", 
       color = "Subgroup")

## 6b. sex specific feature selection
time_fi_female = time_fi_all %>% filter(Sex == "female")
time_fi_female_quant = time_fi_female %>% 
  dplyr::select(unlist(colnames(dis_quant_dat_excl))) ## 781 metabolites for variation test
time_fi_female_meta = time_fi_female %>% dplyr::select(-unlist(colnames(time_fi_female_quant)))
time_fi_female_keep = time_fi_female_quant[,-nearZeroVar(time_fi_female_quant)] ## 774 metabolites
time_fi_male = time_fi_all %>% filter(Sex == "male")
time_fi_male_quant = time_fi_male %>%
  dplyr::select(unlist(colnames(dis_quant_dat_excl)))
time_fi_male_meta = time_fi_male %>% dplyr::select(-unlist(colnames(time_fi_male_quant)))
time_fi_male_keep = time_fi_male_quant[,-nearZeroVar(time_fi_male_quant)] ## 767 metabolites

fi_female_score = time_fi_female_meta %>%
  dplyr::select(fi, devfi)
fi_female_tbls = lapply(1:2, function(x){
  tbl = cbind(time_fi_female_keep, 
              fi = fi_female_score[, x]) %>% as.data.frame()
})

fi_male_score = time_fi_male_meta %>%
  dplyr::select(fi, devfi)
fi_male_tbls = lapply(1:2, function(x){
  tbl = cbind(time_fi_male_keep, 
              fi = fi_male_score[, x]) %>% as.data.frame()
})

female_specific_features_sort = lapply(fi_female_tbls, function(x) {
  sex_specifc_features(x)
})
male_specific_features_sort = lapply(fi_male_tbls, function(x) {
  sex_specifc_features(x)
})

## 6c. selection of top 20% features
female_specific_frequency_lst = lapply(1:2, # female specific features presence ratio summary
                                       function(x) {
                                         feature_freq = Reduce(cbind, female_specific_features_sort[[x]][1:100])
                                         colnames(feature_freq) = paste0("run", 1:100)
                                         freq_sum = rowSums(feature_freq)
                                         return(freq_sum)
                                       })
male_specific_frequency_lst = lapply(1:2, # male specific features presence ratio summary
                                     function(x) {
                                       feature_freq = Reduce(cbind, male_specific_features_sort[[x]][1:100])
                                       colnames(feature_freq) = paste0("run", 1:100)
                                       freq_sum = rowSums(feature_freq)
                                       return(freq_sum)
                                     })
vinter_female = data.frame(xint = sapply(1:2, function(x) {
  quantile(as.numeric(female_specific_frequency_lst[[x]]), prob = 0.8)
}),
variable = c("FI", "devFI"))
vinter_male = data.frame(xint = sapply(1:2, function(x) {
  quantile(as.numeric(male_specific_frequency_lst[[x]]), prob = 0.8)
}),
variable = c("FI", "devFI"))
freq_lst_female = data.frame(metabolites = colnames(time_fi_female_keep),
                             Reduce(cbind, female_specific_frequency_lst)) %>%
  as.data.frame() %>%
  `colnames<-`(c("metabolites", "FI", "devFI")) 
freq_lst_male = data.frame(metabolites = colnames(time_fi_male_keep),
                           Reduce(cbind, male_specific_frequency_lst)) %>%
  as.data.frame() %>%
  `colnames<-`(c("metabolites", "FI", "devFI"))

## 6d. sex specific union features for FI (fig.6a)
female_FI_metabolites = freq_lst_female[freq_lst_female$FI > quantile(as.numeric(freq_lst_female$FI), prob = 0.8), ]$metabolites
female_devFI_metabolites = freq_lst_female[freq_lst_female$devFI > quantile(as.numeric(freq_lst_female$devFI), prob = 0.8), ]$metabolites

male_FI_metabolites = freq_lst_male[freq_lst_male$FI > quantile(as.numeric(freq_lst_male$FI), prob = 0.8), ]$metabolites
male_devFI_metabolites = freq_lst_male[freq_lst_male$devFI > quantile(as.numeric(freq_lst_male$devFI), prob = 0.8),]$metabolites

female_union_features = union_fi_features(female_FI_metabolites, female_devFI_metabolites) #n = 58
male_union_features = union_fi_features(male_FI_metabolites, male_devFI_metabolites) # n = 21
union_features_lst = list(`Whole_cohort` = union_features,
                          `Female` = female_union_features,
                          `Male` = male_union_features)
union_features_m = make_comb_mat(union_features_lst)
ComplexHeatmap::UpSet(union_features_m,
                      comb_col = "#20A387FF",
                      set_order = c("Whole_cohort", "Female", "Male"),
                      top_annotation = upset_top_annotation(union_features_m, add_numbers = TRUE))

## 6e. sex specific features enrichment analysis (suppl_fig.9c and d)
sex_spe_features_lst = list(female_FI_metabolites, 
                            male_devFI_metabolites)
sex_spe_features_name_lst = c("female_FI", "male_devFI")
sex_spe_features_kegg = lapply(1:2, function(x) {
  tbl = get_kegg_res(sex_spe_features_lst[[x]]) %>% arrange(p.value) %>% dplyr::slice(1:15) %>% mutate(level = rep(sex_spe_features_name_lst[[x]], nrow(.)))
  return(tbl)
}) %>% Reduce(rbind, .) %>% as.data.frame() %>%
  arrange(desc(p.value))
sex_spe_features_kegg = within(sex_spe_features_kegg, kegg_path <- ave(as.character(sex_spe_features_kegg$kegg_path), FUN = make.unique))
sex_spe_features_kegg$kegg_path = factor(sex_spe_features_kegg$kegg_path, levels = sex_spe_features_kegg$kegg_path)

ggplot(data = sex_spe_features_kegg,
       aes(x = as.numeric(CompoundHits),
           y = kegg_path,
           color = -log10(p.value),
           size = ratio)) +
  geom_point() +
  scale_color_viridis(direction = -1) +
  labs(x = "Number of metabolite hits",
       y = "",
       color = "-log10(adjusted p-value)",
       size = "Metabolite ratio") +
  facet_wrap(.~ as.factor(level), nrow = 1, ncol = 2, scales = "free")

# 7. sex specific features association

## 7a. abundance association with CURRENT in female (fig.6b)
female_features_assoc = time_fi_female %>%
  dplyr::select(unlist(female_union_features), `Age.at.assesment..days.`,
                `Mouse.ID`, fi, devfi, `Time.label`) %>%
  mutate(age = `Age.at.assesment..days.`/1000)
fic_female_lmer = sapply(1:58, function(x) {
  lmer_res = lmer(fi ~ female_features_assoc[, x] + as.numeric(age) + (1 |`Mouse.ID`), data = female_features_assoc) 
  return(lmer_res)
})

devfic_female_lmer = sapply(1:58, function(x) {
  lmer_res = lmer(devfi ~ female_features_assoc[, x] + (1 |`Mouse.ID`), data = female_features_assoc)
  return(lmer_res)
})
fic_female_coeff = get_coeff(female_union_features, fic_female_lmer)
devfic_female_coeff = get_coeff(female_union_features, devfic_female_lmer)

current_fi_female_bound = base::suppressWarnings(sapply(1:58, function(x) {
  coeff_tbl1 = summary(fic_female_lmer[[x]])$coefficients
  coeff_tbl2 = summary(devfic_female_lmer[[x]])$coefficients
  bound1 = confint(fic_female_lmer[[x]])
  bound2 = confint(devfic_female_lmer[[x]])
  return(c(coeff_tbl1[2, 1], bound1[4, 1], bound1[4, 2], coeff_tbl1[2, 5],
           coeff_tbl2[2, 1], bound2[4, 1], bound2[4, 2], coeff_tbl2[2, 5]))
}))
female_current_overlap = intersect(fic_female_coeff$metab, devfic_female_coeff$metab)
female_current_overlap_dat = current_fi_female_bound %>% t() %>% as.data.frame() %>%
  mutate(metab = unlist(female_union_features)) %>%
  `colnames<-`(c("coeff1", "lb1", "ub1", "p1", 
                 "coeff2", "lb2", "ub2", "p2", "metab")) %>%
  mutate(adjp1 = p.adjust(p1, method = "BH"),
         adjp2 = p.adjust(p2, method = "BH")) %>% 
  filter(metab %in% female_current_overlap) 
female_current_overlap_convert = mapply(c, female_current_overlap_dat[, c("coeff1", "lb1", "ub1", "p1", "metab", "adjp1")] %>% 
                                          mutate(type = "fi") %>% 
                                          arrange(as.numeric(coeff1)),
                                        female_current_overlap_dat[, c("coeff2", "lb2", "ub2", "p2", "metab", "adjp2")] %>% mutate(type = "devfi")) %>%
  as.data.frame() %>%
  mutate(adjp1 = as.numeric(adjp1)) %>%
  mutate(sig = case_when(
    adjp1 < 0.001 ~ "***",
    adjp1 < 0.01 ~ "**",
    adjp1 < 0.05 ~ "*",
    adjp1 < 0.1 ~ "°",
    TRUE ~ ""))

female_current_overlap_convert$metab = factor(female_current_overlap_convert$metab, levels = female_current_overlap_convert$metab[1:6])
female_current_overlap_convert$type = factor(female_current_overlap_convert$type, levels = c("fi", "devfi"))

ggplot(data = female_current_overlap_convert) +
  geom_hline(aes(yintercept = metab), color = "grey95", size = 3) +
  ggstance::geom_pointrangeh(
    aes(x = as.numeric(coeff1),
        y = metab,
        xmin = as.numeric(lb1),
        xmax = as.numeric(ub1),
        color = as.factor(type)),
    size = 0.35, 
    position = position_dodge(width = 0.7)) +
  scale_color_manual(name = "",
                     values = c("#470E61FF", "#2E6D8EFF")) +
  geom_vline(xintercept = 0, lty = 3) +
  labs(x = "Coefficient (95% CI)",
       y = "") +
  xlim(-0.1, 0.1) +
  geom_text(aes(y = metab,
                x = Inf,
                label = sig,
                color = as.factor(type)),
            position = position_dodge(width = 0.7),
            hjust = 1.5, vjust = 0.9, size = 5,
            show.legend = FALSE)

## 7b. abundance association with future change devFI in female (fig.6c)
tp_female_group = data.frame(pre = c(rep("BL", 3), rep("T2", 2), rep("T3", 1)),
                             post = c("T2", "T3", "T4", "T3", "T4", "T4"))
female_prepost_dat = prepost_data_merge(tp_female_group, female_features_assoc, female_union_features)

delta_devfi_female_lmer = sapply(2:59, function(x) {
  lmer_res = lmer(devfi_change ~ female_prepost_dat[, x] + as.numeric(age_change) + (1 |`Mouse.ID`), data = female_prepost_dat)
  return(lmer_res)
})
delta_devfi_female_coeff = get_coeff(female_union_features, delta_devfi_female_lmer) %>% arrange(desc(coeff))

delta_devfi_female_bound = base::suppressWarnings(sapply(1:58, function(x) {
  coeff_tbl = summary(delta_devfi_female_lmer[[x]])$coefficients
  bound = confint(delta_devfi_female_lmer[[x]])
  return(c(coeff_tbl[2, 1], bound[4, 1], bound[4, 2], coeff_tbl[2, 5]))
})) %>% t() %>% as.data.frame() %>%
  mutate(metab = unlist(female_union_features)) %>%
  `colnames<-`(c("coeff", "lb", "ub", "p", "metab")) %>%
  mutate(adjp = p.adjust(p, method = "BH")) %>% 
  filter(adjp < 0.05) %>% 
  arrange(coeff) %>%
  mutate(sig = case_when(
    adjp < 0.001 ~ "***",
    adjp < 0.01 ~ "**",
    adjp < 0.05 ~ "*",
    adjp < 0.1 ~ "°",
    TRUE ~ ""))

delta_devfi_female_bound$metab = factor(delta_devfi_female_bound$metab, levels = delta_devfi_female_bound$metab)
ggplot(data = delta_devfi_female_bound) +
  geom_hline(yintercept = delta_devfi_female_bound$metab[seq(2, 26, 2)], 
             color = "grey95", size = 3) + 
  ggstance::geom_pointrangeh(
    aes(x = as.numeric(coeff),
        y = metab,
        xmin = as.numeric(lb),
        xmax = as.numeric(ub)),
    size = 0.35, 
    position = position_dodge(width = 0.7),
    color = "#470E61FF") +
  geom_vline(xintercept = 0, lty = 3) +
  labs(x = "Coefficient (95% CI)",
       y = "") +
  xlim(-0.1, 0.15) +
  geom_text(aes(y = metab,
                x = Inf,
                label = sig),
            position = position_dodge(width = 0.7),
            hjust = 1.5, vjust = 0.9, size = 5,
            show.legend = FALSE) +
  theme(axis.text.y = element_text(size = 11, face = "bold"))

## 7c. abundance change association with current frailty outcomes
metab_devfic_female_lmer = sapply(2:59, function(x) {
  y = x + 129
  z = x + 63
  lmer_res = lmer(devfi.y ~ female_prepost_dat[, y] + female_prepost_dat[, z] + as.numeric(age_change) + (1 |`Mouse.ID`), data = female_prepost_dat)
  return(lmer_res)
})

metab_devfic_female_coeff = get_coeff(female_union_features, metab_devfic_female_lmer)

## 7d. female metabolite lists combination (suppl_fig.10a)
FI_female_assoc_lst = list(fic_female_coeff$metab, 
                           devfic_female_coeff$metab,
                           delta_devfi_female_coeff$metab,
                           metab_devfic_female_coeff$metab)
FI_female_features_union = Reduce(union, FI_female_assoc_lst)
all_female_FI = Reduce(c, FI_female_assoc_lst)
all_female_FI_df = data.frame(names = names(table(all_female_FI)),
                              count = unname(table(all_female_FI))) %>%
  arrange(desc(count.Freq))
female_test_metab = all_female_FI_df[all_female_FI_df$count.Freq >=2, ]$names
presence_tbl_female = lapply(female_test_metab, function(x) {
  lst = sapply(1:4, function(y){
    if (x %in% FI_female_assoc_lst[[y]]) {
      c = 1} else {c = -1}
  })
}) %>% Reduce(rbind, .) %>% as.data.frame() %>%
  `colnames<-`(c("FIc", "devFIc", "deltadevFI", "metdevFIc")) %>%
  mutate(metab = female_test_metab) %>% 
  melt(.)
order_female = all_female_FI_df[all_female_FI_df$count.Freq >=2, ] %>% arrange(count.Freq)
presence_tbl_female$metab = factor(presence_tbl_female$metab, levels = order_female$names)

ggplot(data = presence_tbl_female,
       aes(x = variable,
           y = metab,
           color = as.factor(value))) +
  geom_hline(yintercept = order_female$names[seq(2, 26, 2)], color = "grey95", size = 4) +
  geom_point(size = 5) +
  scale_color_manual(values = c("grey80", "#D1E11CFF")) +
  theme(legend.position="none",
        axis.text.y = element_text(size = 10, face = "bold"))

## 7e. abundance association with CURRENT in male
male_features_assoc = time_fi_male %>%
  dplyr::select(unlist(male_union_features), `Age.at.assesment..days.`,
                `Mouse.ID`, fi, devfi, `Time.label`) %>%
  mutate(age = `Age.at.assesment..days.`/1000)
fic_male_lmer = sapply(1:21, function(x) {
  lmer_res = lmer(fi ~ male_features_assoc[, x] + as.numeric(age) + (1 |`Mouse.ID`), data = male_features_assoc ) 
  return(lmer_res)
})
devfic_male_lmer = sapply(1:21, function(x) {
  lmer_res = lmer(devfi ~ male_features_assoc[, x] + (1 |`Mouse.ID`), data = male_features_assoc)
  return(lmer_res)
})
fic_male_coeff = get_coeff(male_union_features, fic_male_lmer) %>% arrange(desc(coeff))
devfic_male_coeff = get_coeff(male_union_features, devfic_male_lmer) %>% arrange(desc(coeff))

fic_male_inter = sapply(1:21, function(x) {
  lmer_res = lmer(fi ~ male_features_assoc[, x]*as.numeric(age) + (1 |`Mouse.ID`), data = male_features_assoc) 
  return(lmer_res)
})

fic_male_inter_coeff = sapply(1:21, function(x) {
  coeff_tbl = summary(fic_male_inter[[x]])$coefficients
  return(c(coeff_tbl[2, 1], coeff_tbl[2, 2],
           coeff_tbl[4, 1], coeff_tbl[4, 2], coeff_tbl[4, 5]))
}) %>% t() %>% 
  as.data.frame() %>%
  mutate(metab = unlist(male_union_features)) %>%
  `colnames<-`(c("coeff1", "se1", "coeff2", "se2", "p", "metab")) %>%
  mutate(adjp = p.adjust(p, method = "BH")) %>% 
  filter(adjp < 0.05) %>% arrange(adjp)

## 7f. male metabolite lists combination (suppl_fig.10b)
FI_male_assoc_lst = list(fic_male_remain, 
                         devfic_male_coeff$metab)
FI_male_features_union = Reduce(union, FI_male_assoc_lst)
all_male_FI = Reduce(c, FI_male_assoc_lst)
all_male_FI_df = data.frame(names = names(table(all_male_FI)),
                            count = unname(table(all_male_FI))) %>%
  arrange(desc(count.Freq))

male_test_metab = all_male_FI_df$names
presence_tbl_male = lapply(male_test_metab, function(x) {
  lst = sapply(1:2, function(y){
    if (x %in% FI_male_assoc_lst[[y]]) {
      c = 1} else {c = -1}
  })
}) %>% Reduce(rbind, .) %>% as.data.frame() %>%
  `colnames<-`(c("FIc", "devFIc")) %>%
  mutate(metab = male_test_metab) %>% 
  melt(.)
presence_tbl_male$metab = factor(presence_tbl_male$metab, levels = male_test_metab)

pdf(file = "suppl_fig9b.pdf", height = 5, width = 6)
ggplot(data = presence_tbl_male,
       aes(x = variable,
           y = metab,
           color = as.factor(value))) +
  geom_hline(yintercept = male_test_metab[seq(1, 18, 2)], color = "grey95", size = 4) +
  geom_point(size = 5) +
  scale_color_manual(values = c("grey80", "#443983FF")) +
  theme(legend.position="none",
        axis.text.y = element_text(size = 10, face = "bold"))

## 7g. ergothioneine plot (suppl_fig.7)
ergothioneine_plot = time_fi_all %>%
  dplyr::select(id, Sex, `Age.at.assesment..days.`, fi, devfi,
                ergothioneine) %>%
  dplyr::group_by(Sex) %>% 
  mutate(ergo_scale = scale(ergothioneine)) %>%
  ungroup() %>% as.data.frame() %>% dplyr::select(-ergothioneine) %>%
  melt(., id = c("id", "Sex", "Age.at.assesment..days.")) %>% as.data.frame() %>%
  mutate(value = as.numeric(value))

ggplot(data = ergothioneine_plot,
       aes(x = `Age.at.assesment..days.`,
           y = value,
           color = as.factor(variable))) +
  stat_summary(fun = "median", 
               fun.min = function(z) {quantile(z, 0.25)},
               fun.max = function(z) {quantile(z, 0.75)},
               aes(group = variable, shape = variable)) +
  stat_summary(fun = "median", aes(group = variable), geom = "line", size = 1.5) + 
  scale_color_manual(values = c("#470E61FF", "#2E6D8EFF", "grey40")) +
  facet_wrap(.~ as.factor(Sex), scales = "free")

# 8. association in validation cohort

## 8a. nmn datset association

nmn_dat_assoc = nmn_dat %>%
  mutate(sex_bin = ifelse(Sex == "female", 1, 0)) %>%
  mutate(age = `Age.at.assesment..days.`/1000) %>%
  dplyr::select(unlist(union_features), age,
                `Mouse.ID`, fi, devfi, `Time.label`, sex_bin)

fic_nmn_lmer = sapply(1:104, function(x) {
  lmer_res = lmer(fi ~ nmn_dat_assoc[, x] + as.numeric(age) + as.factor(sex_bin) + (1 |`Mouse.ID`), data = nmn_dat_assoc) 
  return(lmer_res)
})

fic_nmn_coeff = get_coeff(union_features, fic_nmn_lmer)

nmn_prepost_dat = prepost_data_merge(tp_group, nmn_dat_assoc, union_features)

delta_fi_nmn_lmer = sapply(2:105, function(x) {
  lmer_res = lmer(fi_change ~ nmn_prepost_dat[, x] + as.numeric(age.x) + as.factor(sex_bin.x) + as.numeric(age_change) + (1 |`Mouse.ID`), data = nmn_prepost_dat) 
  return(lmer_res)
})

delta_fi_nmn_coeff = get_coeff(union_features, delta_fi_nmn_lmer)

delta_fi_nmn_inter = sapply(2:105, function(x) {
  lmer_res = lmer(fi_change ~ nmn_prepost_dat[, x]*as.numeric(age.x) + nmn_prepost_dat[, x]*as.numeric(age_change) + as.factor(sex_bin.y) + (1 |`Mouse.ID`), data = nmn_prepost_dat) 
  return(lmer_res)
})

nmn_delta_fi_inter_res = sapply(1:104, function(x) {
  coeff_tbl = summary(delta_fi_nmn_inter[[x]])$coefficients
  return(c(coeff_tbl[2, 1],
           coeff_tbl[6, 1], coeff_tbl[6, 2], coeff_tbl[6, 5],
           coeff_tbl[7, 1], coeff_tbl[7, 2], coeff_tbl[7, 5]))
}) %>% t() %>% 
  as.data.frame() %>%
  mutate(metab = unlist(union_features)) %>%
  `colnames<-`(c("coeff", "coeff1", "se1", "p1", "coeff2", "se2", "p2", "metab")) %>%
  mutate(adjp1 = p.adjust(p1, method = "BH"),
         adjp2 = p.adjust(p2, method = "BH")) %>% 
  filter(adjp1 < 0.05 | adjp2 < 0.05) %>%
  arrange(adjp1)

delta_fi_nmn_remain = setdiff(delta_fi_nmn_coeff$metab,
                              intersect(delta_fi_nmn_coeff$metab,
                                        nmn_delta_fi_inter_res$metab))
delta_fi_nmn_keep = intersect(delta_fi_lmer_coeff$metab, delta_fi_nmn_remain) # 7 metabolites

## 8b. female association in nmn
nmn_dat_female_assoc = nmn_dat %>%
  filter(Sex == "female") %>%
  mutate(age = `Age.at.assesment..days.`/1000) %>%
  dplyr::select(unlist(female_union_features), age,
                `Mouse.ID`, fi, devfi, `Time.label`)

tp_female_group = data.frame(pre = c(rep("BL", 3), rep("T2", 2), rep("T3", 1)),
                             post = c("T2", "T3", "T4", "T3", "T4", "T4"))

nmn_prepost_female_dat = prepost_data_merge(tp_female_group, nmn_dat_female_assoc, female_union_features)

## 8c. male association in nmn (suppl_fig.11)
nmn_dat_male_assoc = nmn_dat %>%
  filter(Sex == "male") %>%
  mutate(age = `Age.at.assesment..days.`/1000) %>%
  dplyr::select(unlist(male_union_features), age,
                `Mouse.ID`, fi, devfi, `Time.label`)

fic_male_nmn_lmer = sapply(1:21, function(x) {
  lmer_res = lmer(fi ~ nmn_dat_male_assoc[, x] + as.numeric(age) + (1 |`Mouse.ID`), data = nmn_dat_male_assoc) 
  return(lmer_res)
})

fic_male_nmn_coeff = get_coeff(male_union_features, fic_male_nmn_lmer)

fic_male_nmn_inter = sapply(1:21, function(x) {
  lmer_res = lmer(fi ~ nmn_dat_male_assoc[, x]*as.numeric(age) + (1 |`Mouse.ID`), data = nmn_dat_male_assoc) 
  return(lmer_res)
})

fic_male_nmn_inter_coeff = sapply(1:21, function(x) {
  coeff_tbl = summary(fic_male_nmn_inter[[x]])$coefficients
  return(c(coeff_tbl[2, 1],
           coeff_tbl[4, 1], coeff_tbl[4, 2], coeff_tbl[4, 5]))
}) %>% t() %>% 
  as.data.frame() %>%
  mutate(metab = unlist(male_union_features)) %>%
  `colnames<-`(c("coeff", "coeff1", "se1", "p1", "metab")) %>%
  mutate(adjp1 = p.adjust(p1, method = "BH")) %>% 
  filter(adjp1 < 0.05)

ggplot(data = nmn_dat_male_assoc,
       aes(x = `2-hydroxydecanoate`,
           y = fi,
           color = as.factor(Time.label))) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  stat_cor(aes(label = paste(..p.label..))) +
  scale_color_manual(values = viridis(100)[c(20, 40, 60, 80, 95)])

# 9. metab clock for frailty in mice
total_fea = Reduce(union, list(union_features, female_union_features, male_union_features))
total_fea_rank = freq_lst_tbl[freq_lst_tbl$metabolites %in% total_fea, ] %>% arrange(desc(FI), desc(devFI))
total_fea = total_fea_rank$metabolites

# informativity-based selection to get results_102824_2, code within num_fea_tuning.R
match_element = c(2, 3, 4, 5, 6, 7, 8, 11, 13, 18, 20, 23, 26, 28, 29, 30, 31, 32, 33, 34)
test_df = data.frame(num = c(10, 20, 25, 30, 40, 50, 60, 63, 65, 70, 72, 75, 78, 80, 90, 100, 110, 120, 130, 139),
                     rsqured = sapply(c(1:34)[match_element], function(x) {results_102824_2[[x]][4, 1]}),
                     rmse = sapply(c(1:34)[match_element], function(x) {results_102824_2[[x]][4, 2]}))

num_tuning_p1 = ggplot(data = test_df,
                       aes(x = num,
                           y = rsqured)) +
  geom_point(size = 2.5, color = "grey30") +
  geom_line(color = "#20A387FF", size = 1.2) +
  scale_x_continuous(breaks = c(seq(10, 140, 10))) +
  geom_vline(xintercept = 63, lty = 3) +
  labs(x = "")
num_tuning_p2 = ggplot(data = test_df,
                       aes(x = num,
                           y = rmse)) +
  geom_point(size = 2.5, color = "grey30") +
  geom_line(color = "#20A387FF", size = 1.2) +
  scale_x_continuous(breaks = c(seq(10, 140, 10))) +
  geom_vline(xintercept = 63, lty = 3)

ggarrange(num_tuning_p1, num_tuning_p2, nrow = 2, ncol = 1, align = "hv")

## 9b. building a model with 63 features (fig.7a)
discovery_allfea_dat = time_fi_all %>%
  mutate(age = `Age.at.assesment..days.`,
         sex_bin = ifelse(Sex == "female", 1, 0)) %>%
  dplyr::select(fi, unlist(colnames(dis_quant_dat_excl)), sex_bin, age) %>%
  `colnames<-`(make.names(colnames(.)))

discovery_63fea_dat = time_fi_all %>%
  mutate(age = `Age.at.assesment..days.`,
         sex_bin = ifelse(Sex == "female", 1, 0)) %>%
  dplyr::select(fi, unlist(total_fea[1:63]), sex_bin, age) %>%
  `colnames<-`(make.names(colnames(.)))

set.seed(2024)
rf_all_model = rf_model(discovery_allfea_dat)
rf_63_model = rf_model(discovery_63fea_dat)
sexage_model = lm(fi ~ age + as.factor(sex_bin), data = discovery_63fea_dat)

sexage_pred = predict(sexage_model, discovery_63fea_dat[, c("age", "sex_bin")])
rf_all_pred = predict(rf_all_model, makeX(discovery_allfea_dat[, -1]))
rf_63_pred = predict(rf_63_model, makeX(discovery_63fea_dat[, -1]))

discovery_dat_pred = discovery_63fea_dat %>%
  mutate(agesex_pred = unname(sexage_pred),
         all_pred = unname(rf_all_pred),
         fea63_pred = unname(rf_63_pred)) %>%
  dplyr::select(fi, agesex_pred, all_pred, fea63_pred, sex_bin, age) %>%
  melt(., id = c("fi", "sex_bin", "age")) %>%
  as.data.frame()

plot_fig7a = ggplot(data = discovery_dat_pred,
                    aes(x = fi,
                        y = value,
                        color = as.factor(variable))) +
  geom_point(size = 0.5) + 
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = 1.5) +
  stat_cor(aes(label = paste(..rr.label..))) +
  scale_color_manual(values = viridis(100)[c(10, 50, 80)]) +
  coord_fixed(ratio = 1.5) +
  geom_abline(slope = 1, lty = 3) +
  theme(legend.position = "bottom")

## 9c. model performance in nmn (fig.7b)
nmn_allfea_dat = nmn_dat %>% 
  mutate(age = `Age.at.assesment..days.`,
         sex_bin = ifelse(Sex == "female", 1, 0)) %>%
  dplyr::select(fi, unlist(colnames(dis_quant_dat_excl)), sex_bin, age) %>%
  `colnames<-`(make.names(colnames(.)))
nmn_63fea_dat = nmn_dat %>%
  mutate(age = `Age.at.assesment..days.`,
         sex_bin = ifelse(Sex == "female", 1, 0)) %>%
  dplyr::select(fi, unlist(total_fea[1:63]), sex_bin, age) %>%
  `colnames<-`(make.names(colnames(.)))

nmn_sexage_model = lm(fi ~ age + as.factor(sex_bin), data = nmn_63fea_dat)
nmn_sexage_pred = predict(nmn_sexage_model, nmn_63fea_dat[, c("age", "sex_bin")])
nmn_rf_all_pred = predict(rf_all_model, makeX(nmn_allfea_dat[, -1]))
nmn_rf_63_pred = predict(rf_63_model, makeX(nmn_63fea_dat[, -1]))

nmn_dat_pred = nmn_63fea_dat %>%
  mutate(agesex_pred = unname(nmn_sexage_pred),
         all_pred = unname(nmn_rf_all_pred),
         fea63_pred = unname(nmn_rf_63_pred)) %>%
  dplyr::select(fi, agesex_pred, all_pred, fea63_pred, sex_bin, age) %>%
  melt(., id = c("fi", "sex_bin", "age")) %>%
  as.data.frame()

plot_fig7b = ggplot(data = nmn_dat_pred,
                    aes(x = fi,
                        y = value,
                        color = as.factor(variable))) +
  geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = 1.5) +
  stat_cor(aes(label = paste(..rr.label..))) +
  scale_color_manual(values = viridis(100)[c(10, 50, 80)]) +
  coord_fixed(ratio = 1.5) +
  geom_abline(slope = 1, lty = 3) +
  theme(legend.position = "bottom")

## 9d. model performance in nmn for females and males (fig.7c)
nmn_dat_pred_agetrain = nmn_63fea_dat %>%
  mutate(agesex_pred = unname(nmn_sexage_pred),
         fea63_pred = unname(nmn_rf_63_pred)) %>%
  dplyr::select(fi, agesex_pred, fea63_pred, sex_bin) %>%
  melt(., id = c("fi", "sex_bin")) %>%
  as.data.frame()

plot_fig7c = ggplot(data = nmn_dat_pred_agetrain,
                    aes(x = fi,
                        y = value,
                        color = paste(as.factor(variable), sex_bin))) +
  geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = 1.5) +
  stat_cor(aes(label = paste(..rr.label..))) +
  scale_color_manual(values = viridis(100)[c(20, 50, 70, 90)]) +
  coord_fixed(ratio = 1.5) +
  geom_abline(slope = 1, lty = 3) +
  theme(legend.position = "bottom")

ggarrange(plot_fig7a, plot_fig7b, plot_fig7c,
          nrow = 1, ncol = 3, align = "hv")

sessionInfo()

