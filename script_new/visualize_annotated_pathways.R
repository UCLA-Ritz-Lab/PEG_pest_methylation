## ---------------------------
##
## Script name: visualize_annotated_pathways.R
## Purpose of script: To visualize annotated pathways from the ORA and pathway analysis
##
## Author: Yufan Gong
##
## Date Created: 2025-07-29
##
## Date Modified: 2025-07-29
##
## Copyright (c) Yufan Gong, 2025
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes: 

## load libraries

library(tidyverse)

## load results

list.dirs(here::here(),recursive = FALSE) %>%
  list.files("\\.RData$", full.names = TRUE, recursive = T) %>%
  grep("dmp_count|gsea_champ|gst_kegg|panther_efgr",.,
       value=TRUE, ignore.case = TRUE) %>%
  discard(~str_detect(.x,"old")) %>%
  map(.,load,.GlobalEnv) %>% 
  invisible()

## map cpgs to genes ---------
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

cpg_gene_map <- ann450k %>%
  as.data.frame() %>%
  dplyr::select(CpG = Name, Gene = UCSC_RefGene_Name) %>%
  separate_rows(Gene, sep = ";") %>%  # Some CpGs map to multiple genes
  filter(Gene != "") %>% 
  rename_all(tolower) %>% 
  filter(cpg %in% rownames(dmp_count_total$total))

sig_cpg <- dmp_count_total$total %>% 
  filter(-log10(P.Value) > 6) %>% 
  arrange(P.Value) %>% 
  rownames_to_column("cpg")

# Filter for CpGs in your list
sig_cpg_gene_map <- cpg_gene_map %>%
  filter(cpg %in% sig_cpg$cpg)


## Get the list of genes associated with GO terms -----

# gsea_dmr <- gsea_champ_copper_total$DMR

library(org.Hs.eg.db)

go_list <- c("GO:0007156", "GO:0005654", "GO:0051641", "GO:0016740") %>% 
  set_names(c("Cell adhesion", "nucleoplasm", 
              "cellular_localization", "transferase_activity"))

# Get Entrez IDs annotated to the GO terms

go_list %>% 
  map(function(go_term) {
    AnnotationDbi::select(org.Hs.eg.db, 
                          keys = go_term, 
                          keytype = "GO", 
                          columns = c("ENTREZID", "SYMBOL"))
  }) -> entrez_ids_list

# map genes to the current cpg_gene_map

entrez_ids_list %>% 
  map(function(df) {
    cpg_gene_map %>%
      filter(gene %in% df$SYMBOL)
  }) -> gene_filtered_list


# sig_genes_in_go_term <- intersect(sig_cpg_gene_map$gene, entrez_ids$SYMBOL)


## Get the list of genes associated with Phagocytosis (hsa04666) -----
# library(missMethyl)
# 
# topCpGs_total <- dmp_count_total$total %>% 
#   filter(-log10(P.Value) > 6) %>% 
#   rownames_to_column("cpg") %>%
#   arrange(P.Value)
# 
# 
# sigCpGs_total_all <- topCpGs_total %>% 
#   pull(cpg)
# 
# gst_kegg_total <- gometh(sig.cpg = sigCpGs_total_all, 
#                          all.cpg = rownames(PEG_NOOB_nors_win_filter_total), 
#                          collection = "KEGG", 
#                          plot.bias = FALSE)

# save(gst_kegg_total, file = "gst_kegg_total.RData")

library(clusterProfiler)
# library(gage)
# human_kegg_sets <- kegg.gsets("hsa")

# Get gene list for hsa04666
genes_in_kegg <- enrichKEGG(gene = NULL, organism = 'hsa')  # Dummy run to activate KEGG db
kegg_pathway_genes_Phagocytosis <- as.character(KEGGREST::keggGet("hsa04666")[[1]]$GENE)

gene_symbols <- kegg_pathway_genes_Phagocytosis[seq(2, length(kegg_pathway_genes_Phagocytosis), by = 2)]  # Even indices
entrez_ids <- kegg_pathway_genes_Phagocytosis[seq(1, length(kegg_pathway_genes_Phagocytosis), by = 2)]    # Odd indices

kegg_gene_df_Phagocytosis <- data.frame(ENTREZID = entrez_ids, 
                      SYMBOL = sapply(strsplit(gene_symbols, ";"), `[`, 1))

cpgs_in_kegg_Phagocytosis <- cpg_gene_map %>%
  filter(gene %in% kegg_gene_df_Phagocytosis$SYMBOL)

## Get the list of genes associated with EGFR (P00018) -----

library(PANTHER.db)
pthOrganisms(PANTHER.db) <- "HUMAN"
PANTHER.db
columns(PANTHER.db)
keytypes(PANTHER.db)

go_ids <- keys(PANTHER.db, keytype="PATHWAY_ID")
cols <- "ENTREZ"
panther <- mapIds(PANTHER.db, keys=go_ids, column=cols, 
                  as.numeric(cols), keytype="PATHWAY_ID", multiVals="list")

panther_gene_df_egfr <- data.frame(
  ENTREZID = panther$P00018, 
  SYMBOL = mapIds(org.Hs.eg.db, keys = panther$P00018, 
                  column = "SYMBOL", keytype = "ENTREZID")
)

cpgs_in_panther_egfr <- cpg_gene_map %>%
  filter(gene %in% panther_gene_df_egfr$SYMBOL)

cpg_gene_list_final <- c(
  gene_filtered_list,
  list(
    `Phagocytosis` = cpgs_in_kegg_Phagocytosis,
    `EGFR signaling` = cpgs_in_panther_egfr
  )
)

cpg_gene_pathway <- list(cpg_gene_list_final,
     names(cpg_gene_list_final)) %>% 
  pmap(function(data, pathway){
    data %>% 
      mutate(pathway = pathway)
  }) %>% 
  bind_rows() %>% 
  distinct()


# Visualization -----------------------------------------------------------


## 1. Manhattan plot ----------
all_cpgs <- dmp_count_total$total%>% 
  dplyr::rename(chr = CHR,
                p.value = P.Value) %>% 
  mutate(cpg = rownames(.))

ann_filter <- ann450k %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% rownames(all_cpgs)) %>% 
  dplyr::rename(cpg = Name) %>% 
  dplyr::select(cpg, pos)

all_cpgs_clean <- all_cpgs %>% 
  left_join(ann_filter, by = "cpg") %>% 
  dplyr::rename(position = pos)


chromosomes <- c(1:22)
chromosomes <- intersect(chromosomes, all_cpgs_clean$chr)
chromosome.lengths <- as.numeric(sapply(chromosomes, function(chromosome)
  max(all_cpgs_clean$position[which(all_cpgs_clean$chr == chromosome)])))
chromosome.starts <- c(1,cumsum(chromosome.lengths)+1)

names(chromosome.starts) <- c(chromosomes, "NA")

chromosome.starts_df <- tibble(
  chr = chromosomes,
  chromosome.starts_new = chromosome.starts[1:22]
)



all_cpgs_final <- all_cpgs_clean %>% 
  group_by(chr) %>% 
  mutate(chromosome.lengths = max(position)) %>% 
  ungroup() %>% 
  left_join(chromosome.starts_df, by = "chr") %>% 
  dplyr::rename(chromosome.starts_clean = chromosome.starts_new) %>% 
  mutate(global = position + chromosome.starts_clean - 1) %>% 
  left_join(cpg_gene_pathway %>% 
              rename(gene_filter = gene), by = "cpg") %>% 
  arrange(p.value) %>% 
  mutate(sig = if_else(-log10(p.value) > 6, 1, 0),
         sig_path = if_else(-log10(p.value) > 2.5, 1, 0),
         label = if_else(sig == 1 & pathway %in% c("Cell adhesion", 
                                                   "Phagocytosis", 
                                                   "EGFR signaling"), cpg, ""))

all_cpgs_clean_new <- all_cpgs_final %>% 
  group_by(chr) %>% 
  dplyr::summarize(center = mean(global)) %>% 
  ungroup() %>% 
  arrange(center) %>% 
  mutate(chr = factor(chr, levels = c(1:22)))


ylim <- all_cpgs_final %>% 
  filter(p.value == min(p.value)) %>% 
  mutate(ylim = abs(floor(log10(p.value))) + 1) %>% 
  pull(ylim)

library(ggrepel)
library(ggtext)
library(ggnewscale)
library(Polychrome)
# devtools::install_github("hrbrmstr/hrbrthemes")
library(hrbrthemes)

set.seed(19990121)
# mypal <- c("#2774AE", "#FFB81C", "#2E8B57") #"#d43325"red 
mypal_path <- c("#2b6a99", "#f16c23", "#1b7c3d")
mypal_path_unsig <- c("#6ea1c3", "#f49a60", "#60b283")
mypal_path_unsig_light <- c("#a8c6df", "#f8b98a", "#9cd3ae")


names(mypal_path) <- names(cpg_gene_list_final) %>% 
  discard(~str_detect(.x, "nucleoplasm|cellular|transferase"))

names(mypal_path_unsig) <- names(cpg_gene_list_final) %>% 
  discard(~str_detect(.x, "nucleoplasm|cellular|transferase"))
  

# mypal_grey <- rep(c("#d3d3d3", "#696969"), 11)
mypal_grey <- rep(c("#b7b5b6", "#b7b5b6"), 11)
names(mypal_grey) <- c(1:22)


## Plotting the results

pathway_df_go <- list(
  gsea_champ_copper_total$DMP,
  gsea_champ_copper_total$DMR
) %>% 
  bind_rows() %>% 
  filter(row.names(.) %in% go_list) %>% 
  dplyr::select(TERM, DE, N, P.DE)

pathway_df_kegg <- gst_kegg_total %>% 
  rename(TERM = Description) %>% 
  filter(row.names(.) == "hsa04666") %>% 
  dplyr::select(TERM, DE, N, P.DE)

pathway_df_panther <- methylglm_path_total_op %>% 
  filter(ID == "P00018") %>%
  rename(TERM = pathway) %>% 
  dplyr::select(TERM, DE, N, P.DE)

pathway_df_final <- list(
  pathway_df_go,
  pathway_df_kegg,
  pathway_df_panther
) %>% 
  bind_rows() %>% 
  mutate(gene_ratio = DE/N,
         category = case_when(
           str_detect(TERM, "EGF") ~ "Panther pathway",
           str_detect(TERM, "phagocytosis") ~ "KEGG pathway",
           TRUE ~ "GO terms"
         ))

pathway_df_final_sort <- pathway_df_final %>% 
  arrange(P.DE) %>% 
  mutate(term = fct_reorder(TERM, P.DE))

all_cpgs_final_clean <- all_cpgs_final %>%
  mutate(
    term = case_when(
      pathway == "Cell adhesion" ~ "homophilic cell adhesion via plasma membrane adhesion molecules",
      pathway == "Phagocytosis" ~ "Fc gamma R-mediated phagocytosis",
      pathway == "EGFR signaling" ~ "EGF receptor signaling pathway"
    )) %>% 
  left_join(pathway_df_final_sort, by = "term") 

# %>% 
#   mutate(
#     pvalue_round = signif((P.DE), digits = 2),
#     pathway_new = str_c(pathway, " (", "Gene Ratio=", DE, "/", N, ", p=", pvalue_round, ")"))

pathway_labels <- setNames(
  Filter(Negate(is.na), all_cpgs_final_clean$pathway),  # values shown in the legend
  Filter(Negate(is.na), names(cpg_gene_list_final)) # levels in the color aes
) %>% unique()


all_cpgs_final_clean %>% 
  mutate(label = if_else(label == "cg02313172" 
                         & pathway == "Phagocytosis", "", label)) %>% 
  ggplot(aes(x = global, y = -log10(p.value), size = -log10(p.value))) +
  geom_point(aes(color = as_factor(chr)), alpha = 0.3, 
             size = 2.5, data = . %>% filter(sig == 0)) +
  scale_color_manual(values = mypal_grey, guide = "none") +
  new_scale_color() +
  geom_point(aes(color = as_factor(sig)), alpha = 0.4, 
             size = 2.5, data = . %>% filter(sig == 1)) +
  # scale_color_manual(values = c("1" = "#D1495B"), guide = "none") +
  scale_color_manual(values = c("1" = "#b7b5b6"), guide = "none") +
  new_scale_color() +
  geom_point(aes(color = as_factor(pathway)), data = . %>% 
               filter(pathway %in% c("Cell adhesion", 
                                     "Phagocytosis", 
                                     "EGFR signaling") 
                      & sig_path == 1), size = 6) +
  scale_color_manual(values = mypal_path, 
                     # labels = pathway_labels, 
                     name = "Pathway") +
  new_scale_color() +
  geom_point(aes(color = as_factor(pathway)), data = . %>% 
               filter(pathway %in% c("Cell adhesion", 
                                     "Phagocytosis", 
                                     "EGFR signaling") 
                      & sig_path == 0), size = 6) +
  scale_color_manual(values = mypal_path_unsig, 
                     # labels = pathway_labels, 
                     name = "Pathway (unassociated CpG)",
                     guide = "none") +
  geom_hline(yintercept = -log10(10e-7), color = "red3", linetype = "solid") +
  geom_hline(yintercept = 2.5, color = "#F8766D", linetype = "dashed") +
  geom_text_repel(aes(label = label),
                  size = 10,
                  box.padding = 1,
                  point.padding = 2.5, 
                  force = 20,
                  # force_pull = 2,
                  nudge_x = 1,
                  nudge_y = 1,
                  max.overlaps = Inf) +
  # geom_label_repel(aes(label = label),
  #                 size = 5,
  #                 box.padding = 1,
  #                 nudge_x = 0.25,
  #                 nudge_y = 0.25,
  #                 max.overlaps = Inf) +
  scale_x_continuous(label = all_cpgs_clean_new$chr, breaks = all_cpgs_clean_new$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) + 
  #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis$chr)))) +
  # scale_size_continuous(range = c(3,3), guide = "none") +
  labs(x = "Chromosome", 
       y = "-log<sub>10</sub>(p)") + 
  theme_classic() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.125),
    legend.direction = "vertical", 
    legend.box = "horizontal",
    # legend.box.spacing = unit(0, "mm"),
    legend.key.height = unit(1.1, "cm"),
    # legend.justification = c("right", "top"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.text = element_text(size = 25),
    legend.title = element_text(face = "bold", size = 25),
    axis.title.y = element_markdown(face = "bold", size = 25),
    axis.title.x = element_text(face = "bold", size = 25),
    axis.text.y = element_text(size = 20), 
    axis.text.x = element_text(size = 17, vjust = 0.5)
  )


## 2. lollipop plot ----------



library(ggtext)

color_plate <- colorRampPalette(c("#008bd0", "#eeeeee", "#ffa61d"))(100)

p1 <- pathway_df_final_sort %>% 
  mutate(term = fct_rev(term)) %>%
  # group_by(ONTOLOGY) %>% 
  # slice_head(n = 5) %>% 
  ggplot(aes(x = term, y = -log10(P.DE))) + 
  coord_flip() +
  geom_col(aes(fill = category), 
           width = 0.05) +
  geom_point(aes(color = category, size = gene_ratio), shape = 16) +
  scale_size_continuous(range = c(5,8), name = "Gene Ratio") +
  scale_fill_manual(values = c("GO terms" = "#2774AE", 
                               "KEGG pathway" = "#FFB81C", 
                               "Panther pathway" = "#2E8B57"),
                    name = "Pathway") +
  scale_color_manual(values = c("GO terms" = "#2774AE", 
                                 "KEGG pathway" = "#FFB81C", 
                                 "Panther pathway" = "#2E8B57"),
                     guide = "none") +
  # scale_color_gradientn(colors = color_plate,
  #                       values = scales::rescale(seq(0, 6 , length.out = 100)),
  #                       name = "Gene ratio") +
  # scale_color_gradient(low = "blue", high = "red") +
  
  # facet_wrap(~ONTOLOGY, scales = "free_y") +
  labs(x = NULL, 
       y = bquote("-log<sub>10</sub>(p)")) + 
  theme_ipsum() +
  theme( 
    # legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_markdown(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text.y = element_blank(), 
    axis.text.x = element_text(size = 20, vjust = 0.5),
    axis.ticks = element_blank()
  )

p2 <- pathway_df_final_sort %>% 
  mutate(term = fct_rev(term)) %>%
  # group_by(ONTOLOGY) %>% 
  # slice_head(n = 5) %>% 
  ggplot(aes(x = term, y = 0.05)) + 
  geom_col(
    aes(fill = category),
    width = 1,
    position = position_dodge(width = 0.8)
  )+
  # geom_text(
  #   aes(label = c("GO terms", 
  #                 "KEGG pathway",
  #                 "Panther pathway"), y = 0.025),
  #   color = "white",
  #   angle = -90,
  #   hjust = 0.5,
  #   size = 5
  # ) +
  coord_flip() +
  scale_fill_manual(
    values = c("GO terms" = "#2774AE", 
               "KEGG pathway" = "#FFB81C", 
               "Panther pathway" = "#2E8B57")
  ) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

library(patchwork)

p2 + 
  p1 + 
  plot_layout(widths = c(1, 8))
