test <- dmp_count_total$total%>% 
  dplyr::rename(chr = CHR,
                p.value = P.Value) %>% 
  mutate(cpg = rownames(.))


library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

annotation_filter <- annotation %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% rownames(test)) %>% 
  dplyr::rename(cpg = Name) %>% 
  dplyr::select(cpg, pos)

test_clean <- test %>% 
  left_join(annotation_filter, by = "cpg") %>% 
  dplyr::rename(position = pos)

chromosomes <- c(1:22)
chromosomes <- intersect(chromosomes, test_clean$chr)
chromosome.lengths <- as.numeric(sapply(chromosomes, function(chromosome)
  max(test_clean$position[which(test_clean$chr == chromosome)])))
chromosome.starts <- c(1,cumsum(chromosome.lengths)+1)

names(chromosome.starts) <- c(chromosomes, "NA")

chromosome.starts_df <- tibble(
  chr = chromosomes,
  chromosome.starts_new = chromosome.starts[1:22]
)

test_clean <- test_clean %>% 
  group_by(chr) %>% 
  mutate(chromosome.lengths = max(position)) %>% 
  ungroup() %>% 
  left_join(chromosome.starts_df, by = "chr") %>% 
  dplyr::rename(chromosome.starts_clean = chromosome.starts_new) %>% 
  mutate(global = position + chromosome.starts_clean - 1)

test_clean_new <- test_clean %>% 
  group_by(chr) %>% 
  dplyr::summarize(center = mean(global)) %>% 
  ungroup() %>% 
  arrange(center) %>% 
  mutate(chr = factor(chr, levels = c(1:22)))


ylim <- test_clean %>% 
  filter(p.value == min(p.value)) %>% 
  mutate(ylim = abs(floor(log10(p.value))) + 1) %>% 
  pull(ylim)

library(ggrepel)
library(ggtext)
library(ggnewscale)
library(Polychrome)

set.seed(19990121)
mypal <- Polychrome::createPalette(22,  c("#999999", "#E69F00", "#56B4E9"))
names(mypal) <- c(1:22)

test_clean_final <- test_clean %>% 
  arrange(p.value) %>% 
  mutate(sig = if_else(-log10(p.value) > 6, 1, 0),
         
         label = if_else(row_number() <= 10, cpg, ""))


test_clean_final %>%
  ggplot(aes(x = global, y = -log10(p.value), size = -log10(p.value))) +
  geom_point(aes(color = as_factor(chr)), alpha = 0.75, data = . %>% filter(sig == 0)) +
  scale_color_manual(values = mypal) +
  new_scale_color() +
  geom_point(aes(color = as_factor(sig)), data = . %>% filter(sig == 1)) +
  scale_color_manual(values = c("1" = "red3")) +
  geom_hline(yintercept = -log10(10e-7), color = "red", linetype = "dashed") + 
  geom_label_repel(aes(label = label), 
                   size = 4,
                   box.padding = 1,
                   nudge_x = 0.25,
                   nudge_y = 0.25) +
  scale_x_continuous(label = test_clean_new$chr, breaks = test_clean_new$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome", 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(face = "bold", size = 25),
    axis.title.x = element_text(face = "bold", size = 25),
    axis.text.y = element_text(size = 20), 
    axis.text.x = element_text(size = 20, vjust = 0.5)
  )
