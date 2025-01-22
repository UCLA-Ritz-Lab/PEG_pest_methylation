test <- gsea_champ_copper_total$DMP %>% 
  mutate(gene_ratio = DE/N)

test_sort <- test %>% 
  arrange(P.DE) %>% 
  mutate(term = fct_reorder(TERM, P.DE))


test_sort %>% 
  mutate(term = fct_rev(term)) %>%
  group_by(ONTOLOGY) %>% 
  ggplot(aes(x = -log10(P.DE), y = term, color = gene_ratio)) + 
  scale_color_gradient(low = "blue", high = "red") +
  geom_point() +
  facet_grid(rows = vars(ONTOLOGY), scales = "free_y") +
  labs(y = "GO terms", 
       x = "-log<sub>10</sub>(p)",
       color = "Gene ratio") + 
  theme_bw()+
  theme( 
    # legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_markdown(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text.y = element_text(size = 8), 
    axis.text.x = element_text(size = 15, vjust = 0.5)
  )
