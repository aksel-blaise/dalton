## Cluster analysis

```{r cluster, out.width = "100%", dpi = 300, echo=TRUE, warning=FALSE}
# attributes
pcdata <- as_df(pca.outlines)
row.names(pcdata) <- pcdata$File_Name

# cluster analysis
cluster <- CLUST(pca.outlines, 
                 dist_method = "euclidean", 
                 hclust_method = "complete",
                 retain = 1:10,
                 tip_labels = "File_Name")
clusterdata <- select(pcdata, heart.out, bev) %>% 
  droplevels()

# plot 1
ggtree(cluster, layout = "circular") %>% 
  clusterdata +
  geom_tippoint(aes(colour = heart.reg)) +
  geom_tiplab2(aes(label = heart.reg, colour = heart.out, align = TRUE),
               offset = 0.005, size = 2.5) +
  labs(colour = "Heartland") +
  scale_colour_manual(values = wes_palette("Moonrise2")) +
  guides(col = guide_legend(ncol = 3, byrow = TRUE)) +
  theme(legend.position = "bottom")
```