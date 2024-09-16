# Install gplots package if not installed
if(!require(gplots)) install.packages("gplots")
library(gplots)
# Remove non-numeric columns and scale the data
data_for_clustering <- glioblastoma_data[, sapply(glioblastoma_data, is.numeric)]
scaled_data <- t(scale(t(data_for_clustering)))  # Mean-center and scale by rows (genes)
# Heatmap using a diverging color palette (blue-white-red)
heatmap.2(scaled_data, 
          Rowv = TRUE, 
          Colv = TRUE, 
          dendrogram = "both", 
          #col = bluered(100),  # Diverging color palette
          col = colorRampPalette(c("blue", "white", "red"))(100),  # Blue to white to red
          trace = "none",
          margins = c(10, 10), 
          key = TRUE, 
          keysize = 1.5,
          density.info = "none",
          cexCol = 0.7,
          main = "Heatmap with Diverging Palette")


# Heatmap using a sequential color palette (light blue to dark blue)
x11()
heatmap.2(scaled_data, 
          Rowv = TRUE, 
          Colv = TRUE, 
          dendrogram = "both", 
          col = colorRampPalette(c("lightblue", "blue", "darkblue"))(100),  # Sequential palette
          trace = "none",
          margins = c(10, 10), 
          key = TRUE, 
          keysize = 1.5,
          density.info = "none",
          cexCol = 0.7,
          main = "Heatmap with Sequential Palette")

# Cluster only rows (genes)
x11()
heatmap.2(scaled_data, 
          Rowv = TRUE,         # Cluster rows (genes)
          Colv = FALSE,        # Do not cluster columns (samples)
          dendrogram = "row",  # Show dendrogram for rows only
          col = bluered(100),  # Diverging color palette
          trace = "none",
          margins = c(10, 10), 
          key = TRUE, 
          cexCol = 0.7,
          density.info = "none",
          main = "Clustered Genes (Rows) Only")


# Cluster only columns (samples)
x11()
ggcorrplot(scaled_data, 
          Rowv = FALSE,        # Do not cluster rows (genes)
          Colv = TRUE,         # Cluster columns (samples)
          dendrogram = "column",# Show dendrogram for columns only
          col = bluered(100),  # Diverging color palette
          trace = "none",
          margins = c(10, 10), 
          key = TRUE, 
          cexCol = 0.7,
          density.info = "none",
          main = "Clustered Samples (Columns) Only")


install.packages("ggcorrplot")
library(ggcorrplot)

# Cluster both rows (genes) and columns (samples)
x11()
heatmap.2(scaled_data, 
          Rowv = TRUE,         # Cluster rows (genes)
          Colv = TRUE,         # Cluster columns (samples)
          dendrogram = "both", # Show dendrogram for both rows and columns
          col = bluered(100),  # Diverging color palette
          trace = "none",
          margins = c(10, 10), 
          key = TRUE, 
          cexCol = 0.7,
          density.info = "none",
          main = "Clustered Genes and Samples Together")


