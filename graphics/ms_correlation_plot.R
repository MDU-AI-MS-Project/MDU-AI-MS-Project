# Load required libraries
library(corrplot)
library(readr)

# Define file paths
root <- "/mnt/c/MDU/GNPSnew"
csv_file <- file.path(root, "predictions_correlation_spearman.csv")
heatmap_svg_file <- file.path(root, "predictions_correlation_spearman_heatmap_diagram_R.svg")

# Read the correlation matrix from the CSV file
correlation_matrix <- read_csv(csv_file, show_col_types = FALSE)

# Convert to matrix format (assuming first column contains row names)
row_names <- correlation_matrix[[1]]
correlation_matrix <- correlation_matrix[, -1]
correlation_matrix <- as.matrix(correlation_matrix)
rownames(correlation_matrix) <- row_names

# Check for and handle missing/infinite values
print("Checking for missing/infinite values:")
print(paste("Missing values:", sum(is.na(correlation_matrix))))
print(paste("Infinite values:", sum(is.infinite(correlation_matrix))))

# Replace missing/infinite values with 0 (or use na.omit if preferred)
correlation_matrix[is.na(correlation_matrix)] <- 0
correlation_matrix[is.infinite(correlation_matrix)] <- 0

# Print actual data range
print(paste("Data range: min =", round(min(correlation_matrix), 3), 
           "max =", round(max(correlation_matrix), 3)))

# Define color palette from white to dark blue
color_palette1 <- colorRampPalette(c("#f7f7f7", "#d1e5f0", "#92c5de", 
                                   "#4393c3", "#2166ac", "#053061"))(10)

# Define color palette RBG
GradBlueCompressed <- function(n = 200, exponent = 0.37) {
  # n: number of colors in the palette
  # exponent: compression exponent (<1 = more colors in low range)
  
  # full linear sequence
  x <- seq(0, 1, length.out = n)
  
  # compressed sequence
  x_compressed <- x^exponent
  
  # base gradient from lightblue to darkblue
  base_cols <- colorRampPalette(c("white", "lightblue", "blue", "purple"))(n)
  
  # reorder according to compressed mapping
  compressed_cols <- base_cols[ceiling(x_compressed * n)]
  print(head(compressed_cols))
  print(tail(compressed_cols))
  return(compressed_cols)
}

color_palette2 <- colorRampPalette(c("#670000", "#B20000", "#D60000", "#F40000", "#FD0000", "#FFFFFF"))(600)
my.col <- colorRampPalette(c("#ADD8E6","#A0CFE0","#93C6DB","#86BDD6","#79B4D0","#6CAACB","#5FA1C6","#5298C0","#458FBB","#3886B6","#2B7CB0","#1E73AB","#116AA6","#0461A0","#00599B","#005094","#00488D","#004086","#00387F","#003078","#002871","#00206A","#001863","#00105C","#000855","#00004E","#000047","#000040","#000039","#00008B"
))  # your colors 
set.seed(0)
# Open SVG device
svg(heatmap_svg_file, width = 12, height = 10)

# First plot: circles in lower triangle only (no color legend)
corrplot(correlation_matrix/2, 
         method = "circle", 
         type = "lower",
         is.corr = FALSE,
         col = GradBlueCompressed(), # Use built-in blue-white-red palette
         col.lim = c(0.0, 1.4),     # Adjust color limits to match data range
         #outline = TRUE,
         tl.pos = "n",              # No text labels on this plot
         cl.pos = "n",              # No color legend
         #itle = "Correlation Between Individual Predictors and Reference Data",                # No title here
         mar = c(2, 1, 3, 2))

# Second plot: numbers in upper triangle only (add to existing plot)
corrplot(correlation_matrix, 
         method = "number", 
         type = "upper", 
         add = TRUE,                # Add to existing plot
         tl.pos = "d",              # Text labels on diagonal
         tl.col = "black",          # Black diagonal labels
         tl.srt = 45,               # Rotate diagonal labels 45 degrees
         tl.cex = 1.2,              # Larger diagonal text
         number.cex = 1.2,          # Larger numbers
         number.font = 2,           # Bold numbers
         number.digits = 3,         # 3 decimal places
         col = '#096bde',           # Black numbers
         cl.pos = "n")              # No color legend

# Close the SVG device
dev.off()

cat("Heatmap diagram saved as SVG to", heatmap_svg_file, "\n")