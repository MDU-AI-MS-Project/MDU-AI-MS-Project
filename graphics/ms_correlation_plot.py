import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Filepath for the correlation matrix CSV
root = Path("/mnt/c/MDU/GNPSnew")
csv_file = root / "predictions_correlation_spearman.csv"  # Replace with your CSV file path
heatmap_svg_file = root / "predictions_correlation_spearman_heatmap_diagram.svg"  # Output SVG file path

# Read the correlation matrix from the CSV file
correlation_matrix = pd.read_csv(csv_file, index_col=0)

# Generate the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", fmt=".2f", square=True)

# Tilt the x-axis labels
plt.xticks(rotation=45, ha="right")  # Rotate labels 45 degrees and align them to the right

# Save the heatmap as an SVG file
plt.savefig(heatmap_svg_file, format="svg", bbox_inches="tight")

print(f"Heatmap diagram saved as SVG to {heatmap_svg_file}")