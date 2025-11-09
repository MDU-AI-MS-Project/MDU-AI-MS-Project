from pathlib import Path
import re

# Filepaths for the SVG files
root = Path("/mnt/c/MDU/GNPSnew")
wide_bottom_file3 = root / "ms_shap_tree_fragment.svg"  # Replace with your wide SVG file path
square_rightup_file3 = root / "ms_shap_beeswarm_plot_with_impact_in_labels.svg"  # Replace with your first square SVG file path
square_leftup_file3 = root / "predictions_correlation_spearman_heatmap_diagram_R.svg"  # Replace with your second square SVG file path
output_svg_file3 = root / "merged_figure3.svg"  # Output merged SVG file path

square_left_file2 = root / "model_performance_charts" / "model_performance_chart_alldesc_allpred.svg"
square_right_file2 = root / "model_performance_charts" / "descriptor_impact_chart.svg"
output_svg_file2 = root / "merged_figure2.svg"  # Output merged SVG file path


# Function to preprocess SVG content and remove XML declarations, <!DOCTYPE>, and <svg> tags
def preprocess_svg(svg_content, font_scale=1.0):
    # Remove XML declaration if present
    if svg_content.startswith('<?xml'):
        svg_content = svg_content.split('?>', 1)[1]
    # Remove <!DOCTYPE> declaration using regex to handle arbitrary whitespace
    svg_content = re.sub(r'<!DOCTYPE\s+svg\s+PUBLIC\s+".*?"\s+".*?">', '', svg_content, flags=re.DOTALL)
    # Remove <svg> tags
    svg_content = svg_content.replace('<svg', '<g').replace('</svg>', '</g>')
    # Remove any invalid or conflicting namespaces
    svg_content = svg_content.replace('xmlns="http://www.w3.org/2000/svg"', '')
    svg_content = svg_content.replace('xmlns:xlink="http://www.w3.org/1999/xlink"', '')
    # Scale fonts if needed
    svg_content = re.sub(r'font-size="(\d+\.?\d*)"', lambda m: f'font-size="{float(m.group(1)) * font_scale}"', svg_content)
    return svg_content.strip()

# Read and preprocess the content of the SVG files
with open(wide_bottom_file3, "r") as f:
    wide_bottom_content3 = preprocess_svg(f.read())

with open(square_rightup_file3, "r") as f:
    square_rightup_content3 = preprocess_svg(f.read(), font_scale=1.2)

with open(square_leftup_file3, "r") as f:
    square_leftup_content3 = preprocess_svg(f.read(), font_scale=1.5)

# Create the merged SVG content
merged_svg_content3 = f"""
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="2000px" height="2200px">
    <!-- Second square image (left-top) -->
    <g transform="translate(0, 0) scale(1.2)">
        {square_leftup_content3}
    </g>
    <!-- First square image (right-top) -->
    <g transform="translate(970, 40) scale(1.6)">
        {square_rightup_content3}
    </g>
    <!-- Wide image at bottom -->
    <g transform="translate(100, 900)">
        {wide_bottom_content3}
    </g>
</svg>
"""

# Save the merged SVG file
with open(output_svg_file3, "w") as f:
    f.write(merged_svg_content3)
print(f"Merged SVG figure saved to {output_svg_file3}")

with open(square_right_file2, "r") as f:
    square_right_content2 = preprocess_svg(f.read(), font_scale=1.2)

with open(square_left_file2, "r") as f:
    square_left_content2 = preprocess_svg(f.read(), font_scale=1.2)

# Create the merged SVG content
merged_svg_content2 = f"""\
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="2000px" height="1000px">
    <!-- image (left) -->
    <g transform="translate(30, 100) scale(1.05)">
        {square_left_content2}
    </g>
    <!-- image (right) -->
    <g transform="translate(1000, 100) scale(1.05)">
        {square_right_content2}
    </g>
</svg>
"""
# Save the merged SVG file
with open(output_svg_file2, "w") as f:
    f.write(merged_svg_content2)
print(f"Merged SVG figure saved to {output_svg_file2}")

