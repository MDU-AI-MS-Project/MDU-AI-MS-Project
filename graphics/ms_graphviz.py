

from pathlib import Path
from graphviz import Source, Digraph

root = Path("/mnt/c/MDU/GNPSnew/shap")
dot_file = root / "ms_shap_tree.dot"
numbered_dot_file = root / "ms_shap_tree_with_numbers.dot"
svg_file = root / "ms_shap_tree_with_numbers.svg"

# Preprocess the .dot file to include node names/numbers in labels
with open(dot_file, "r") as f:
    dot_content = f.read()

'''
# Add node names/numbers to labels
updated_dot_content = []
for line in dot_content.splitlines():
    if "[" in line and "label=" in line:
        node_number = line.split("[")[0].strip()  # Extract node number
        label_start = line.find("label=") + len("label=")
        label_end = line.find("]", label_start)
        label = line[label_start:label_end].strip('"')
        updated_line = line.replace(label, f"{node_number}: {label}")
        updated_dot_content.append(updated_line)
    else:
        updated_dot_content.append(line)

# Save the updated .dot file
with open(numbered_dot_file, "w") as f:
    f.write("\n".join(updated_dot_content))
'''
# Render the updated .dot file
#graph = Source.from_file(numbered_dot_file)
# Create a Digraph object and set attributes
graph = Digraph()
graph.attr(ratio="compress")  # Set graph attributes
graph.attr(rankdir="TB")  # Change layout to top-to-bottom
#graph.graph_attr['nodesep'] = '0.1'  # Reduce spacing between nodes

# Render the updated .dot file
graph.body.extend(dot_content)  # Add the updated .dot content

graph.render(svg_file, format="svg", cleanup=True)

print(f"Graph rendered and saved to {svg_file}")