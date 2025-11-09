import shap
import pandas as pd
import numpy as np
import xgboost as xgb
import pickle
import json
import graphviz
import matplotlib.pyplot as plt
from pathlib import Path

root = Path("/mnt/c/MDU/GNPSnew")
datafile = root / "GNPSnew_canonical_output_split_test15_binary_1000lines.txt"
model_file = root / "2107-XGBoost_Paramset449_valid-AllDescAllPred.json"
shap_csv_file = root / "ms_shap_values.csv"
importance_csv_file = root / "ms_shap_importance.csv"
bar_plot_file = root / "ms_shap_bar_plot.png"
force_plot_file = root / "ms_shap_force_plot.html"
beeswarm_plot_file = root / "ms_shap_beeswarm_plot.png"
tree_dot_file = root / "ms_shap_tree.dot"
tree_png_file = root / "ms_shap_tree.png"

df = pd.read_csv(datafile)
test_df = df.drop(columns=['mol_id', 'primary_mol_id', 'reference_spectra', 'smiles', 'canonical_smiles', 'is_primary', 'iceberg_probability'])

if test_df.isnull().values.any():
    nan_rows_test = test_df[test_df.isnull().any(axis=1)]
    print("Dropped rows with Nans tesing set:")
    print(nan_rows_test.index.tolist())
    test_df = test_df.dropna()
else:
    print("No NaNs found in the testing data.")

print(f"Removing rows with non-binary values in predictions and results columns...")
columns_to_check = ['cfmid_spectra','rassp_spectra','scarf_spectra','massformer_spectra','iceberg_spectra']
test_df = test_df[test_df[columns_to_check].isin([0, 1]).all(axis=1)]
print(f"Number of rows in test after dropping non-binary: {len(test_df)}")

# Load the trained model
xgb_model = xgb.Booster()
xgb_model.load_model(model_file)

# Plot the decision tree
tree_dot = xgb.to_graphviz(xgb_model, num_trees=0)
tree_dot.save(tree_dot_file)
# Save the dot file as a PNG image
graph = graphviz.Source(tree_dot)
graphviz.render('dot',format='png', filepath=tree_dot_file, outfile=tree_png_file)

# SHAP TreeExplainer is optimized for XGBoost
explainer = shap.TreeExplainer(xgb_model)
shap_values = explainer.shap_values(test_df)

# Convert to DataFrame for easier analysis
shap_df = pd.DataFrame(shap_values, columns= test_df.columns)
# Save the SHAP DataFrame to a CSV file
shap_df.to_csv(shap_csv_file, index=False)

feature_importance = np.abs(shap_values).mean(axis=0)
importance_df = pd.DataFrame({
    'feature': test_df.columns,
    'importance': feature_importance
}).sort_values(by='importance', ascending=False)
importance_df.to_csv(importance_csv_file, index=False)

# SHAP plots

# Create a SHAP Explanation object for visualization
explanation = shap.Explanation(
    values=shap_values,  # SHAP values for single prediction
    base_values=explainer.expected_value,  # baseline
    data=test_df,  # feature values
    feature_names=test_df.columns.tolist()  # feature names
)

# The bar plot shows the impact of each feature on the model output
# shap.summary_plot(shap_values, test_df, plot_type="bar")
plt.figure(2)
shap.plots.bar(explanation, max_display=14, show=False)
plt.savefig(bar_plot_file, bbox_inches='tight')

# The force plot shows the SHAP values for predictions interactively with various grouping options.
# Note: Force plot require JavaScript to be rendered in a Jupyter notebook or HTML file
shap.initjs()
force_plt = shap.plots.force(explainer.expected_value, shap_values, test_df)
with open(force_plot_file, "w") as f:
    f.write("<html><head>{}</head><body>{}</body></html>".format(shap.getjs(), force_plt.html()))

# The beeswarm plot shows the distribution of the SHAP values of the predictions for each feature,
# sorted by feature impact on the model output (limited to top 10 features by default)

# Calculate mean impact values for the top features
mean_impacts = np.abs(shap_values).mean(axis=0)

# Update feature names to include impact values
updated_feature_names = [
    f"{feature} ({mean_impacts[i]:.2f})"
    #f"{feature} (Impact: {mean_impacts[i]:.2f})"
    for i, feature in enumerate(explanation.feature_names)
]

# Update the explanation object with the new feature names
explanation.feature_names = updated_feature_names

# Generate the beeswarm plot with updated feature labels
plt.figure(3)
shap.plots.beeswarm(explanation, max_display=14, show=False)

# Save the beeswarm plot as an SVG file
beeswarm_plot_svg_file = root / "ms_shap_beeswarm_plot_with_impact_in_labels.svg"
plt.savefig(beeswarm_plot_svg_file, format="svg", bbox_inches='tight')

print(f"Beeswarm plot with impact values in labels saved as SVG to {beeswarm_plot_svg_file}")

# The waterfall plot shows the SHAP values for a single prediction
# Uncomment the following lines to visualize the waterfall plot for a specific prediction
'''
explanation_0 = shap.Explanation(
    values=shap_values[0],  # SHAP values for single prediction
    base_values=explainer.expected_value,  # baseline
    data=test_df.iloc[0],  # feature values
    feature_names=test_df.columns.tolist()  # feature names
)
shap.waterfall_plot(explanation_0)
'''

# The dependence plot shows the relationship between a feature and the SHAP values
# Uncomment the following lines to visualize the dependence plot for a specific feature
'''
plt.figure(4)
shap.dependence_plot("scarf_spectra", shap_values, test_df, interaction_index="massformer_spectra")
'''