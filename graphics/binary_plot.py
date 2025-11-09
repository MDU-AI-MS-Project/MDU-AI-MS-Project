import plotly.graph_objects as go
import pandas as pd
import numpy as np

# Create DataFrame from the provided data
data = [
  {"Model": "QDA", "Metric": "Accuracy", "Value": 0.8998, "Error": 0.0001},
  {"Model": "QDA", "Metric": "F1", "Value": 0.9072, "Error": 0.0001},
  {"Model": "Naive-Bayes", "Metric": "Accuracy", "Value": 0.9029, "Error": 0.0016},
  {"Model": "Naive-Bayes", "Metric": "F1", "Value": 0.9089, "Error": 0.0008},
  {"Model": "Nearest-Neighbors", "Metric": "Accuracy", "Value": 0.9452, "Error": 0.0002},
  {"Model": "Nearest-Neighbors", "Metric": "F1", "Value": 0.9244, "Error": 0.0006},
  {"Model": "Random-Forest", "Metric": "Accuracy", "Value": 0.9463, "Error": 0.0001},
  {"Model": "Random-Forest", "Metric": "F1", "Value": 0.9239, "Error": 0.0003},
  {"Model": "Neural-Net", "Metric": "Accuracy", "Value": 0.9463, "Error": 0.0001},
  {"Model": "Neural-Net", "Metric": "F1", "Value": 0.9245, "Error": 0.0004},
  {"Model": "AdaBoost", "Metric": "Accuracy", "Value": 0.9464, "Error": 0.0001},
  {"Model": "AdaBoost", "Metric": "F1", "Value": 0.9251, "Error": 0.0001},
  {"Model": "RBF-SVM", "Metric": "Accuracy", "Value": 0.9469, "Error": 0.0001},
  {"Model": "RBF-SVM", "Metric": "F1", "Value": 0.9268, "Error": 0.0001},
  {"Model": "Decision-Tree", "Metric": "Accuracy", "Value": 0.9469, "Error": 0.0001},
  {"Model": "Decision-Tree", "Metric": "F1", "Value": 0.9269, "Error": 0.0002},
  {"Model": "XGBoost", "Metric": "Accuracy", "Value": 0.9471, "Error": 0.0001},
  {"Model": "XGBoost", "Metric": "F1", "Value": 0.9267, "Error": 0.0002}
]
pred_bin_plot_data = [
    {'Model': 'AdaBoost', 'Metric': 'Accuracy', 'Value': 0.5031, 'Error': 0.0065},
    {'Model': 'AdaBoost', 'Metric': 'F1', 'Value': 0.0053, 'Error': 0.0019},
    {'Model': 'Decision-Tree', 'Metric': 'Accuracy', 'Value': 0.6326, 'Error': 0.0029},
    {'Model': 'Decision-Tree', 'Metric': 'F1', 'Value': 0.4199, 'Error': 0.0074},
    {'Model': 'Naive-Bayes', 'Metric': 'Accuracy', 'Value': 0.5440, 'Error': 0.0014},
    {'Model': 'Naive-Bayes', 'Metric': 'F1', 'Value': 0.2426, 'Error': 0.0146},
    {'Model': 'Nearest-Neighbors', 'Metric': 'Accuracy', 'Value': 0.6419, 'Error': 0.0128},
    {'Model': 'Nearest-Neighbors', 'Metric': 'F1', 'Value': 0.4756, 'Error': 0.0328},
    {'Model': 'Neural-Net', 'Metric': 'Accuracy', 'Value': 0.6269, 'Error': 0.0049},
    {'Model': 'Neural-Net', 'Metric': 'F1', 'Value': 0.4031, 'Error': 0.0017},
    {'Model': 'QDA', 'Metric': 'Accuracy', 'Value': 0.5563, 'Error': 0.0004},
    {'Model': 'QDA', 'Metric': 'F1', 'Value': 0.2873, 'Error': 0.0176},
    {'Model': 'Random-Forest', 'Metric': 'Accuracy', 'Value': 0.6326, 'Error': 0.0029},
    {'Model': 'Random-Forest', 'Metric': 'F1', 'Value': 0.4199, 'Error': 0.0074},
    {'Model': 'Vote-2', 'Metric': 'Accuracy', 'Value': 0.5589, 'Error': 0.0050},
    {'Model': 'Vote-2', 'Metric': 'F1', 'Value': 0.5487, 'Error': 0.0123},
    {'Model': 'Vote-cfmid_spectra', 'Metric': 'Accuracy', 'Value': 0.5408, 'Error': 0.0033},
    {'Model': 'Vote-cfmid_spectra', 'Metric': 'F1', 'Value': 0.3113, 'Error': 0.0233},
    {'Model': 'Vote-iceberg_spectra', 'Metric': 'Accuracy', 'Value': 0.5435, 'Error': 0.0015},
    {'Model': 'Vote-iceberg_spectra', 'Metric': 'F1', 'Value': 0.2394, 'Error': 0.0135},
    {'Model': 'Vote-massformer_spectra', 'Metric': 'Accuracy', 'Value': 0.5287, 'Error': 0.0019},
    {'Model': 'Vote-massformer_spectra', 'Metric': 'F1', 'Value': 0.4776, 'Error': 0.0095},
    {'Model': 'Vote-rassp_spectra', 'Metric': 'Accuracy', 'Value': 0.5276, 'Error': 0.0006},
    {'Model': 'Vote-rassp_spectra', 'Metric': 'F1', 'Value': 0.3948, 'Error': 0.0090},
    {'Model': 'Vote-scarf_spectra', 'Metric': 'Accuracy', 'Value': 0.4386, 'Error': 0.0073},
    {'Model': 'Vote-scarf_spectra', 'Metric': 'F1', 'Value': 0.5281, 'Error': 0.0103},
    {'Model': 'XGBoost', 'Metric': 'Accuracy', 'Value': 0.6269, 'Error': 0.0049},
    {'Model': 'XGBoost', 'Metric': 'F1', 'Value': 0.4031, 'Error': 0.0017}
]

df = pd.DataFrame(pred_bin_plot_data)

# Sort models by accuracy values (lowest to highest)
accuracy_df = df[df['Metric'] == 'Accuracy'].sort_values('Value')
model_order = accuracy_df['Model'].tolist()

# Prepare data for plotting
accuracy_data = df[df['Metric'] == 'Accuracy'].set_index('Model').reindex(model_order)
f1_data = df[df['Metric'] == 'F1'].set_index('Model').reindex(model_order)
colors_accuracy = ['#4472C4' if model == 'RBF-SVM' else 'crimson' for model in accuracy_data.index]  # Blue for Accuracy, Orange for F1
colors_f1 = ['#FF8C00', '#4472C4']  # Orange for F1, Blue for Accuracy
# Create the figure
fig = go.Figure()

# Add Accuracy bars
fig.add_trace(go.Bar(
    y=accuracy_data.index,
    x=accuracy_data['Value'],
    error_x=dict(type='data', array=accuracy_data['Error'], color='#333333', thickness=2, width=4),
    orientation='h',
    name='Accuracy',
    marker_color= colors_accuracy,  # Blue color
    offsetgroup=1,
    cliponaxis=False
))

# Add F1 bars
fig.add_trace(go.Bar(
    y=f1_data.index,
    x=f1_data['Value'],
    error_x=dict(type='data', array=f1_data['Error'], color='#333333', thickness=2, width=4),
    orientation='h',
    name='F1',
    marker_color='#FF8C00',  # Orange color
    offsetgroup=2,
    cliponaxis=False
))

# Add text annotations at the right end of error bars
for i, model in enumerate(model_order):
    acc_val = accuracy_data.loc[model, 'Value']
    acc_err = accuracy_data.loc[model, 'Error']
    f1_val = f1_data.loc[model, 'Value']
    f1_err = f1_data.loc[model, 'Error']
    
    # Accuracy label at right end of error bar
    fig.add_annotation(
        x = 0.98,  # Fixed position for F1 labels
        #x=acc_val + acc_err + 0.005,
        y=i - 0.22,  # Slight offset down for accuracy bar
        text=f"{acc_val:.4f} ± {acc_err:.4f}",
        showarrow=False,
        font=dict(size=12, color='#000000'),
        xanchor='left'
    )
    
    # F1 label at right end of error bar
    fig.add_annotation(
        #x=f1_val + f1_err + 0.005,
        x = 0.98,  # Fixed position for F1 labels
        y=i + 0.22,  # Slight offset up for F1 bar
        text=f"{f1_val:.4f} ± {f1_err:.4f}",
        showarrow=False,
        font=dict(size=12, color='#000000'),
        xanchor='left'
    )

# Update layout
fig.update_layout(
    autosize=False,
    title='AllPred Model Performance',
    xaxis_title='Perf Score',
    yaxis_title='ML Models',
    barmode='group',
    bargap=0.15,
    bargroupgap=0.1,
    height=1600,
    width=900,
    xaxis=dict(
        range=[0.00, 1.02], 
        showgrid=True, 
        gridcolor='#D3D3D3', 
        gridwidth=1
    ),
    yaxis=dict(
        showgrid=True, 
        gridcolor='#D3D3D3', 
        gridwidth=1
    ),
    legend=dict(orientation='h', yanchor='bottom', y=1.05, xanchor='center', x=0.5)
)

# Save the chart
fig.write_image("model_performance_chart.svg", format="svg",height=800, width=900)
#fig.show()