import pandas as pd
import matplotlib.pyplot as plt

def plot_3d(filename):
    # Load the CSV file
    file_path = filename
    data = pd.read_csv(file_path)
    
    # Specify the columns to plot
    columns_to_plot = ['massformer 100', 'ms-pred scarf', 'ms-pred iceberg']
    
    # Set the transparency of the bars (0.0 to 1.0, where 1.0 is completely opaque)
    bar_transparency = 1
    # Set the bar width
    bar_width = 5
    
    
    # Extract the first column title
    first_column_title = data.columns[0]
    
    #define threshold
    threshold=10
    
    #remove all values below threshold
    for column in columns_to_plot:
        data[column] = data[column].apply(lambda x: x if x >= threshold else 0)
        
    
    # Create a 3D bar graph
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Number of bars
    num_bars = len(data[data.columns[0]])
    
    # Positions for bars on x and y axes
    x_pos = range(num_bars)
    y_pos = range(len(columns_to_plot))
    
    # Create bars
    for i, column in enumerate(columns_to_plot):
        ax.bar(data[first_column_title], data[column], zs=i, zdir='y', alpha=bar_transparency,width=bar_width, label=column)
    
    # Set labels and title
    ax.set_xlabel('m/z')
    #ax.set_ylabel('Columns')
    ax.set_zlabel('Intensity')
    ax.set_title(first_column_title)
    
    # Set y-axis ticks to column names
    ax.set_yticks(y_pos)
    ax.set_yticklabels(columns_to_plot)
    
    # Show legend
    ax.legend()
    
    # Show the plot
    plt.show()


plot_3d('masteroutput.csv')