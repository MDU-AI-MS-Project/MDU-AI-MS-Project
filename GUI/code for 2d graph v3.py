
import pandas as pd
import matplotlib.pyplot as plt

def plot_2d(filename):
    # Load the CSV file
    file_path = filename
    data = pd.read_csv(file_path)
    
    # Specify the columns to plot
    columns_to_plot = ['massformer 100', 'CFM-ID 20', 'CFM-ID 40']
    
    #define threshold
    threshold=5
    
    #remove all values below threshold
    for column in columns_to_plot:
        data[column] = data[column].apply(lambda x: x if x >= threshold else 0)
    
    # Extract the first column title
    first_column_title = data.columns[0]
    
    # Set the transparency of the bars (0.0 to 1.0, where 1.0 is completely opaque)
    bar_transparency = 0.3
    
    # Create a bar graph
    plt.figure(figsize=(10, 6))
    for column in columns_to_plot:
        plt.bar(data[first_column_title], data[column], alpha=bar_transparency, label=column)
    
    # Add labels and title
    plt.xlabel('m/z')
    plt.ylabel('Intensity')
    plt.title(data.columns[0])
    #plt.title('Bar Graph of Selected Columns')
    plt.legend()
    
    # Show the plot
    plt.show()

plot_2d('masteroutput.csv')