import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import importlib.util
import threading
import traceback

# Load external plotting functions
def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

module_2d = load_module("code_for_2d_graph_nv6", "/home/greg/MDU_outputs/GUI/code_for_2d_graph_nv6.py")
module_3d = load_module("code_for_3d_graph_nv6", "/home/greg/MDU_outputs/GUI/code_for_3d_graph_nv6.py")
run_models_module = load_module("gui_run_all_nv4", "/home/greg/MDU_outputs/GUI/gui_run_all_nv4.py")

# Save the figures to files
def save_file():
    if not current_fig_2d or not current_fig_3d:
        messagebox.showerror("Save File", "No graph to save!")
        return
    
    save_path_2d = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[
        ("PNG files", "*.png"), ("JPEG files", "*.jpeg"), ("TIFF files", "*.tiff"), ("PDF files", "*.pdf"), ("All files", "*.*")
    ], initialfile="2D_plot.png")
    if save_path_2d:
        current_fig_2d.savefig(save_path_2d)
    
    save_path_3d = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[
        ("PNG files", "*.png"), ("JPEG files", "*.jpeg"), ("TIFF files", "*.tiff"), ("PDF files", "*.pdf"), ("All files", "*.*")
    ], initialfile="3D_plot.png")
    if save_path_3d:
        current_fig_3d.savefig(save_path_3d)
    
    messagebox.showinfo("Save File", "Files saved successfully!")

# Generate and display the graphs
def plot_graph():
    global current_fig_2d, current_fig_3d
    try:
        threshold = float(threshold_entry.get())
    except ValueError:
        messagebox.showerror("Plot Graph", "Please enter a valid threshold value.")
        return
    
    selected_columns = [tool for tool, var in tool_vars.items() if var.get()]
    
    if not selected_columns:
        messagebox.showerror("Plot Graph", "Please select at least one tool/column to plot.")
        return
    
    # Plot 2D graph
    try:
        fig_2d, ax_2d = module_2d.plot_2d(threshold, selected_columns)
        fig_3d, ax_3d = module_3d.plot_3d(threshold, selected_columns)
    except ValueError as e:
        messagebox.showerror("Plot Graph", str(e))
        return
    
    # Display 2D graph
    for widget in graph_frame.winfo_children():
        widget.destroy()
    
    canvas_2d = FigureCanvasTkAgg(fig_2d, master=graph_frame)
    canvas_2d.draw()
    canvas_2d.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    # Display 3D graph
    canvas_3d = FigureCanvasTkAgg(fig_3d, master=graph_frame)
    canvas_3d.draw()
    canvas_3d.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    current_fig_2d = fig_2d  # Update current figure for saving functionality
    current_fig_3d = fig_3d  # Update current figure for saving functionality

def run_models():
    smiles_string = molecule_entry.get()
    if not smiles_string:
        messagebox.showerror("Input Error", "Please enter a valid SMILES string.")
        return

    def run():
        loading_label.config(text="Running models...")
        try:
            print("Running process_smiles function...")
            output_message = run_models_module.process_smiles(smiles_string)
            print("process_smiles completed successfully.")
            loading_label.config(text="Models run successfully.")
            messagebox.showinfo("Output Data", output_message)
        except Exception as e:
            # Capture the full traceback
            error_message = ''.join(traceback.format_exception(None, e, e.__traceback__))
            print(f"Error during processing:\n{error_message}")
            loading_label.config(text=f"Error: {error_message}")
        finally:
            # Ensure the error message doesn't disappear
            pass

    thread = threading.Thread(target=run)
    thread.start()

def on_closing():
    plt.close('all')  # Close all plot windows
    app.destroy()  # Destroy the main application window

app = tk.Tk()
app.title("MS/MS Graph Generator")
app.geometry("1000x800")  # Increase the size of the panel

# Display startup message
messagebox.showinfo("Startup Message", "Start Docker Desktop as root or with sudo to use CFM-ID tool")

# Create a canvas and scrollbars
canvas = tk.Canvas(app)
scroll_y = tk.Scrollbar(app, orient="vertical", command=canvas.yview)
scroll_x = tk.Scrollbar(app, orient="horizontal", command=canvas.xview)
scrollable_frame = tk.Frame(canvas)

scrollable_frame.bind(
    "<Configure>",
    lambda e: canvas.configure(
        scrollregion=canvas.bbox("all")
    )
)

canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
canvas.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)

# Pack the canvas and scrollbars
canvas.pack(side="left", fill="both", expand=True)
scroll_y.pack(side="right", fill="y")
scroll_x.pack(side="bottom", fill="x")

# Common style for labels and buttons
title_font = ("Arial", 35, "bold")
button_font = ("Arial", 25, "bold")
label_font = ("Arial", 20)
input_font = ("Arial", 20)
button_bg = "#A9A9A9"

# Define tool_vars dictionary
tools = ["CFM-ID 10", "CFM-ID 20", "CFM-ID 40", "ms-pred scarf", "ms-pred iceberg", 
         "massformer 20", "massformer 40", "massformer 60", "massformer 80", "rassp"]
tool_vars = {tool: tk.BooleanVar() for tool in tools}

# Frame for Molecule MS/MS calculation
frame1 = tk.Frame(scrollable_frame, relief=tk.RAISED, borderwidth=2)
frame1.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
tk.Label(frame1, text="Molecule MS/MS Calculation", font=title_font).pack(pady=10)

molecule_frame = tk.Frame(frame1)
molecule_frame.pack(pady=5)
tk.Label(molecule_frame, text="Molecule Name (SMILES):", font=input_font).pack(side=tk.TOP, padx=5)
molecule_entry = tk.Entry(molecule_frame, font=input_font)
molecule_entry.pack(side=tk.TOP, padx=5)

run_button = tk.Button(frame1, text="Run Models", font=button_font, bg=button_bg, command=run_models)
run_button.pack(pady=5)

loading_label = tk.Label(frame1, text="", font=input_font)
loading_label.pack(pady=5)

# Frame for Tools
frame2 = tk.Frame(scrollable_frame, relief=tk.RAISED, borderwidth=2)
frame2.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
tk.Label(frame2, text="Tools Selection", font=title_font).pack(pady=10)

threshold_frame = tk.Frame(frame2)
threshold_frame.pack(pady=5)
tk.Label(threshold_frame, text="Threshold Value:", font=input_font).pack(side=tk.TOP, padx=5)
threshold_entry = tk.Entry(threshold_frame, font=input_font)
threshold_entry.pack(side=tk.TOP, padx=5)

tools_frame = tk.Frame(frame2)
tools_frame.pack(pady=5)

# Frame for CFM-ID Tools
cfm_id_frame = tk.Frame(tools_frame, bg="green", borderwidth=2, relief=tk.RAISED)
cfm_id_frame.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.BOTH, expand=True)
tk.Label(cfm_id_frame, text="CFM-ID", font=label_font, bg="green").pack(pady=5)
cfm_id_tools = ["CFM-ID 10", "CFM-ID 20", "CFM-ID 40"]
for tool in cfm_id_tools:
    tk.Checkbutton(cfm_id_frame, text=tool, variable=tool_vars[tool], font=input_font, bg="white").pack(anchor=tk.W, padx=5, pady=2)

# Frame for MS-PRED Tools
mspred_frame = tk.Frame(tools_frame, bg="blue", borderwidth=2, relief=tk.RAISED)
mspred_frame.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.BOTH, expand=True)
tk.Label(mspred_frame, text="MS-PRED", font=label_font, bg="blue").pack(pady=5)
mspred_tools = ["ms-pred scarf", "ms-pred iceberg"]
for tool in mspred_tools:
    tk.Checkbutton(mspred_frame, text=tool, variable=tool_vars[tool], font=input_font, bg="white").pack(anchor=tk.W, padx=5, pady=2)

# Frame for MassFormer
massformer_frame = tk.Frame(tools_frame, bg="red", borderwidth=2, relief=tk.RAISED)
massformer_frame.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.BOTH, expand=True)
tk.Label(massformer_frame, text="MassFormer", font=label_font, bg="red").pack(pady=5)
massformer_tools = ["massformer 20", "massformer 40", "massformer 60", "massformer 80"]
for tool in massformer_tools:
    tk.Checkbutton(massformer_frame, text=tool, variable=tool_vars[tool], font=input_font, bg="white").pack(anchor=tk.W, padx=5, pady=2)

# Frame for Rassp
rassp_frame = tk.Frame(tools_frame, bg="yellow", borderwidth=2, relief=tk.RAISED)
rassp_frame.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.BOTH, expand=True)
tk.Label(rassp_frame, text="Rassp", font=label_font, bg="yellow").pack(pady=5)
rassp_tools = ["rassp"]
for tool in rassp_tools:
    tk.Checkbutton(rassp_frame, text=tool, variable=tool_vars[tool], font=input_font, bg="white").pack(anchor=tk.W, padx=5, pady=2)

plot_button = tk.Button(frame2, text="Plot Graph", command=plot_graph, font=button_font, bg=button_bg)
plot_button.pack(pady=10)

# Frame for Graph Display
frame3 = tk.Frame(scrollable_frame, relief=tk.RAISED, borderwidth=2)
frame3.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

graph_frame = tk.Frame(frame3)
graph_frame.pack(fill=tk.BOTH, expand=True)

save_button = tk.Button(frame3, text="Save Graph", command=save_file, font=button_font, bg=button_bg)
save_button.pack(side=tk.BOTTOM, pady=10)

current_fig_2d = None
current_fig_3d = None

app.protocol("WM_DELETE_WINDOW", on_closing)  # Handle window close event

app.mainloop()
