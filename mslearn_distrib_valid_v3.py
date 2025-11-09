import argparse
import pandas as pd
import numpy as np
from dataclasses import dataclass, asdict
import json
from pathlib import Path
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import xgboost as xgb
from xgboost import XGBClassifier
from sklearn.utils.class_weight import compute_class_weight
import joblib

# Custom Voting Classifier
class VoteClassifier(BaseEstimator, ClassifierMixin):
    def __init__(self, boundary: int, column=None):
        super().__init__()
        self.boundary = boundary
        self.single_column = column

    def fit(self, X=None, y=None):
        # No fitting necessary for the vote
        return self

    def predict(self, X):
        # Apply the vote function to each row
        return pd.DataFrame(X).apply(self.vote, axis=1)

    def vote(self, row):
        if self.single_column:
            return row[self.single_column]
        # Count the number of positive values in the row
        positive = (row > 0).sum()
        # If enough columns have a positive value, return 1 (or True)
        if positive >= self.boundary:
            return 1
        else:
            return 0

    def score(self, X, y):
        # Predict the values
        y_pred = self.predict(X)
        # Calculate the accuracy score
        return accuracy_score(y, y_pred)

@dataclass(frozen=True)
class InputParameters:
    id: int
    train_path: str
    test_path: str
    classification_mode: str
    fitting_set: str
    classifier: str

# Function to compare two InputParameters instances, ignoring the id field
def same_input_parameters(param1: InputParameters, param2: InputParameters) -> bool:
    return (
        param1.train_path == param2.train_path and
        param1.test_path == param2.test_path and
        param1.classification_mode == param2.classification_mode and
        param1.fitting_set == param2.fitting_set and
        param1.classifier == param2.classifier
    )

# Function to iteratively remove highly correlated columns and return the list of dropped columns
def remove_highly_correlated(df, threshold=0.7):
    to_drop = []
    while True:
        corr_matrix = df.corr().abs()
        upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        max_corr = upper.max().max()
        if max_corr < threshold:
            break
        # Find the column with the highest correlation
        to_remove = upper.stack().idxmax()[1]
        to_drop.append(to_remove)
        df = df.drop(columns=[to_remove])
    return df, to_drop

def get_all_desc_all_pred(train_df, test_df, y_train, y_test):
    return (
        train_df.drop(columns=['reference_spectra', 'smiles', 'canonical_smiles', 'is_primary', 'iceberg_probability']),
        y_train,
        test_df.drop(columns=['reference_spectra', 'smiles', 'canonical_smiles', 'is_primary', 'iceberg_probability']),
        y_test
    )

def get_all_pred(train_df, test_df, y_train, y_test):
    return (
        train_df[['cfmid_spectra', 'rassp_spectra', 'scarf_spectra', 'massformer_spectra', 'iceberg_spectra']],
        y_train,
        test_df[['cfmid_spectra', 'rassp_spectra', 'scarf_spectra', 'massformer_spectra', 'iceberg_spectra']],
        y_test
    )

def get_all_desc(train_df, test_df, y_train, y_test):
    return (
        train_df.drop(columns=['cfmid_spectra', 'rassp_spectra', 'scarf_spectra', 'massformer_spectra', 'iceberg_spectra', 'iceberg_probability', 'reference_spectra', 'smiles', 'canonical_smiles', 'is_primary']),
        y_train,
        test_df.drop(columns=['cfmid_spectra', 'rassp_spectra', 'scarf_spectra', 'massformer_spectra', 'iceberg_spectra', 'iceberg_probability', 'reference_spectra', 'smiles', 'canonical_smiles', 'is_primary']),
        y_test
    )

def get_all_desc_all_pred_uncorr(train_df, test_df, y_train, y_test):
    train_df_all_desc_all_pred, to_drop_all_desc_all_pred = remove_highly_correlated(train_df.drop(columns=['reference_spectra','smiles','canonical_smiles','is_primary','iceberg_probability']))
    to_drop_all_desc_all_pred.extend(['reference_spectra','smiles','canonical_smiles','is_primary','iceberg_probability'])
    test_df_all_desc_all_pred = test_df.drop(columns=to_drop_all_desc_all_pred)
    return (
        train_df_all_desc_all_pred,
        y_train,
        test_df_all_desc_all_pred,
        y_test
    )

def get_all_desc_uncorr(train_df, test_df, y_train, y_test):
    train_df_all_desc, to_drop_all_desc = remove_highly_correlated(train_df.drop(columns=['cfmid_spectra','rassp_spectra','scarf_spectra','massformer_spectra','iceberg_spectra','iceberg_probability','reference_spectra','smiles','canonical_smiles','is_primary']))
    to_drop_all_desc.extend(['cfmid_spectra','rassp_spectra','scarf_spectra','massformer_spectra','iceberg_spectra','iceberg_probability','reference_spectra','smiles','canonical_smiles','is_primary'])
    test_df_all_desc = test_df.drop(columns=to_drop_all_desc)
    return (
        train_df_all_desc,
        y_train,
        test_df_all_desc,
        y_test
    )

# Generate json template file for the input parameters

def generate_input_file():
    # Path to the training and testing data files
    workdir = Path('/mnt/c/MDU')
    train_path_list = [
        workdir / 'GNPSnew/GNPSnew_canonical_output_split_train70_binary_1000lines.txt',
        workdir / 'GNPSnew/GNPSnew_canonical_output_split_train70_bin4_1000lines.txt'
    ]
    test_path_list = [
        workdir / 'GNPSnew/GNPSnew_canonical_output_split_test15_binary_1000lines.txt',
        workdir / 'GNPSnew/GNPSnew_canonical_output_split_test15_bin4_1000lines.txt'
    ]
    classification_modes = [
        'binary',
        'binned'
    ]
    parameter_set_list = []
    id = 1
    for (train_path, test_path, mode) in zip(train_path_list, test_path_list, classification_modes):
        for fitting_set in fitting_sets.keys():
            for classifier in classifier_dict.keys():
                if 'Vote-' in classifier and (fitting_set != 'AllPred' or mode != 'binary'):
                    continue
                parameter_set_list.append(InputParameters(id, str(train_path), str(test_path), mode, fitting_set, classifier))
                id += 1
    with open(json_file_path, 'w') as f:
        json.dump([asdict(p) for p in parameter_set_list], f, indent=4)

def load_list_from_json(filename: str) -> list[InputParameters]:
    # Load the list from the JSON file
    with open(filename, 'r') as f:
        data_list = json.load(f)
    
    # Create a list to store unique SpectraData instances
    param_list = []
    unique_ids = set()
    
    for item in data_list:
        # Create a SpectraData instance
        param_data = InputParameters(**item)
        
        # Ensure no duplicate IDs
        if param_data.id in unique_ids:
            raise ValueError(f"Duplicate ID detected: {param_data.id}")
        
        # Ensure no duplicate items
        if any(same_input_parameters(param_data, existing_param) for existing_param in param_list):
            raise ValueError(f"Duplicate item detected: {param_data}")
        
        param_list.append(param_data)
        unique_ids.add(param_data.id)
    
    return param_list

def compute_metrics_per_molecule(test_df, y_true, y_pred, classification_mode, output_path):
    # Combine predictions and canonical_smiles
    df = test_df[["canonical_smiles"]].copy()
    df["true_label"] = y_true.values
    df["predicted_label"] = y_pred

    grouped = df.groupby("canonical_smiles")
    results = []

    for mol, group in grouped:
        if classification_mode == "binary":
            average = "binary"
        else:
            average = "weighted"

        precision = precision_score(group["true_label"], group["predicted_label"], average=average, zero_division=0)
        recall = recall_score(group["true_label"], group["predicted_label"], average=average, zero_division=0)
        accuracy = accuracy_score(group["true_label"], group["predicted_label"])
        f1 = f1_score(group["true_label"], group["predicted_label"], average=average, zero_division=0)

        results.append({
            "canonical_smiles": mol,
            "n_samples": len(group),
            "accuracy": accuracy,
            "precision": precision,
            "recall": recall,
            "f1": f1
        })

    result_df = pd.DataFrame(results)
    result_df.to_csv(output_path, index=False)
    print(f"Saved per-molecule metrics to: {output_path}")


# Define the classifiers with their parameters
classifier_dict = {
    "Nearest-Neighbors": KNeighborsClassifier(n_neighbors=3),
    "RBF-SVM": SVC(gamma=2, C=1, random_state=42),
    "Gaussian-Process": GaussianProcessClassifier(1.0 * RBF(1.0), random_state=42),
    "Decision-Tree": DecisionTreeClassifier(max_depth=5, random_state=42),
    "Random-Forest": RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1, random_state=42),
    "Neural-Net": MLPClassifier(alpha=1, max_iter=1000, random_state=42),
    "AdaBoost": AdaBoostClassifier(random_state=42),
    "Naive-Bayes": GaussianNB(),
    "QDA": QuadraticDiscriminantAnalysis(),
    "XGBoost": xgb.XGBClassifier(tree_method="hist", early_stopping_rounds=3),
    "Vote-2": VoteClassifier(2),
    "Vote-cfmid_spectra": VoteClassifier(1, column='cfmid_spectra'),
    "Vote-rassp_spectra": VoteClassifier(1, column='rassp_spectra'),
    "Vote-scarf_spectra": VoteClassifier(1, column='scarf_spectra'),
    "Vote-massformer_spectra": VoteClassifier(1, column='massformer_spectra'),
    "Vote-iceberg_spectra": VoteClassifier(1, column='iceberg_spectra'),
    "XGBoost_Paramset449_valid": XGBClassifier(objective='binary:logistic', base_score=None, booster=None, colsample_bylevel=None, colsample_bynode=None, colsample_bytree=1.0, device=None, eval_metric='auc', gamma=1, grow_policy=None, interaction_constraints=None, learning_rate=0.1, max_bin=None, max_cat_threshold=None, max_cat_to_onehot=None, max_delta_step=None, max_depth=7, max_leaves=None, min_child_weight=None, monotone_constraints=None, multi_strategy=None, n_jobs=None, num_parallel_tree=None, random_state=None, reg_alpha=1, reg_lambda=1, sampling_method=None, scale_pos_weight=5.3, subsample=0.85, tree_method='hist', validate_parameters=None, verbosity=0, num_class=None, use_label_encoder=False, early_stopping_rounds=3),
    "XGBoost_Paramset1025_valid": XGBClassifier(objective='multi:softprob', base_score=None, booster=None, colsample_bylevel=None, colsample_bynode=None, colsample_bytree=0.85, device=None, eval_metric='mlogloss', gamma=0, grow_policy=None, interaction_constraints=None, learning_rate=0.1, max_bin=None, max_cat_threshold=None, max_cat_to_onehot=None, max_delta_step=None, max_depth=7, max_leaves=None, min_child_weight=None, monotone_constraints=None, multi_strategy=None, n_jobs=None, num_parallel_tree=None, random_state=None, reg_alpha=0, reg_lambda=1, sampling_method=None, scale_pos_weight=1.0, subsample=0.7, tree_method='hist', validate_parameters=None, verbosity=0, num_class=4, use_label_encoder=False, early_stopping_rounds=3),
    "XGBoost_Paramset1541_valid": XGBClassifier(objective='multi:softprob', base_score=None, booster=None, colsample_bylevel=None, colsample_bynode=None, colsample_bytree=1.0, device=None, eval_metric='mlogloss', gamma=1, grow_policy=None, interaction_constraints=None, learning_rate=0.1, max_bin=None, max_cat_threshold=None, max_cat_to_onehot=None, max_delta_step=None, max_depth=7, max_leaves=None, min_child_weight=None, monotone_constraints=None, multi_strategy=None, n_jobs=None, num_parallel_tree=None, random_state=None, reg_alpha=1, reg_lambda=0.5, sampling_method=None, scale_pos_weight=1.0, subsample=1.0, tree_method='hist', validate_parameters=None, verbosity=0, num_class=4, use_label_encoder=False, early_stopping_rounds=3),
    "XGBoost_Paramset656_valid": XGBClassifier(objective='binary:logistic', base_score=None, booster=None, colsample_bylevel=None, colsample_bynode=None, colsample_bytree=1.0, device=None, eval_metric='auc', gamma=0, grow_policy=None, interaction_constraints=None, learning_rate=0.1, max_bin=None, max_cat_threshold=None, max_cat_to_onehot=None, max_delta_step=None, max_depth=7, max_leaves=None, min_child_weight=None, monotone_constraints=None, multi_strategy=None, n_jobs=None, num_parallel_tree=None, random_state=None, reg_alpha=0, reg_lambda=2, sampling_method=None, scale_pos_weight=5.3, subsample=0.7, tree_method='hist', validate_parameters=None, verbosity=0, num_class=None, use_label_encoder=False, early_stopping_rounds=3)
}

# Define the feature sets to use for fitting the classifiers
fitting_sets = {
    'AllDescAllPred': lambda: get_all_desc_all_pred(train_df, test_df, y_train, y_test),
    'AllPred': lambda: get_all_pred(train_df, test_df, y_train, y_test),
    'AllDesc': lambda: get_all_desc(train_df, test_df, y_train, y_test),
    'AllDescAllPredUncorr': lambda: get_all_desc_all_pred_uncorr(train_df, test_df, y_train, y_test),
    'AllDescUncorr': lambda: get_all_desc_uncorr(train_df, test_df, y_train, y_test)
}


#
# Main
#

# Set up argument parser
parser = argparse.ArgumentParser(description="Generate input file for mslearn.")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-g", "--generate", action="store_true", help="Generate input file")
group.add_argument("-i", "--input-id", type=int, help="Input ID integer (mandatory)")

parser.add_argument("-w", "--workdir", type=str, default=".", help="Working directory (default: current directory)")
parser.add_argument("-j", "--json-file", type=str, default="mslearn_input.json", help="JSON file name (default: mslearn_input.json)")
parser.add_argument("--save-model", action="store_true", help="If set, save trained model to disk")

# Parse arguments
args = parser.parse_args()
# Set variables based on parsed arguments
workdir = Path(args.workdir)

models_dir = workdir / "models"
if args.save_model:
    models_dir.mkdir(exist_ok=True)


json_file_path = workdir / args.json_file
input_id = args.input_id

# Verify that workdir is a valid directory
if not workdir.is_dir():
    print(f"Error: The specified workdir '{workdir}' is not a valid directory.")
    exit(1)

# Verify that json_file is readable

if not args.generate and not json_file_path.is_file():
    print(f"Error: The directory for the specified json file '{json_file_path}' does not exist.")
    exit(1)

# Generate the input file
if args.generate:
    generate_input_file()
    exit(0)

# Read the input file
with open(json_file_path, 'r') as f:
    input_parameters_list = load_list_from_json(json_file_path)

# Find the input parameters with the specified ID
selected_input_parameters = None
for input_parameters in input_parameters_list:
    if input_parameters.id == input_id:
        selected_input_parameters = input_parameters
        break

# Check if the input ID was found
if selected_input_parameters is None:
    print(f"Error: Input ID {input_id} not found in the input file.")
    exit(1)


# Load the DataFrames from CSV files
train_df = pd.read_csv(Path(selected_input_parameters.train_path))
test_df = pd.read_csv(Path(selected_input_parameters.test_path))

# Print the number of unique values in the "canonical_smiles" column
print(f"Number of unique values in 'canonical_smiles' in train_df: {train_df['canonical_smiles'].nunique()}")
print(f"Number of unique values in 'canonical_smiles' in test_df: {test_df['canonical_smiles'].nunique()}")



# Check for NaNs in the training data and drop rows with NaNs
if train_df.isnull().values.any():
    nan_rows_train = train_df[train_df.isnull().any(axis=1)]
    print("Dropped rows with NaNs in training set:")
    print(nan_rows_train.index.tolist())
    train_df = train_df.dropna()
else:
    print("No NaNs found in the training data.")

# Check for NaNs in the testing data and drop rows with NaNs
if test_df.isnull().values.any():
    nan_rows_test = test_df[test_df.isnull().any(axis=1)]
    print("Dropped rows with Nans tesing set:")
    print(nan_rows_test.index.tolist())
    test_df = test_df.dropna()
else:
    print("No NaNs found in the testing data.")


# Drop rows where the values in spectra columns are not 0 or 1
if selected_input_parameters.classification_mode=="binary":
    print(f"Removing rows with non-binary values in predictions and results columns...")
    columns_to_check = ['cfmid_spectra','rassp_spectra','scarf_spectra','massformer_spectra','iceberg_spectra','iceberg_probability','reference_spectra']
    train_df = train_df[train_df[columns_to_check].isin([0, 1]).all(axis=1)]
    test_df = test_df[test_df[columns_to_check].isin([0, 1]).all(axis=1)]

    # Print the number of rows after dropping
    print(f"Number of rows in training after dropping non-binary: {len(train_df)}")
    print(f"Number of rows in test after dropping non-binary: {len(test_df)}")

# Separate the DataFrame into features (X) and target (y)
y_train = train_df['reference_spectra']
y_test = test_df['reference_spectra']

selected_fitting_sets = {selected_input_parameters.fitting_set: fitting_sets[selected_input_parameters.fitting_set]}
selected_classifiers = {selected_input_parameters.classifier: classifier_dict[selected_input_parameters.classifier]}

# Initialize results DataFrame
results_df = pd.DataFrame(columns=["Input", "Model", "Accuracy", "Precision", "Recall", "F1"])

# Fit classifiers and store scores
for input_name, get_fitting_set in selected_fitting_sets.items():
    X_train, y_train, X_test, y_test = get_fitting_set()
    print(f"Fitting classifiers for input: {input_name}...")
    for name, clf in selected_classifiers.items():
        if name.startswith("XGBoost"):
            clf.fit(X_train, y_train, eval_set=[(X_test, y_test)], verbose=False)
        else:
            clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)
        if selected_input_parameters.classification_mode=="binary":
            precision = precision_score(y_test, y_pred, average='binary')
            recall = recall_score(y_test, y_pred, average='binary')
            f1 = f1_score(y_test, y_pred, average='binary')
        else:
            precision = precision_score(y_test, y_pred, average='weighted')
            recall = recall_score(y_test, y_pred, average='weighted')
            f1 = f1_score(y_test, y_pred, average='weighted')

        # Combine X_test, y_test, and y_pred into one DataFrame
        prediction_df = X_test.copy()
        prediction_df["true_label"] = y_test.values
        prediction_df["predicted_label"] = y_pred

        # Save to CSV
        pred_outfile = Path(workdir / f"{input_id}-{name}-{input_name}-predictions.csv")
        prediction_df.to_csv(pred_outfile, index=False)
        print(f"Saved predictions to: {pred_outfile}")
        
        # Create a temporary DataFrame for the current classifier
        temp_df = pd.DataFrame([{
                "Input": input_name,
                "Model": name,
                "Accuracy": accuracy,
                "Precision": precision,
                "Recall": recall,
                "F1": f1
        }])
        
        # Save the temporary DataFrame to a CSV file
        temp_df.to_csv(Path(models_dir / f"{input_id}-{name}-{input_name}-output.csv"), index=False)

        # Save per-molecule metrics
        metrics_outfile = Path(workdir / f"{input_id}-{name}-{input_name}-molecule_metrics.csv")
        compute_metrics_per_molecule(
            test_df.loc[X_test.index],  # get canonical_smiles from aligned indices
            y_test,
            y_pred,
            selected_input_parameters.classification_mode,
            metrics_outfile
        )
        
        # Append the results to the main results DataFrame
        results_df = results_df._append(temp_df, ignore_index=True)
        if args.save_model:
            model_filename = f"{input_id}-{name}-{input_name}.json"
            # model_filename = f"{input_id}-{name}-{input_name}.pkl"
            model_path = models_dir / model_filename
            clf.save_model(model_path)
            # joblib.dump(clf, model_path)
            print(f"Saved model to: {model_path}")

# Create separate pivot DataFrames for each metric
pivot_accuracy_df = results_df.pivot(index="Model", columns="Input", values="Accuracy")
pivot_precision_df = results_df.pivot(index="Model", columns="Input", values="Precision")
pivot_recall_df = results_df.pivot(index="Model", columns="Input", values="Recall")
pivot_f1_df = results_df.pivot(index="Model", columns="Input", values="F1")

# Print each pivot DataFrame
print("\n")
print("Accuracy Results:")
print(pivot_accuracy_df)
print("\nPrecision Results:")
print(pivot_precision_df)
print("\nRecall Results:")
print(pivot_recall_df)
print("\nF1 Results:")
print(pivot_f1_df)

'''
# Save each pivot DataFrame to a separate CSV file
pivot_accuracy_df.to_csv('binary_classifier_accuracy_results.csv')
pivot_precision_df.to_csv('binary_classifier_precision_results.csv')
pivot_recall_df.to_csv('binary_classifier_recall_results.csv')
pivot_f1_df.to_csv('binary_classifier_f1_results.csv')


# Test the single column votes
if "Vote-" in selected_input_parameters.classifier:
    test_columns = ['cfmid_spectra','rassp_spectra','scarf_spectra','massformer_spectra','iceberg_probability']
    x_test = test_df[test_columns]
    print("\n")
    for column in ['cfmid_spectra','rassp_spectra','scarf_spectra','massformer_spectra','iceberg_probability']:
        single_column_clf = VoteClassifier(1, column)
        # no fitting necessary for our custom vote classifier
        print(f"Accuracy for column: {column} is: {single_column_clf.score(test_df, y_test)}")
'''