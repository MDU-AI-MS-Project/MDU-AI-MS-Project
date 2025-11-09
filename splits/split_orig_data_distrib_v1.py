import pandas as pd
import argparse

def bin_and_save(df, columns, bin_type, filename):
    df_binned = df.copy(deep=True)

    if bin_type == "binary":
        for col in columns:
            assert all(isinstance(value, float) for value in df_binned[col]), f"{col} must contain only floats"
            df_binned[col] = df_binned[col].apply(lambda x: 1 if x > 0.0 else 0)

    elif bin_type == "bin4":
        bin_width = 0.33
        for col in columns:
            assert all(isinstance(value, float) for value in df_binned[col]), f"{col} must contain only floats"
            df_binned[col] = df_binned[col].apply(lambda x: min(int(x / bin_width), 3) if x > 0.0 else 0)

    else:
        raise ValueError("Unsupported bin_type")

    df_binned.to_csv(filename, index=False)
    return df_binned

def summarize_from_df(name, df, label):
    size = len(df)
    dist = df["reference_spectra"].value_counts(normalize=True).sort_index()
    print(f"{name} set ({label}) size: {size}")
    for k, v in dist.items():
        print(f"  class {int(k)}: {v:.2%}")
    print()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("split_id", type=int, help="Split number for seeding randomness")
    args = parser.parse_args()
    split_id = args.split_id
    seed = split_id

    data = pd.read_csv("/RG/compbio/michaelel/data/GNPSnew_canonical_output_split2.csv")
    columns = ['iceberg_spectra', 'iceberg_probability', 'scarf_spectra', 'cfmid_spectra', 'rassp_spectra', 'massformer_spectra', 'reference_spectra']

    # Step 1: Validation split
    non_zeros = data[data["reference_spectra"] > 0]
    zeros = data[data["reference_spectra"] == 0]
    val_size = int(0.15 * len(data))
    val_nz = val_size // 2
    val_z = val_size - val_nz
    val_set = pd.concat([
        non_zeros.sample(n=val_nz, random_state=seed),
        zeros.sample(n=val_z, random_state=seed + 1)
    ])
    train_set = data.drop(val_set.index)
    train_set = train_set.sample(frac=1, random_state=seed + 2).reset_index(drop=True)
    val_set = val_set.sample(frac=1, random_state=seed + 3).reset_index(drop=True)

    val_path = f"/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split_valid15_split{split_id}.csv"
    val_set.to_csv(val_path, index=False)

    # Step 2: Test split from train
    test_size = int(0.15 * len(train_set))
    train_nz = train_set[train_set["reference_spectra"] > 0]
    train_z = train_set[train_set["reference_spectra"] == 0]
    test_nz = test_size // 2
    test_z = test_size - test_nz
    test_set = pd.concat([
        train_nz.sample(n=test_nz, random_state=seed + 4),
        train_z.sample(n=test_z, random_state=seed + 5)
    ])
    final_train_set = train_set.drop(test_set.index)
    final_train_set = final_train_set.sample(frac=1, random_state=seed + 6).reset_index(drop=True)
    test_set = test_set.sample(frac=1, random_state=seed + 7).reset_index(drop=True)

    # Save raw splits
    train_path = f"/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split_train70_split{split_id}.csv"
    test_path = f"/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split_test15_split{split_id}.csv"
    final_train_set.to_csv(train_path, index=False)
    test_set.to_csv(test_path, index=False)

    # Step 3: Binning and capture both versions
    binary_val = bin_and_save(val_set, columns, "binary", f"/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split_valid15_binary_split{split_id}.csv")
    bin4_val   = bin_and_save(val_set, columns, "bin4",   f"/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split_valid15_4bin_split{split_id}.csv")
    
    binary_test = bin_and_save(test_set, columns, "binary", f"/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split_test15_binary_split{split_id}.csv")
    bin4_test   = bin_and_save(test_set, columns, "bin4",   f"/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split_test15_bin4_split{split_id}.csv")
    
    binary_train = bin_and_save(final_train_set, columns, "binary", f"/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split_train70_binary_split{split_id}.csv")
    bin4_train   = bin_and_save(final_train_set, columns, "bin4",   f"/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split_train70_bin4_split{split_id}.csv")

    # Step 4: Summary from in-memory binned data
    print(f"\n=== Summary for split {split_id} ===\n")

    summarize_from_df("Train", binary_train, "binary")
    summarize_from_df("Train", bin4_train,   "bin4")

    summarize_from_df("Validation", binary_val, "binary")
    summarize_from_df("Validation", bin4_val,   "bin4")

    summarize_from_df("Test", binary_test, "binary")
    summarize_from_df("Test", bin4_test,   "bin4")

if __name__ == "__main__":
    main()
