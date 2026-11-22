import glob
import pandas as pd

# Automatically find the CSV file in the current directory
csv_files = glob.glob("*.csv")
if not csv_files:
    raise FileNotFoundError("No CSV file found in the current folder.")

file_path = csv_files[0]
print(f"Loading dataset: {file_path}\n")

df = pd.read_csv(file_path)

print("=== 1. SHAPE & DUPLICATES ===")
print(f"Rows: {df.shape[0]} | Columns: {df.shape[1]}")
print(f"Duplicate rows: {df.duplicated().sum()}\n")

print("=== 2. COLUMNS, DTYPES & MISSING VALUES ===")
summary = pd.DataFrame({
    "Dtype": df.dtypes,
    "Non-Null": df.notnull().sum(),
    "Missing": df.isnull().sum(),
    "Missing_%": (df.isnull().mean() * 100).round(2),
    "Unique_Vals": df.nunique()
})
print(summary.to_string(), "\n")

print("=== 3. SAMPLE UNIQUE VALUES (CATEGORICAL / LOW-CARDINALITY) ===")
for col in df.columns:
    if df[col].dtype == "object" or df[col].nunique() <= 15:
        print(f"- {col}: {df[col].dropna().unique()[:10]}")

print("\n=== 4. FIRST 5 ROWS ===")
print(df.head().to_string())