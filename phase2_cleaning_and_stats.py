import re
import numpy as np
import pandas as pd
from scipy import stats

# 1. Load raw dataset
df = pd.read_csv("Bacteria_dataset_Multiresictance.csv")

# Drop PII and unstructured noisy notes
df = df.drop(columns=["Name", "Email", "Address", "Notes"], errors="ignore")

# 2. Split 'age/gender' into numeric Age and standardized Gender
def parse_age_gender(val):
    if pd.isna(val):
        return np.nan, np.nan
    parts = str(val).strip().split("/")
    age = pd.to_numeric(parts[0], errors="coerce") if len(parts) > 0 else np.nan
    gender = parts[1].strip().upper() if len(parts) > 1 else np.nan
    if gender not in ["M", "F"]:
        gender = np.nan
    if pd.notna(age) and (age < 0 or age > 110):
        age = np.nan
    return age, gender

df[["Age", "Gender"]] = df["age/gender"].apply(lambda x: pd.Series(parse_age_gender(x)))

# 3. Clean 'Souches' (strip leading isolate ID like 'S290 ' and fix species typos)
def clean_species(val):
    if pd.isna(val):
        return np.nan
    cleaned = re.sub(r"^S\d+\s+", "", str(val).strip())
    return cleaned if cleaned.lower() not in ["missing", "?", "nan", ""] else np.nan

species_typo_map = {
    "E.coi": "Escherichia coli",
    "E.cli": "Escherichia coli",
    "E. coli": "Escherichia coli",
    "Enter.bacteria spp.": "Enterobacteria spp.",
    "Enteobacteria spp.": "Enterobacteria spp.",
    "Klbsiella pneumoniae": "Klebsiella pneumoniae",
    "Klebsie.lla pneumoniae": "Klebsiella pneumoniae",
    "Protus mirabilis": "Proteus mirabilis",
    "Proeus mirabilis": "Proteus mirabilis",
    "Prot.eus mirabilis": "Proteus mirabilis",
}

df["Species"] = df["Souches"].apply(clean_species).replace(species_typo_map)

# 4. Standardize clinical binary columns & Infection_Freq
binary_map = {"yes": 1, "true": 1, "1": 1, "no": 0, "false": 0, "0": 0}
for col in ["Diabetes", "Hypertension", "Hospital_before"]:
    df[col] = df[col].astype(str).str.strip().str.lower().map(binary_map)

df["Infection_Freq"] = pd.to_numeric(df["Infection_Freq"], errors="coerce")

# 5. Standardize all 15 Antibiotic Susceptibility Test (AST) columns
ast_cols = [
    "AMX/AMP", "AMC", "CZ", "FOX", "CTX/CRO", "IPM", "GEN", "AN",
    "Acide nalidixique", "ofx", "CIP", "C", "Co-trimoxazole", "Furanes", "colistine"
]

ast_map = {
    "r": "R", "s": "S", "i": "I", "intermediate": "I",
    "?": np.nan, "missing": np.nan, "nan": np.nan
}

for col in ast_cols:
    df[col] = df[col].astype(str).str.strip().str.lower().map(ast_map)

# 6. Feature Engineering: Resistance Count, Tested Count, and Resistance Rate (%)
df["Tested_Drugs"] = df[ast_cols].notnull().sum(axis=1)
df["Resistance_Count"] = (df[ast_cols] == "R").sum(axis=1)
df["Resistance_Rate_Pct"] = np.where(
    df["Tested_Drugs"] > 0,
    (df["Resistance_Count"] / df["Tested_Drugs"]) * 100,
    np.nan
)

# Drop rows with 0 tested drugs (completely empty clinical records)
df_clean = df[df["Tested_Drugs"] > 0].copy()
df_clean.to_csv("cleaned_amr_data.csv", index=False)

print(f"=== CLEANED DATASET SHAPE: {df_clean.shape[0]} rows x {df_clean.shape[1]} cols ===\n")

print("=== STANDARDIZED BACTERIAL SPECIES ===")
print(df_clean["Species"].value_counts(dropna=False), "\n")

print("=== UNIT 1: LOCATION & VARIABILITY METRICS ===")
num_cols = ["Age", "Infection_Freq", "Resistance_Count", "Resistance_Rate_Pct"]
stats_rows = []
for col in num_cols:
    s = df_clean[col].dropna()
    n = len(s)
    mean_val = s.mean()
    med_val = s.median()
    trim_mean = stats.trim_mean(s, 0.1)  # 10% trimmed mean
    std_val = s.std(ddof=1)
    var_val = s.var(ddof=1)
    q1, q3 = s.quantile(0.25), s.quantile(0.75)
    iqr_val = q3 - q1
    se_val = std_val / np.sqrt(n)
    ci_low, ci_high = mean_val - 1.96 * se_val, mean_val + 1.96 * se_val
    stats_rows.append({
        "Variable": col, "N": n, "Mean": round(mean_val, 2),
        "Median": round(med_val, 2), "TrimMean_10%": round(trim_mean, 2),
        "SD": round(std_val, 2), "Variance": round(var_val, 2),
        "IQR": round(iqr_val, 2), "95%_CI": f"[{ci_low:.2f}, {ci_high:.2f}]"
    })
print(pd.DataFrame(stats_rows).to_string(index=False), "\n")

print("=== UNIT 1: COVARIANCE MATRIX (NUMERIC FEATURES) ===")
print(df_clean[num_cols].cov().round(2), "\n")

print("=== UNIT 1: PEARSON CORRELATION MATRIX ===")
print(df_clean[num_cols + ["Diabetes", "Hypertension", "Hospital_before"]].corr().round(3))