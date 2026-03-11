"""
simulate_sequences.py
---------------------
Generates synthetic patient sequences that mimic the real data format used in
the comorbidity network analysis pipeline.

Output files (in --output_dir):
  sequence_main_simulated.csv      — main disease events
  sequence_parallel_simulated.csv  — parallel risk-factor events
  simulation_summary.txt           — simulation metadata

Usage:
  python simulate_sequences.py --n_patients 5000 --seed 42 --output_dir data_simulated

Column format (both output files):
  ID                    — patient identifier
  age_in_days_since_born — age at event occurrence (days)
  event_name            — MeSH disease or risk-factor name
"""

import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import truncnorm

# ---------------------------------------------------------------------------
# 1.  Disease definitions
# ---------------------------------------------------------------------------

# Real-data incidence rates (proportion of patients with ≥1 event)
DISEASE_INCIDENCE = {
    "Hypertension":                          0.963,
    "Dyslipidemias":                         0.752,
    "Coronary Disease":                      0.479,
    "Diabetes Mellitus":                     0.363,
    "Atrial Fibrillation":                   0.300,
    "Acute Kidney Injury":                   0.198,
    "Cerebrovascular Disorders":             0.190,
    "Kidney Diseases, Chronic, Stage 3":     0.175,
    "Heart Failure":                         0.171,
    "Peripheral Vascular Diseases":          0.086,
    "Pulmonary Embolism":                    0.085,
    "Mitral Valve Insufficiency":            0.075,
    "Pacemaker, Artificial":                 0.097,
    "Stroke, Hemorrhagic":                   0.032,
    "Heart Arrest":                          0.030,
    "Tachycardia, Ventricular":              0.020,
    "Kidney Diseases, Chronic, Stage 4":     0.021,
    "Kidney Failure, Chronic, Endstage":     0.012,
    "Cardiomyopathy, Dilated":               0.014,
    "Kidney Diseases, Chronic, Stage 5":     0.014,
    "Tricuspid Valve Insufficiency":         0.008,
    "Stroke":                                0.030,
    "Heart Failure, Endstage":               0.003,
    "Mitral Valve Stenosis":                 0.005,
}

# Mean onset age (years) per disease — used to centre truncated Normal draws
DISEASE_ONSET_MEAN_YR = {
    "Hypertension":                          55.0,
    "Dyslipidemias":                         55.0,
    "Coronary Disease":                      62.0,
    "Diabetes Mellitus":                     57.0,
    "Atrial Fibrillation":                   65.0,
    "Acute Kidney Injury":                   67.0,
    "Cerebrovascular Disorders":             65.0,
    "Kidney Diseases, Chronic, Stage 3":     64.0,
    "Heart Failure":                         68.0,
    "Peripheral Vascular Diseases":          63.0,
    "Pulmonary Embolism":                    63.0,
    "Mitral Valve Insufficiency":            62.0,
    "Pacemaker, Artificial":                 70.0,
    "Stroke, Hemorrhagic":                   67.0,
    "Heart Arrest":                          68.0,
    "Tachycardia, Ventricular":              66.0,
    "Kidney Diseases, Chronic, Stage 4":     67.0,
    "Kidney Failure, Chronic, Endstage":     70.0,
    "Cardiomyopathy, Dilated":               60.0,
    "Kidney Diseases, Chronic, Stage 5":     69.0,
    "Tricuspid Valve Insufficiency":         65.0,
    "Stroke":                                66.0,
    "Heart Failure, Endstage":               72.0,
    "Mitral Valve Stenosis":                 60.0,
}

# Disease ordering constraints: (prerequisite, consequence)
# The consequence must occur AFTER the prerequisite if both are sampled.
DISEASE_ORDER_CONSTRAINTS = [
    ("Kidney Diseases, Chronic, Stage 3", "Kidney Diseases, Chronic, Stage 4"),
    ("Kidney Diseases, Chronic, Stage 4", "Kidney Diseases, Chronic, Stage 5"),
    ("Kidney Diseases, Chronic, Stage 4", "Kidney Failure, Chronic, Endstage"),
    ("Heart Failure",                      "Heart Failure, Endstage"),
]

# Standard deviation for onset age draws (years)
ONSET_STD_YR = 8.0
# Age bounds (years)
AGE_MIN_YR, AGE_MAX_YR = 35.0, 85.0

# ---------------------------------------------------------------------------
# 2.  Risk-factor definitions (parallel sequence)
# ---------------------------------------------------------------------------

# Each entry: (event_name, base_probability)
# Static factors are assigned at a baseline age draw; dynamic ones near disease events.
STATIC_RISK_FACTORS = [
    # Sex
    ("Female", 0.50),
    ("Male",   0.50),
    # Ethnicity (approximately mutually exclusive — handled specially)
    # Family history
    ("Mother Diabetes Disease",        0.20),
    ("Mother Heart Disease",           0.25),
    ("Mother Stroke Disease",          0.12),
    ("Father Diabetes Disease",        0.18),
    ("Father Heart Disease",           0.28),
    ("Father Stroke Disease",          0.10),
    ("Siblings Diabetes Disease",      0.15),
    ("Siblings Heart Disease",         0.20),
    ("Siblings Stroke Disease",        0.08),
    # BMI
    ("Body Mass Index, Low",           0.05),
    ("Body Mass Index, Medium",        0.45),
    ("Body Mass Index, High",          0.50),
    # Socioeconomic
    ("Deprivation Index, Low",         0.33),
    ("Deprivation Index, Middle",      0.34),
    ("Deprivation Index, High",        0.33),
    ("Employment Status, Active",      0.50),
    ("Employment Status, Inactive",    0.35),
    ("Employment Status, No",          0.15),
    ("Educational Attainment, High",   0.30),
    ("Educational Attainment, Medium", 0.40),
    ("Educational Attainment, Low",    0.30),
    ("Marital Status, Married",        0.55),
    ("Marital Status, Single",         0.25),
    ("Marital Status, Divorced",       0.10),
    ("Marital Status, Widowed",        0.10),
    ("Living Alone",                   0.30),
    ("Not Living Alone",               0.70),
    ("Smoking, Yes",                   0.25),
    ("Smoking, No",                    0.75),
]

ETHNICITY_OPTIONS = [
    "White", "Asian", "Black", "Mixed", "Other"
]
ETHNICITY_PROBS = [0.87, 0.06, 0.03, 0.02, 0.02]

# Dynamic risk factors: assigned near disease events
DYNAMIC_RISK_FACTORS = [
    # Blood pressure
    ("Systolic Blood Pressure, High",        0.45),
    ("Systolic Blood Pressure, Normal",      0.55),
    ("Diastolic Blood Pressure, High",       0.35),
    ("Diastolic Blood Pressure, Normal",     0.65),
    # Biomarkers
    ("HbA1c, High",                          0.30),
    ("HbA1c, Normal",                        0.70),
    ("BNP, High",                            0.25),
    ("BNP, Normal",                          0.75),
    ("NT-proBNP, High",                      0.22),
    ("NT-proBNP, Normal",                    0.78),
    ("Hemoglobin, Low",                      0.20),
    ("Hemoglobin, Normal",                   0.80),
    ("Platelet Count, Low",                  0.10),
    ("Platelet Count, Normal",               0.90),
    ("White Blood Cell Count, High",         0.12),
    ("White Blood Cell Count, Normal",       0.88),
    # Ejection fraction
    ("Ejection Fraction, Reduced",           0.30),
    ("Ejection Fraction, Mid-range",         0.20),
    ("Ejection Fraction, Preserved",         0.50),
    # Heart rate
    ("Heart Rate, High",                     0.15),
    ("Heart Rate, Normal",                   0.75),
    ("Heart Rate, Low",                      0.10),
    # Comorbidities
    ("Chronic Obstructive Pulmonary Disease", 0.15),
    ("Deep Vein Thrombosis",                 0.08),
    ("Anemia",                               0.18),
    ("Thyroid Diseases",                     0.12),
    ("Sepsis",                               0.07),
    ("Pneumonia",                            0.10),
    ("HIV Infections",                       0.02),
    ("Liver Diseases",                       0.08),
    ("Rheumatic Diseases",                   0.07),
    ("Peptic Ulcer",                         0.05),
    ("Urinary Tract Infections",             0.10),
    ("Alcoholic Cardiomyopathy",             0.04),
    ("Amyloidosis",                          0.02),
    ("Bradycardia",                          0.08),
    ("Pericarditis",                         0.03),
    ("Myocarditis",                          0.03),
    ("Heart Block, First Degree",            0.08),
    ("Heart Block, Complete",                0.03),
    ("Thrombocytopenia",                     0.07),
    ("Hypertension Control, Good",           0.50),
    ("Hypertension Control, Bad",            0.50),
    ("Diabetes Control, Good",               0.55),
    ("Diabetes Control, Bad",                0.45),
]

# ---------------------------------------------------------------------------
# 3.  Helper functions
# ---------------------------------------------------------------------------

def _truncated_normal_days(mean_yr: float, std_yr: float,
                            low_yr: float, high_yr: float,
                            rng: np.random.Generator) -> int:
    """Draw an age in days from a truncated Normal distribution."""
    mean_d = mean_yr * 365
    std_d  = std_yr  * 365
    low_d  = low_yr  * 365
    high_d = high_yr * 365
    a = (low_d  - mean_d) / std_d
    b = (high_d - mean_d) / std_d
    # Draw using scipy, but with a numpy rng via the bits workaround
    sample = truncnorm.rvs(a, b, loc=mean_d, scale=std_d,
                           random_state=int(rng.integers(0, 2**31)))
    return int(np.clip(round(sample), low_d, high_d))


def _sample_age_entry(rng: np.random.Generator) -> int:
    """Draw a patient's age at study entry (40–70 years) in days."""
    return int(rng.integers(40 * 365, 70 * 365))


# ---------------------------------------------------------------------------
# 4.  Per-patient simulation
# ---------------------------------------------------------------------------

def simulate_patient_main(pid: int, rng: np.random.Generator) -> list[dict]:
    """Return list of main-sequence disease event dicts for one patient."""
    events = []

    # Sample which diseases occur for this patient
    present_diseases = {}
    for disease, incidence in DISEASE_INCIDENCE.items():
        if rng.random() < incidence:
            onset_days = _truncated_normal_days(
                DISEASE_ONSET_MEAN_YR[disease], ONSET_STD_YR,
                AGE_MIN_YR, AGE_MAX_YR, rng
            )
            present_diseases[disease] = onset_days

    # Enforce ordering constraints
    for prereq, consequence in DISEASE_ORDER_CONSTRAINTS:
        if prereq in present_diseases and consequence in present_diseases:
            if present_diseases[consequence] <= present_diseases[prereq]:
                gap_days = int(rng.integers(180, 3 * 365))  # 6 months to 3 years
                present_diseases[consequence] = present_diseases[prereq] + gap_days
        elif consequence in present_diseases and prereq not in present_diseases:
            # Remove consequence if prerequisite is absent
            del present_diseases[consequence]

    # Generate 1–4 events per disease (repeat visits)
    for disease, first_onset in present_diseases.items():
        n_events = int(rng.choice([1, 2, 3, 4], p=[0.55, 0.25, 0.12, 0.08]))
        current_age = first_onset
        for _ in range(n_events):
            if current_age > AGE_MAX_YR * 365:
                break
            events.append({
                "ID":                    pid,
                "age_in_days_since_born": current_age,
                "event_name":            disease,
            })
            # Gap until next event: Exponential with mean 2 years
            gap = int(rng.exponential(2 * 365))
            current_age += gap

    return events


def simulate_patient_parallel(pid: int, main_events: list[dict],
                               rng: np.random.Generator) -> list[dict]:
    """Return list of parallel risk-factor event dicts for one patient."""
    events = []
    baseline_age = _sample_age_entry(rng)

    # --- Sex (mutually exclusive) ---
    sex = rng.choice(["Female", "Male"])
    events.append({
        "ID": pid,
        "age_in_days_since_born": baseline_age,
        "event_name": sex,
    })

    # --- Ethnicity (mutually exclusive) ---
    ethnicity = rng.choice(ETHNICITY_OPTIONS, p=ETHNICITY_PROBS)
    events.append({
        "ID": pid,
        "age_in_days_since_born": baseline_age,
        "event_name": ethnicity,
    })

    # --- BMI (mutually exclusive) ---
    bmi = rng.choice(
        ["Body Mass Index, Low", "Body Mass Index, Medium", "Body Mass Index, High"],
        p=[0.05, 0.45, 0.50]
    )
    events.append({"ID": pid, "age_in_days_since_born": baseline_age, "event_name": bmi})

    # --- IMD (mutually exclusive) ---
    imd = rng.choice(
        ["Deprivation Index, Low", "Deprivation Index, Middle", "Deprivation Index, High"],
        p=[0.33, 0.34, 0.33]
    )
    events.append({"ID": pid, "age_in_days_since_born": baseline_age, "event_name": imd})

    # --- Employment (mutually exclusive) ---
    emp = rng.choice(
        ["Employment Status, Active", "Employment Status, Inactive", "Employment Status, No"],
        p=[0.50, 0.35, 0.15]
    )
    events.append({"ID": pid, "age_in_days_since_born": baseline_age, "event_name": emp})

    # --- Education (mutually exclusive) ---
    edu = rng.choice(
        ["Educational Attainment, High", "Educational Attainment, Medium",
         "Educational Attainment, Low"],
        p=[0.30, 0.40, 0.30]
    )
    events.append({"ID": pid, "age_in_days_since_born": baseline_age, "event_name": edu})

    # --- Marital status (mutually exclusive) ---
    mar = rng.choice(
        ["Marital Status, Married", "Marital Status, Single",
         "Marital Status, Divorced", "Marital Status, Widowed"],
        p=[0.55, 0.25, 0.10, 0.10]
    )
    events.append({"ID": pid, "age_in_days_since_born": baseline_age, "event_name": mar})

    # --- Living situation (mutually exclusive) ---
    liv = rng.choice(["Living Alone", "Not Living Alone"], p=[0.30, 0.70])
    events.append({"ID": pid, "age_in_days_since_born": baseline_age, "event_name": liv})

    # --- Smoking ---
    smk = rng.choice(["Smoking, Yes", "Smoking, No"], p=[0.25, 0.75])
    events.append({"ID": pid, "age_in_days_since_born": baseline_age, "event_name": smk})

    # --- Bernoulli static risk factors ---
    bernoulli_static = [
        ("Mother Diabetes Disease",   0.20),
        ("Mother Heart Disease",      0.25),
        ("Mother Stroke Disease",     0.12),
        ("Father Diabetes Disease",   0.18),
        ("Father Heart Disease",      0.28),
        ("Father Stroke Disease",     0.10),
        ("Siblings Diabetes Disease", 0.15),
        ("Siblings Heart Disease",    0.20),
        ("Siblings Stroke Disease",   0.08),
    ]
    for rf_name, prob in bernoulli_static:
        if rng.random() < prob:
            events.append({
                "ID": pid,
                "age_in_days_since_born": baseline_age,
                "event_name": rf_name,
            })

    # --- Dynamic risk factors: assigned near each disease event ---
    for ev in main_events:
        age = ev["age_in_days_since_born"]
        jitter = int(rng.integers(-90, 90))  # ±3 months
        rf_age = int(np.clip(age + jitter, AGE_MIN_YR * 365, AGE_MAX_YR * 365))

        for rf_name, prob in DYNAMIC_RISK_FACTORS:
            if rng.random() < prob:
                events.append({
                    "ID": pid,
                    "age_in_days_since_born": rf_age,
                    "event_name": rf_name,
                })

    return events


# ---------------------------------------------------------------------------
# 5.  Main pipeline
# ---------------------------------------------------------------------------

def simulate(n_patients: int, seed: int, output_dir: str,
             min_disease_events: int = 2) -> None:
    os.makedirs(output_dir, exist_ok=True)
    rng = np.random.default_rng(seed)

    all_main     = []
    all_parallel = []
    kept = 0
    attempted = 0

    print(f"Simulating sequences for {n_patients} patients (seed={seed})…")

    while kept < n_patients:
        attempted += 1
        pid = attempted  # use attempt number as raw ID

        main_events = simulate_patient_main(pid, rng)

        # Filter: require at least `min_disease_events` distinct diseases
        distinct = len({ev["event_name"] for ev in main_events})
        if distinct < min_disease_events:
            continue

        # Remove events at age 0 (shouldn't occur with current code, but mirrors
        # the real preprocessing step)
        main_events = [ev for ev in main_events if ev["age_in_days_since_born"] > 0]
        if not main_events:
            continue

        parallel_events = simulate_patient_parallel(pid, main_events, rng)

        all_main.extend(main_events)
        all_parallel.extend(parallel_events)
        kept += 1

        if kept % 500 == 0:
            print(f"  {kept}/{n_patients} patients generated …")

    # Build DataFrames and sort chronologically
    df_main = (
        pd.DataFrame(all_main)
        .sort_values(["ID", "age_in_days_since_born"])
        .reset_index(drop=True)
    )
    df_parallel = (
        pd.DataFrame(all_parallel)
        .sort_values(["ID", "age_in_days_since_born"])
        .reset_index(drop=True)
    )

    # Save
    main_path = os.path.join(output_dir, "sequence_main_simulated.csv")
    par_path  = os.path.join(output_dir, "sequence_parallel_simulated.csv")
    df_main.to_csv(main_path,     index=False)
    df_parallel.to_csv(par_path,  index=False)

    # Summary
    n_events_main = len(df_main)
    n_events_par  = len(df_parallel)
    disease_counts = df_main.groupby("event_name")["ID"].nunique()
    incidence_pct  = (disease_counts / n_patients * 100).sort_values(ascending=False)

    summary_lines = [
        f"Simulation summary",
        f"==================",
        f"Seed:            {seed}",
        f"Patients kept:   {n_patients}",
        f"Patients tried:  {attempted}",
        f"Main events:     {n_events_main}",
        f"Parallel events: {n_events_par}",
        f"",
        f"Disease incidence (% of patients with ≥1 event):",
    ]
    for d, pct in incidence_pct.items():
        summary_lines.append(f"  {d:<45} {pct:.1f}%")

    summary_path = os.path.join(output_dir, "simulation_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary_lines))

    print(f"\nDone.")
    print(f"  {main_path}  ({n_events_main} rows)")
    print(f"  {par_path}   ({n_events_par} rows)")
    print(f"  {summary_path}")


# ---------------------------------------------------------------------------
# 6.  CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate synthetic comorbidity sequences for pipeline demonstration."
    )
    parser.add_argument("--n_patients",  type=int,   default=5000,
                        help="Number of patients to generate (default: 5000)")
    parser.add_argument("--seed",        type=int,   default=42,
                        help="Random seed for reproducibility (default: 42)")
    parser.add_argument("--output_dir",  type=str,   default="data_simulated",
                        help="Directory for output CSV files (default: data_simulated)")
    parser.add_argument("--min_events",  type=int,   default=2,
                        help="Minimum distinct diseases per patient (default: 2)")
    args = parser.parse_args()

    simulate(
        n_patients=args.n_patients,
        seed=args.seed,
        output_dir=args.output_dir,
        min_disease_events=args.min_events,
    )
