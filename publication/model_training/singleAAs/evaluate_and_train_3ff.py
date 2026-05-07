"""Train the packaged isolated-amino-acid multi-forcefield Fdewet model.

This script uses the processed isolated-amino-acid training tables in this
directory and writes the model artifact used by HydroMap:

    hydromap/models/Fdewet_isolated_aa_multi_forcefield.joblib
"""

from __future__ import annotations

from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
FEATURES_CSV = SCRIPT_DIR / "all_features.csv"
TARGET_CSV = SCRIPT_DIR / "INDUS_dewet_FE_singleAAs.csv"
OUTPUT_MODEL = REPO_ROOT / "hydromap" / "models" / "Fdewet_isolated_aa_multi_forcefield.joblib"
OUTPUT_FIGURE = SCRIPT_DIR / "predicted_actual_3ff.png"

SELECTED_FEATURES = ["U_pw", "PC1_tri", "bulk_80-90_tri"]
AA_LABELS = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def load_training_data() -> tuple[pd.DataFrame, pd.Series]:
    features_df = pd.read_csv(FEATURES_CSV, index_col=[0, 1])
    target_df = pd.read_csv(TARGET_CSV, index_col=[0, 1])

    features_df = features_df.rename(
        columns={
            "total_pot": "U_pw",
            "PC1": "PC1_tri",
            "bulk_80-90": "bulk_80-90_tri",
        }
    )

    return features_df[SELECTED_FEATURES], target_df["free_energy"]


def leave_one_out_predictions(features: pd.DataFrame, target: pd.Series) -> tuple[list[float], list[float]]:
    scaled_features = StandardScaler().fit_transform(features)
    model = LinearRegression()

    actual_values: list[float] = []
    predicted_values: list[float] = []
    for train_idx, test_idx in LeaveOneOut().split(scaled_features):
        x_train = scaled_features[train_idx]
        x_test = scaled_features[test_idx]
        y_train = target.values[train_idx]
        y_test = target.values[test_idx]

        model.fit(x_train, y_train)
        y_pred = model.predict(x_test)

        actual_values.append(float(y_test[0]))
        predicted_values.append(float(y_pred[0]))

    return actual_values, predicted_values


def print_metrics(actual_values: list[float], predicted_values: list[float]) -> None:
    actual = np.array(actual_values)
    predicted = np.array(predicted_values)
    rmse = np.sqrt(np.mean((actual - predicted) ** 2))

    print(f"LOOCV R^2: {r2_score(actual_values, predicted_values):.3f}")
    print(f"LOOCV Pearson r: {pearsonr(actual_values, predicted_values)[0]:.3f}")
    print(f"LOOCV Spearman r: {spearmanr(actual_values, predicted_values)[0]:.3f}")
    print(f"LOOCV RMSE: {rmse:.3f} kJ/mol")


def save_loocv_plot(features: pd.DataFrame, actual_values: list[float], predicted_values: list[float]) -> None:
    amino_acids = features.index.get_level_values(0)
    forcefields = features.index.get_level_values(1)
    color_map = {"a03ws": "red", "a99SBdisp": "blue", "C36m": "green"}

    plt.figure(figsize=(6, 6))
    for amino_acid, forcefield, actual, predicted in zip(amino_acids, forcefields, actual_values, predicted_values):
        plt.scatter(actual, predicted, color=color_map[forcefield], s=40, alpha=0.5, label=forcefield)
        plt.text(
            actual,
            predicted + 0.02,
            AA_LABELS[amino_acid],
            ha="center",
            va="bottom",
            fontsize=12,
            color=color_map[forcefield],
        )

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), fontsize=14)

    plt.plot([2.8, 5.6], [2.8, 5.6], "k--", linewidth=1)
    plt.xlim(2.8, 5.6)
    plt.ylim(2.8, 5.6)
    plt.title(r"$F_{dewet}$ Predicted vs Actual: isolated amino acids (3 force fields)", fontsize=18)
    plt.xlabel("Actual Free Energy (kJ/mol)", fontsize=15)
    plt.ylabel("Predicted Free Energy (kJ/mol)", fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    actual_array = np.array(actual_values)
    predicted_array = np.array(predicted_values)
    rmse = np.sqrt(np.mean((actual_array - predicted_array) ** 2))
    performance_text = (
        f"$R^2$ = {r2_score(actual_values, predicted_values):.2f}\n"
        f"Spearman's r = {spearmanr(actual_values, predicted_values)[0]:.2f}\n"
        f"RMSE = {rmse:.2f} kJ/mol"
    )
    plt.text(0.7, 0.1, performance_text, ha="center", va="bottom", fontsize=15, transform=plt.gca().transAxes)
    plt.tight_layout()
    plt.savefig(OUTPUT_FIGURE, dpi=300)
    plt.close()


def train_final_model(features: pd.DataFrame, target: pd.Series) -> None:
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)

    final_model = LinearRegression()
    final_model.fit(scaled_features, target)

    print("Final model coefficients:")
    for feature, coefficient in zip(SELECTED_FEATURES, final_model.coef_):
        print(f"{feature}: {coefficient:.4f}")

    OUTPUT_MODEL.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": final_model,
            "scaler": scaler,
            "feature_order": SELECTED_FEATURES,
        },
        OUTPUT_MODEL,
    )
    print(f"Final model + scaler saved to '{OUTPUT_MODEL}'")


def main() -> None:
    print("Selected features:", SELECTED_FEATURES)
    features, target = load_training_data()
    actual_values, predicted_values = leave_one_out_predictions(features, target)
    print_metrics(actual_values, predicted_values)
    save_loocv_plot(features, actual_values, predicted_values)
    train_final_model(features, target)


if __name__ == "__main__":
    main()
