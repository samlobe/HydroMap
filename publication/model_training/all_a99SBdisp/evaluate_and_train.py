"""Train the packaged HydroMap Fdewet model.

This script uses the processed training table in this directory and writes the
model artifact used by HydroMap:

    hydromap/models/Fdewet.joblib
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
INPUT_CSV = SCRIPT_DIR / "all_features.csv"
OUTPUT_MODEL = REPO_ROOT / "hydromap" / "models" / "Fdewet.joblib"
OUTPUT_FIGURE = SCRIPT_DIR / "protein_context_LOOCV.png"

SELECTED_FEATURES = ["U_pw", "40-50_tri", "140-150_tri"]
KT_KJMOL_AT_300K = 0.008314 * 300.0


def load_training_data() -> tuple[pd.DataFrame, pd.Series]:
    features_df = pd.read_csv(INPUT_CSV, index_col=[0, 1])

    features_df.index = features_df.index.set_levels(
        features_df.index.levels[0].str.replace("singleAA", "capped AA"),
        level=0,
    )
    features_df = features_df.rename(
        columns={
            "total_pot": "U_pw",
            "40-50": "40-50_tri",
            "140-150": "140-150_tri",
        }
    )

    category_renames = {
        "SARS2": "SARS-CoV-2 spike",
        "alpha": "alpha variant",
        "beta": "beta variant",
        "gamma": "gamma variant",
        "omicron": "omicron variant",
    }
    for old_name, new_name in category_renames.items():
        features_df.index = features_df.index.set_levels(
            features_df.index.levels[0].str.replace(old_name, new_name),
            level=0,
        )

    target = features_df["INDUS_dewet_FE"] * KT_KJMOL_AT_300K
    return features_df[SELECTED_FEATURES], target


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
    categories = features.index.get_level_values(0)
    residues = features.index.get_level_values(1)
    color_map = {
        "capped AA": "blue",
        "hydrophobin": "gray",
        "SARS-CoV-2 spike": "green",
        "alpha variant": "red",
        "beta variant": "cyan",
        "gamma variant": "magenta",
        "omicron variant": "gold",
    }

    plt.figure(figsize=(6, 6))
    for category, residue, actual, predicted in zip(categories, residues, actual_values, predicted_values):
        plt.scatter(actual, predicted, color=color_map[category], s=40, alpha=0.5, label=category)
        plt.text(actual, predicted + 0.06, residue, ha="center", va="bottom", fontsize=12, color=color_map[category])

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), fontsize=14)

    plt.plot([0, 10], [0, 10], "k--", linewidth=1)
    plt.xlim(3.5, 9)
    plt.ylim(3.5, 9)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel(r"Actual $F_{dewet}$ (kJ/mol)", fontsize=15)
    plt.ylabel(r"Predicted $F_{dewet}$ (kJ/mol)", fontsize=15)
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
