"""Sanity-check the committed model export before it is published.

Confirms that the JSON the browser needs is present, that every dense layer is
shaped consistently, and that the presentation settings carry the fields the
page reads. Uses only the standard library so the workflow needs no setup.
"""

import json
import sys

LAYERS = ("hidden_1", "hidden_2", "output")
SITE_KEYS = ("app_name", "tagline", "attribution", "thresholds", "scenarios")


def check_model(path):
    with open(path) as handle:
        model = json.load(handle)

    previous_units = None
    for name in LAYERS:
        layer = model["layers"][name]
        units, inputs = layer["units"], layer["inputs"]
        if len(layer["weight"]) != units * inputs:
            raise SystemExit(f"{name}: weight has {len(layer['weight'])} values, expected {units * inputs}")
        if len(layer["bias"]) != units:
            raise SystemExit(f"{name}: bias has {len(layer['bias'])} values, expected {units}")
        if previous_units is not None and inputs != previous_units:
            raise SystemExit(f"{name}: expects {inputs} inputs but the previous layer emits {previous_units}")
        previous_units = units
    if previous_units != 1:
        raise SystemExit("The output layer must emit a single value.")

    preprocessor = model["preprocessor"]
    if preprocessor["x_sd"] <= 0:
        raise SystemExit("The feature standard deviation must be positive.")
    if not preprocessor["scenario_levels"]:
        raise SystemExit("No scenario encoding was exported.")
    print(f"{path}: {len(LAYERS)} layers, scenarios {sorted(preprocessor['scenario_levels'])}")
    return set(preprocessor["scenario_levels"])


def check_site(path, model_scenarios):
    with open(path) as handle:
        site = json.load(handle)
    missing = [key for key in SITE_KEYS if key not in site]
    if missing:
        raise SystemExit(f"{path} is missing: {', '.join(missing)}")
    scenario_ids = {scenario["id"] for scenario in site["scenarios"]}
    unknown = scenario_ids - model_scenarios
    if unknown:
        raise SystemExit(f"{path} offers scenarios the model cannot encode: {', '.join(sorted(unknown))}")
    print(f"{path}: scenarios {sorted(scenario_ids)}")


if __name__ == "__main__":
    check_site(sys.argv[2], check_model(sys.argv[1]))
