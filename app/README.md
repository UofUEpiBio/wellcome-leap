# Prototype application

A static, dependency-free page that carries the fitted AMPLIFY surrogate and
evaluates it in the browser. No build step, no framework, no network calls.

| File | Role |
| --- | --- |
| `index.html` | Page structure and all narrative text. |
| `styles.css` | Single stylesheet, light theme, no build step. |
| `model.js` | The forward pass, mirroring `masked_count_net` in `R/torch_model.R`. |
| `app.js` | Tabs, input handling, the feature mapping, flagging, and the response curve. |
| `model.json` | Exported weights and preprocessing constants. **Generated.** |
| `site.json` | Editable presentation settings. |

## Running it locally

```sh
python3 -m http.server --directory app 8000
```

Then open <http://localhost:8000>. Opening `index.html` straight from disk does
not work: the page fetches `model.json` and `site.json` at runtime, and browsers
block that over `file://`.

## Regenerating the model

`model.json` is written from the generated torch weights:

```sh
Rscript scripts/05_export_web_model.R
```

The script checks a `torch`-free forward pass against `torch` before writing, so
a mismatch fails loudly rather than shipping wrong numbers. Re-run it after any
retraining; nothing else in the app needs to change.

## Page layout

The page has two tabs, driven entirely by `app.js`: **Estimator**, which carries
the tool itself and a short statement of what it is, and **How it works**, which
carries the three-piece description, the input mapping, and the scope of the
full version. The tab is reflected in the URL fragment (`#estimator`,
`#method`), so either view can be linked directly.

## Editing the wording

`site.json` holds the app name, tagline, attribution line, partner
institutions, repository link, flag thresholds, organism labels, and the carrier
profile. Edit it and redeploy; no code changes are required. The flag thresholds
are:

```json
"thresholds": { "outbreak": 1.1, "self_limiting": 0.9 }
```

Estimates above `outbreak` are flagged as outbreak potential, below
`self_limiting` as self-limiting, and in between as borderline.

## The two indicators are a presentation layer

The fitted surrogate takes one standardized individual feature and the organism
label, because that is what the simulation produces. The page exposes exactly
those two, one indicator each, so every estimate maps back to a row of the
training data:

| Indicator in the app | Model input | Mapping |
| --- | --- | --- |
| Colonisation density, log₁₀ CFU/g stool | `x` | `z = (density - population_mean) / population_sd`, capped at ±3 |
| Organism identified | `scenario` | organism label to `scenario` id, encoded by `preprocessor.scenario_levels` |

The slider range, the population mean and standard deviation behind the z-score,
and the organism names all live in `site.json`; `app.js` applies them and the
**How it works** tab tabulates the result, so the mapping is visible rather than
implied. It is illustrative — the simulation's feature is standard normal and has
no clinical units, and the organism names stand in for its lower- and
higher-transmission scenarios — but it is a one-to-one relabelling, not a derived
score.

Nothing in this layer touches `model.json`: change the labels or the constants
freely without retraining. Keep the `id` of each entry in `scenarios` matching a
key in the model's `preprocessor.scenario_levels`, or the export check in CI will
fail.
