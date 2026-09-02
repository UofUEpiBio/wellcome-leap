# Prototype application

A static, dependency-free page with two functional views: a transparent
multiscale pARG mechanism explorer and the original fitted AMPLIFY surrogate.
Both run in the browser. There is no build step, framework, server component,
or external network call.

| File | Role |
| --- | --- |
| `index.html` | Page structure and all narrative text. |
| `styles.css` | Single stylesheet, light theme, no build step. |
| `model.js` | The forward pass, mirroring `masked_count_net` in `R/torch_model.R`. |
| `app.js` | Shared tabs plus legacy-surrogate inputs, flagging, and response curve. |
| `model.json` | Exported weights and preprocessing constants. **Generated.** |
| `site.json` | Editable presentation settings. |
| `multiscale-model.js` | Dependency-free ODE and reproduction-number calculations mirroring the R implementation. |
| `multiscale-app.js` | Controls and rendering for the multiscale explorer. |
| `multiscale.json` | Deliberately small synthetic mechanism, site, and intervention configuration. |
| `ml-diagram.jpeg` | Architecture figure shown on the **How it works** tab. |

## Running it locally

```sh
python3 -m http.server --directory app 8000
```

Then open <http://localhost:8000>. Opening `index.html` straight from disk does
not work: the page fetches its JSON configuration at runtime, and browsers block
that over `file://`.

## Multiscale explorer

The default tab evaluates the same simplified equations as `R/within_host.R`
and `R/reproduction.R`. It reports four targets:

- within-host reference and effective R, from a seed-normalized first-generation
  matrix of successful pARG acquisition on new bacterial backgrounds; and
- between-host reference and effective R, from the time integral of the
  within-host resistant fraction, contact rate, establishment, and current
  susceptible fraction.

Reference values remove antibiotic exposure and set the between-host
susceptible fraction to one. The explorer does not calculate population-level R
or R(t). It is synthetic and explanatory, not a clinical decision tool. The R
pipeline remains the source for simulation, fitting, ABM experiments, native-R
torch emulation, and rendered reports; the browser evaluates the mechanism
directly so no generated multiscale data or weights need to be committed.

`multiscale.json` mirrors `default_multiscale_config()` and
`multiscale_site_table()`. A test guards their shared values. Keep the JavaScript
equations synchronized when the corresponding R equations change.

## Regenerating the model

`model.json` is written from the generated torch weights:

```sh
Rscript scripts/05_export_web_model.R
```

The script checks a `torch`-free forward pass against `torch` before writing, so
a mismatch fails loudly rather than shipping wrong numbers. Re-run it after any
retraining; nothing else in the app needs to change.

## Page layout

The page has three tabs. **Multiscale explorer** is the default mechanistic
view, **Legacy surrogate** preserves the originally published two-input model,
and **How it works** explains how omics, the ODE, ML, and the ABM connect. The tab
is reflected in the URL fragment (`#multiscale`, `#estimator`, `#method`), so any
view can be linked directly.

`ml-diagram.jpeg` is a compressed copy of `fig/ml-diagram.jpeg`, kept inside
`app/` because the Pages workflow publishes this directory alone. Regenerate it
after editing the source figure:

```sh
sips -s format jpeg -s formatOptions 72 --resampleWidth 1400 \
  fig/ml-diagram.jpeg --out app/ml-diagram.jpeg
```

Keep the `width` and `height` attributes on the `<img>` in step with whatever
that writes, so the page reserves the right space before the image loads.

## Editing the wording

`site.json` holds the app name, tagline, attribution line, partner
institutions, repository link, flag thresholds, organism labels, and the
indicator mapping. Edit it and redeploy; no code changes are required. The flag
thresholds are:

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
