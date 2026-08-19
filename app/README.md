# Prototype application

A static, dependency-free page that carries the fitted AMPLIFY surrogate and
evaluates it in the browser. No build step, no framework, no network calls.

| File | Role |
| --- | --- |
| `index.html` | Page structure and all narrative text. |
| `styles.css` | Single stylesheet; light and dark themes follow the browser. |
| `model.js` | The forward pass, mirroring `masked_count_net` in `R/torch_model.R`. |
| `app.js` | Input handling, flagging, and the inline SVG response curve. |
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

## Editing the wording

`site.json` holds the app name, tagline, attribution line, partner
institutions, repository link, flag thresholds, and scenario labels. Edit it and
redeploy; no code changes are required. The flag thresholds are:

```json
"thresholds": { "outbreak": 1.1, "self_limiting": 0.9 }
```

Estimates above `outbreak` are flagged as outbreak potential, below
`self_limiting` as self-limiting, and in between as borderline.

## Deployment

`.github/workflows/pages.yml` validates `model.json` and `site.json`, then
publishes this directory to GitHub Pages on every push to `main` that touches
`app/`. Pages must be enabled for the repository with **Source: GitHub Actions**.
