# Sharing CascadeScope with colleagues

Two ways to let colleagues try the app without cloning it and installing R.

## Option 1 — Static in-browser bundle (no install)

Compile the app to a self-contained static site with
[shinylive](https://posit-dev.github.io/shinylive/); it runs R in the browser
via webR, so colleagues need only a modern browser. First load fetches the R
packages (~30-60 MB) and is slow (20-40 s), then caches.

From R, in this folder:

```r
install.packages("shinylive")   # once
shinylive::export(".", "site")  # builds ./site
```

Preview it locally before sharing:

```bash
cd site && python3 -m http.server 8000   # open http://localhost:8000
```

Then host the `site/` folder on any static host and share the link:

- **Netlify Drop** (no account): drag the `site` folder onto <https://app.netlify.com/drop>.
- **GitHub Pages:** put `site/` on a `gh-pages` branch (or in `docs/`), enable Pages.
- **Any web server / S3 bucket:** upload `site/`, serve `index.html`.

Opening `index.html` from disk (file://) will not work; it must be served over http(s).
Building the bundle needs the system `libarchive` library (for the R `archive` /
`shinylive` packages).

## Option 2 — shinyapps.io (a durable hosted link)

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name="<you>", token="<token>", secret="<secret>")
rsconnect::deployApp(".", appFiles = c(
  "app.R", "ici_locked_reference.R", "score_cohort.R",
  "references/ICI_v1.rds", "example_cohort.csv", "www/logo.svg", "www/icon.svg"))
```

This runs on a real R server, so every feature works exactly as it does locally.
(`make_app_assets.R` is left out on purpose: it only rebuilds the reference.)

## What colleagues can do

- **Patient report:** pick a patient (or upload a one-sample file) and read the
  archetype call, the standardized-score gauge, the cascade drill-down, and the
  per-patient 98-gene cascade tree.
- **Cohort explorer:** load a cohort (CSV / TSV / Excel, gene x sample) or the
  bundled example; see the archetype donut, the layer-score distributions, a
  click-through beeswarm, the table, and a cohort PDF.
- Use **Load example cohort** to try it instantly with no file.

**Research use only.** The app is investigational and ships only public example
data (an IMvigor210 subset); it contains no patient data.
