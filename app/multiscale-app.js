/* Browser controller for the synthetic multiscale mechanism explorer. */
(function () {
  "use strict";

  const state = {
    config: null,
    surrogate: null,
    site: null,
    intervention: "baseline",
    h: 0,
    gamma: 0,
    delta: 0,
    antibioticDays: 0,
    contactRate: 0,
    susceptibleFraction: 0,
    timer: null,
    mlTimer: null,
    useQuantitative: true,
    useGenomic: true,
    useClinical: true
  };

  const el = function (id) { return document.getElementById(id); };

  function loadJson(path) {
    return fetch(path, { cache: "no-cache" }).then(function (response) {
      if (!response.ok) throw new Error(path + ": " + response.status);
      return response.json();
    });
  }

  Promise.all([loadJson("multiscale.json"), loadJson("model.json")])
    .then(function (loaded) {
      state.config = loaded[0];
      state.surrogate = loaded[1].multiscale;
      if (!state.surrogate) throw new Error("model.json has no multiscale export.");
      buildSiteChoice();
      configureSliders();
      buildInterventionChoice();
      bindControls();
      initializeSurrogate();
      chooseSite(state.config.defaults.site);
    })
    .catch(function (error) {
      el("computing-status").textContent = "Model unavailable";
      el("mechanism-results").setAttribute("aria-busy", "false");
      console.error(error);
    });

  function buildSiteChoice() {
    state.config.sites.forEach(function (site) {
      ["site-select", "ml-site-select"].forEach(function (id) {
        const option = document.createElement("option");
        option.value = site.id;
        option.textContent = site.label;
        el(id).appendChild(option);
      });
    });
  }

  function configureSliders() {
    const controls = [
      ["h", "h-slider"],
      ["gamma", "gamma-slider"],
      ["delta", "delta-slider"],
      ["antibiotic_days", "antibiotic-slider"],
      ["contact_rate", "contact-slider"],
      ["susceptible_fraction", "susceptible-slider"]
    ];
    controls.forEach(function (item) {
      const range = state.config.ranges[item[0]];
      const slider = el(item[1]);
      slider.min = range.min;
      slider.max = range.max;
      slider.step = range.step;
    });
  }

  function buildInterventionChoice() {
    const host = el("intervention-choice");
    state.config.interventions.forEach(function (intervention) {
      const button = document.createElement("button");
      button.type = "button";
      button.role = "radio";
      button.dataset.intervention = intervention.id;
      button.textContent = intervention.label;
      button.addEventListener("click", function () {
        state.intervention = intervention.id;
        updateInterventionChoice();
        scheduleCalculation();
      });
      host.appendChild(button);
      const option = document.createElement("option");
      option.value = intervention.id;
      option.textContent = intervention.label;
      el("ml-intervention-select").appendChild(option);
    });
  }

  function bindControls() {
    el("site-select").addEventListener("change", function (event) {
      chooseSite(event.target.value);
    });
    [
      ["h-slider", "h"],
      ["gamma-slider", "gamma"],
      ["delta-slider", "delta"],
      ["antibiotic-slider", "antibioticDays"],
      ["contact-slider", "contactRate"],
      ["susceptible-slider", "susceptibleFraction"]
    ].forEach(function (item) {
      el(item[0]).addEventListener("input", function (event) {
        state[item[1]] = parseFloat(event.target.value);
        updateReadouts();
        scheduleCalculation();
      });
    });
  }

  function clamped(value, range) {
    return Math.max(range.min, Math.min(range.max, value));
  }

  function chooseSite(id) {
    const site = state.config.sites.find(function (candidate) {
      return candidate.id === id;
    }) || state.config.sites[0];
    const defaults = state.config.defaults;
    const ranges = state.config.ranges;
    state.site = site;
    state.h = clamped(defaults.h * site.h_multiplier, ranges.h);
    state.gamma = clamped(defaults.gamma * site.gamma_multiplier, ranges.gamma);
    state.delta = clamped(defaults.delta * site.delta_multiplier, ranges.delta);
    state.antibioticDays = defaults.antibiotic_days;
    state.contactRate = site.contact_rate;
    state.susceptibleFraction = site.susceptible_fraction;
    el("site-select").value = site.id;
    el("site-modalities").textContent = site.modalities;
    updateReadouts();
    updateInterventionChoice();
    scheduleCalculation();
  }

  function updateReadouts() {
    el("h-slider").value = state.h;
    el("gamma-slider").value = state.gamma;
    el("delta-slider").value = state.delta;
    el("antibiotic-slider").value = state.antibioticDays;
    el("contact-slider").value = state.contactRate;
    el("susceptible-slider").value = state.susceptibleFraction;
    el("h-readout").textContent = state.h.toFixed(3);
    el("gamma-readout").textContent = state.gamma.toFixed(3);
    el("delta-readout").textContent = state.delta.toFixed(3);
    el("antibiotic-readout").textContent = state.antibioticDays.toFixed(0) + " days";
    el("contact-readout").textContent = state.contactRate.toFixed(2);
    el("susceptible-readout").textContent =
      (state.susceptibleFraction * 100).toFixed(0) + "%";
  }

  function interventionById(id) {
    return state.config.interventions.find(function (intervention) {
      return intervention.id === id;
    });
  }

  function updateInterventionChoice() {
    document.querySelectorAll("#intervention-choice button").forEach(function (button) {
      button.setAttribute(
        "aria-checked",
        button.dataset.intervention === state.intervention ? "true" : "false"
      );
    });
    const intervention = interventionById(state.intervention);
    el("intervention-note").textContent = intervention.description;
    el("selected-scenario").textContent = intervention.label;
  }

  function parametersFor(intervention) {
    return {
      h: state.h * (intervention === "conjugation_inhibition" ? 0.5 : 1),
      gamma: state.gamma,
      delta: state.delta
    };
  }

  function conditionsFor(intervention) {
    return {
      antibioticDays: intervention === "shorter_antibiotic"
        ? Math.min(3, state.antibioticDays)
        : state.antibioticDays,
      contactRate: state.contactRate,
      susceptibleFraction: state.susceptibleFraction
    };
  }

  function scheduleCalculation() {
    if (state.timer) window.clearTimeout(state.timer);
    el("computing-status").textContent = "Calculating…";
    el("mechanism-results").setAttribute("aria-busy", "true");
    state.timer = window.setTimeout(calculate, 45);
  }

  function calculate() {
    try {
      const baseline = AmplifyMultiscale.targets(
        parametersFor("baseline"),
        conditionsFor("baseline"),
        state.config
      );
      const selected = state.intervention === "baseline"
        ? baseline
        : AmplifyMultiscale.targets(
          parametersFor(state.intervention),
          conditionsFor(state.intervention),
          state.config
        );
      renderMetrics(selected, baseline);
      renderMatrix(selected.matrix);
      el("product-comparator").textContent = format(selected.productBetween);
      el("carriage-duration").textContent = selected.carriageDuration.toFixed(1);
      el("computing-status").textContent = "Updated";
    } catch (error) {
      el("computing-status").textContent = "Calculation failed";
      console.error(error);
    }
    el("mechanism-results").setAttribute("aria-busy", "false");
  }

  function format(value) {
    if (value >= 10) return value.toFixed(1);
    return value.toFixed(2);
  }

  function renderMetrics(selected, baseline) {
    ["r0Within", "reWithin", "r0Between", "reBetween"].forEach(function (name) {
      const value = selected[name];
      const card = document.querySelector("[data-metric=" + name + "]").closest(".metric-card");
      document.querySelector("[data-metric=" + name + "]").textContent = format(value);
      const deltaNode = document.querySelector("[data-delta=" + name + "]");
      card.classList.remove("is-growing", "is-declining");
      card.classList.add(value > 1 ? "is-growing" : "is-declining");
      if (state.intervention === "baseline") {
        deltaNode.textContent = value > 1 ? "above replacement" : "below replacement";
      } else {
        const delta = value - baseline[name];
        const sign = delta > 0 ? "+" : delta < 0 ? "−" : "";
        deltaNode.textContent = sign + Math.abs(delta).toFixed(2) + " vs current conditions";
      }
    });
    const within = selected.reWithin > 1 ? "can spread to additional backgrounds" :
      "is not expected to replace itself across backgrounds";
    const between = selected.reBetween > 1 ? "can sustain secondary colonisations" :
      "is not expected to replace one carrier between hosts";
    el("threshold-note").innerHTML =
      "At the selected conditions, the pARG <strong>" + within +
      "</strong> within a host, and one carrier <strong>" + between +
      "</strong>. These thresholds refer to different infectious units.";
  }

  function renderMatrix(matrix) {
    const backgrounds = state.config.backgrounds;
    let html = "<thead><tr><th>Recipient ↓ / donor →</th>";
    backgrounds.forEach(function (background) {
      html += "<th>" + background.id + "</th>";
    });
    html += "</tr></thead><tbody>";
    matrix.forEach(function (row, recipient) {
      html += "<tr><th>" + backgrounds[recipient].id + "<small>" +
        backgrounds[recipient].species + "</small></th>";
      row.forEach(function (value, donor) {
        html += recipient === donor
          ? "<td class=\"matrix-zero\">—</td>"
          : "<td>" + value.toFixed(3) + "</td>";
      });
      html += "</tr>";
    });
    el("matrix-table").innerHTML = html + "</tbody>";
  }

  /* ---------- fitted observation-to-R surrogate ---------- */

  function initializeSurrogate() {
    const center = state.surrogate.preprocessor.center;
    const inputFields = {
      "ml-qpcr-baseline": "qpcr_baseline",
      "ml-qpcr-peak": "qpcr_peak",
      "ml-qpcr-day30": "qpcr_day30",
      "ml-ecoli-day30": "ecoli_day30",
      "ml-klebsiella-day30": "klebsiella_day30",
      "ml-linked-backgrounds": "linked_backgrounds",
      "ml-linkage-observations": "linkage_observations"
    };
    Object.keys(inputFields).forEach(function (id) {
      const input = el(id);
      const value = center[inputFields[id]];
      input.value = Number.isFinite(value) ? value.toFixed(input.step === "1" ? 0 : 2) : "0";
      input.addEventListener("input", scheduleSurrogatePrediction);
    });
    ["quantitative", "genomic", "clinical"].forEach(function (block) {
      el("ml-use-" + block).addEventListener("change", function (event) {
        state["use" + block.charAt(0).toUpperCase() + block.slice(1)] = event.target.checked;
        updateSurrogateControls();
        scheduleSurrogatePrediction();
      });
    });
    el("ml-site-select").value = state.config.defaults.site;
    el("ml-site-select").addEventListener("change", function (event) {
      applySurrogateSite(event.target.value);
    });
    el("ml-intervention-select").value = "baseline";
    el("ml-intervention-select").addEventListener("change", scheduleSurrogatePrediction);
    ["ml-antibiotic-days", "ml-contact-rate", "ml-susceptible-fraction"].forEach(
      function (id) { el(id).addEventListener("input", scheduleSurrogatePrediction); }
    );
    const training = state.surrogate.generated.training;
    el("ml-training-profiles").textContent = training.profiles;
    el("ml-training-rows").textContent = training.profile_scenarios;
    el("ml-model-date").textContent = "Weights " + state.surrogate.generated.date;
    applySurrogateSite(state.config.defaults.site);
    updateSurrogateControls();
  }

  function applySurrogateSite(id) {
    const site = state.config.sites.find(function (candidate) {
      return candidate.id === id;
    }) || state.config.sites[0];
    el("ml-site-select").value = site.id;
    el("ml-antibiotic-days").value = state.config.defaults.antibiotic_days;
    el("ml-contact-rate").value = site.contact_rate;
    el("ml-susceptible-fraction").value = site.susceptible_fraction;
    scheduleSurrogatePrediction();
  }

  function updateSurrogateControls() {
    ["quantitative", "genomic", "clinical"].forEach(function (block) {
      const key = "use" + block.charAt(0).toUpperCase() + block.slice(1);
      el("ml-control-" + block).classList.toggle("is-off", !state[key]);
    });
  }

  function numericValue(id) {
    const value = parseFloat(el(id).value);
    return Number.isFinite(value) ? value : NaN;
  }

  function surrogateData() {
    const intervention = el("ml-intervention-select").value;
    return {
      qpcr_baseline: numericValue("ml-qpcr-baseline"),
      qpcr_peak: numericValue("ml-qpcr-peak"),
      qpcr_day30: numericValue("ml-qpcr-day30"),
      ecoli_day30: numericValue("ml-ecoli-day30"),
      klebsiella_day30: numericValue("ml-klebsiella-day30"),
      linked_backgrounds: numericValue("ml-linked-backgrounds"),
      linkage_observations: numericValue("ml-linkage-observations"),
      antibiotic_days: intervention === "shorter_antibiotic"
        ? Math.min(3, numericValue("ml-antibiotic-days"))
        : numericValue("ml-antibiotic-days"),
      contact_rate: numericValue("ml-contact-rate"),
      susceptible_fraction: numericValue("ml-susceptible-fraction"),
      conjugation_multiplier: intervention === "conjugation_inhibition" ? 0.5 : 1,
      site_id: el("ml-site-select").value
    };
  }

  function observedBlocks() {
    const blocks = [];
    if (state.useQuantitative) blocks.push("quantitative");
    if (state.useGenomic) blocks.push("genomic");
    if (state.useClinical) blocks.push("clinical");
    return blocks;
  }

  function scheduleSurrogatePrediction() {
    if (state.mlTimer) window.clearTimeout(state.mlTimer);
    state.mlTimer = window.setTimeout(predictSurrogate, 20);
  }

  function patternLabel(pattern) {
    return {
      all: "All three modalities",
      no_quantitative: "Genomic + clinical/context",
      no_genomic: "Quantitative + clinical/context",
      no_clinical: "Quantitative + genomic",
      quantitative_only: "Quantitative only",
      genomic_only: "Genomic only",
      clinical_only: "Clinical/context only"
    }[pattern] || "No modalities";
  }

  function predictSurrogate() {
    const blocks = observedBlocks();
    const error = el("ml-form-error");
    if (!blocks.length) {
      error.hidden = false;
      document.querySelectorAll("[data-surrogate-metric]").forEach(function (node) {
        node.textContent = "—";
        node.closest(".metric-card").classList.remove("is-growing", "is-declining");
      });
      el("ml-pattern-label").textContent = "Waiting for observations";
      el("ml-threshold-note").textContent = "";
      el("ml-heldout-rmse").textContent = "—";
      return;
    }
    error.hidden = true;
    const prediction = AmplifyMultiscaleSurrogate.predict(
      state.surrogate,
      surrogateData(),
      blocks
    );
    state.surrogate.target_names.forEach(function (target) {
      const value = prediction[target];
      const node = document.querySelector("[data-surrogate-metric=" + target + "]");
      const card = node.closest(".metric-card");
      node.textContent = format(value);
      card.classList.remove("is-growing", "is-declining");
      card.classList.add(value > 1 ? "is-growing" : "is-declining");
      document.querySelector("[data-surrogate-status=" + target + "]").textContent =
        value > 1 ? "above replacement" : "below replacement";
    });
    el("ml-pattern-label").textContent = patternLabel(prediction.pattern);
    const within = prediction.re_within > 1 ? "above" : "below";
    const between = prediction.re_between > 1 ? "above" : "below";
    el("ml-threshold-note").innerHTML =
      "Estimated effective R is <strong>" + within +
      " replacement within host</strong> and <strong>" + between +
      " replacement between hosts</strong>. The units are distinct.";
    const evaluation = AmplifyMultiscaleSurrogate.evaluationFor(
      state.surrogate,
      prediction.pattern
    );
    el("ml-heldout-rmse").textContent = evaluation.map(function (row) {
      const label = {
        r0_within: "R₀ WH",
        re_within: "Rₑ WH",
        r0_between: "R₀ BH",
        re_between: "Rₑ BH"
      }[row.target];
      return label + " " + Number(row.rmse).toFixed(2);
    }).join(" · ");
  }
})();
