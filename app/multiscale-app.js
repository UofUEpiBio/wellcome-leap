/* Browser controller for the synthetic multiscale mechanism explorer. */
(function () {
  "use strict";

  const state = {
    config: null,
    site: null,
    intervention: "baseline",
    h: 0,
    gamma: 0,
    delta: 0,
    antibioticDays: 0,
    contactRate: 0,
    susceptibleFraction: 0,
    timer: null
  };

  const el = function (id) { return document.getElementById(id); };

  fetch("multiscale.json", { cache: "no-cache" })
    .then(function (response) {
      if (!response.ok) throw new Error("multiscale.json: " + response.status);
      return response.json();
    })
    .then(function (config) {
      state.config = config;
      buildSiteChoice();
      configureSliders();
      buildInterventionChoice();
      bindControls();
      chooseSite(config.defaults.site);
    })
    .catch(function (error) {
      el("computing-status").textContent = "Model unavailable";
      el("mechanism-results").setAttribute("aria-busy", "false");
      console.error(error);
    });

  function buildSiteChoice() {
    const select = el("site-select");
    state.config.sites.forEach(function (site) {
      const option = document.createElement("option");
      option.value = site.id;
      option.textContent = site.label;
      select.appendChild(option);
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
})();
