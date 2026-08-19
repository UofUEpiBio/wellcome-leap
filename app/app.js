/*
 * AMPLIFY prototype interface.
 *
 * Loads the exported weights (model.json) and the editable presentation
 * settings (site.json), then re-evaluates the surrogate on every input change.
 * Everything runs locally; there is no server component.
 */
(function () {
  "use strict";

  const state = {
    model: null,
    site: null,
    useX: true,
    useScenario: true,
    x: 0,
    scenario: "higher"
  };

  const el = function (id) { return document.getElementById(id); };

  /* ---------- loading ---------- */

  function loadJson(path) {
    return fetch(path, { cache: "no-cache" }).then(function (response) {
      if (!response.ok) throw new Error(path + ": " + response.status);
      return response.json();
    });
  }

  Promise.all([loadJson("model.json"), loadJson("site.json")])
    .then(function (loaded) {
      state.model = loaded[0];
      state.site = loaded[1];
      applySite();
      buildScenarioChoice();
      bindControls();
      update();
    })
    .catch(function (error) {
      el("load-error").hidden = false;
      el("estimator").hidden = true;
      console.error(error);
    });

  /* ---------- presentation settings ---------- */

  function applySite() {
    const site = state.site;
    document.querySelectorAll("[data-site]").forEach(function (node) {
      const value = site[node.dataset.site];
      if (typeof value === "string" && value.length) node.textContent = value;
    });
    document.querySelectorAll("[data-threshold]").forEach(function (node) {
      node.textContent = formatNumber(site.thresholds[node.dataset.threshold], 2);
    });
    document.title = site.consortium + " — " + site.app_name;
    if (site.repository_url) el("repo-link").href = site.repository_url;
    if (Array.isArray(site.partner_institutions)) {
      el("footer-partners").textContent = site.partner_institutions.join(" · ");
    }
    const generated = state.model.generated;
    if (generated && generated.date) el("model-date").textContent = generated.date;
  }

  function scenarios() {
    return state.site.scenarios;
  }

  function scenarioById(id) {
    return scenarios().find(function (item) { return item.id === id; }) || null;
  }

  /* ---------- controls ---------- */

  function buildScenarioChoice() {
    const container = el("scenario-choice");
    container.innerHTML = "";
    scenarios().forEach(function (scenario) {
      const button = document.createElement("button");
      button.type = "button";
      button.role = "radio";
      button.textContent = scenario.label;
      button.dataset.scenario = scenario.id;
      button.addEventListener("click", function () {
        state.scenario = scenario.id;
        update();
      });
      container.appendChild(button);
    });
  }

  function bindControls() {
    el("use-x").addEventListener("change", function (event) {
      state.useX = event.target.checked;
      update();
    });
    el("use-scenario").addEventListener("change", function (event) {
      state.useScenario = event.target.checked;
      update();
    });
    el("x-slider").addEventListener("input", function (event) {
      state.x = parseFloat(event.target.value);
      update();
    });
    // The chart scales with its viewBox, so resizing needs no redraw.
  }

  /* ---------- formatting ---------- */

  function formatNumber(value, digits) {
    return value.toFixed(typeof digits === "number" ? digits : 2);
  }

  function band(value) {
    const thresholds = state.site.thresholds;
    if (value > thresholds.outbreak) {
      return {
        key: "is-outbreak",
        label: "Outbreak potential",
        meaning:
          "Above one expected secondary infection: this individual is expected to more than " +
          "replace themselves, so transmission can grow."
      };
    }
    if (value < thresholds.self_limiting) {
      return {
        key: "is-safe",
        label: "Self-limiting",
        meaning:
          "Below one expected secondary infection: chains started by this individual are " +
          "expected to die out."
      };
    }
    return {
      key: "is-borderline",
      label: "Borderline",
      meaning:
        "Close enough to one that the model cannot separate growth from decline for this " +
        "individual. Treat it as unresolved rather than safe."
    };
  }

  function patternLabel(pattern) {
    if (pattern === "both") return "Feature + scenario";
    if (pattern === "x_only") return "Feature only";
    return "Scenario only";
  }

  /* ---------- update cycle ---------- */

  function update() {
    el("control-x").classList.toggle("is-off", !state.useX);
    el("control-scenario").classList.toggle("is-off", !state.useScenario);
    el("x-readout").textContent = formatNumber(state.x);
    el("x-slider").value = String(state.x);

    document.querySelectorAll("#scenario-choice button").forEach(function (button) {
      const active = button.dataset.scenario === state.scenario;
      button.setAttribute("aria-checked", active ? "true" : "false");
    });
    const active = scenarioById(state.scenario);
    el("scenario-note").textContent = active ? active.note : "";

    const result = el("result");
    result.classList.remove("is-outbreak", "is-borderline", "is-safe");

    if (!state.useX && !state.useScenario) {
      el("form-error").hidden = false;
      el("result-value").textContent = "—";
      el("result-flag").textContent = "Waiting for input";
      el("result-meaning").textContent = "";
      el("meta-pattern").textContent = "—";
      el("meta-rmse").textContent = "—";
      renderChart();
      return;
    }
    el("form-error").hidden = true;

    const x = state.useX ? state.x : null;
    const scenario = state.useScenario ? state.scenario : null;
    const prediction = AmplifyModel.predict(state.model, x, scenario);
    const flag = band(prediction.value);

    result.classList.add(flag.key);
    el("result-value").textContent = formatNumber(prediction.value);
    el("result-flag").textContent = flag.label;
    el("result-meaning").textContent = flag.meaning;
    el("meta-pattern").textContent = patternLabel(prediction.pattern);

    const rmse = AmplifyModel.testRmse(state.model, prediction.pattern, scenario);
    el("meta-rmse").textContent = rmse === null ? "—" : "± " + formatNumber(rmse);

    renderChart(prediction, flag);
  }

  /* ---------- chart ---------- */

  const CHART = {
    width: 640,
    height: 230,
    left: 42,
    right: 12,
    top: 14,
    bottom: 30
  };

  function xGrid() {
    const values = [];
    for (let x = -3; x <= 3.0001; x += 0.05) values.push(Math.round(x * 100) / 100);
    return values;
  }

  function renderChart(prediction, flag) {
    const host = el("chart");
    if (!state.model) return;

    if (!state.useX && !state.useScenario) {
      host.innerHTML = svgShell(
        '<text class="empty-text" x="' + (CHART.width / 2) +
        '" y="' + (CHART.height / 2) + '" text-anchor="middle">' +
        "Provide at least one attribute to see the response.</text>"
      );
      el("chart-caption").textContent =
        "Predicted individual reproduction number across the feature range.";
      return;
    }

    const scenario = state.useScenario ? state.scenario : null;
    const grid = xGrid();
    const constant = !state.useX;
    const values = constant
      ? grid.map(function () { return prediction.value; })
      : AmplifyModel.predictCurve(state.model, grid, scenario);

    const maxValue = Math.max(1.35, Math.max.apply(null, values) * 1.12);
    const plotWidth = CHART.width - CHART.left - CHART.right;
    const plotHeight = CHART.height - CHART.top - CHART.bottom;
    const px = function (x) { return CHART.left + ((x + 3) / 6) * plotWidth; };
    const py = function (v) { return CHART.top + plotHeight - (v / maxValue) * plotHeight; };

    let parts = "";

    // y grid and labels
    const ticks = niceTicks(maxValue);
    ticks.forEach(function (tick) {
      const y = py(tick);
      parts += '<line class="grid-line" x1="' + CHART.left + '" y1="' + y +
        '" x2="' + (CHART.width - CHART.right) + '" y2="' + y + '"/>';
      parts += '<text class="axis-text" x="' + (CHART.left - 8) + '" y="' + (y + 3.5) +
        '" text-anchor="end">' + tick + "</text>";
    });

    // x labels
    [-3, -2, -1, 0, 1, 2, 3].forEach(function (tick) {
      parts += '<text class="axis-text" x="' + px(tick) + '" y="' +
        (CHART.height - CHART.bottom + 16) + '" text-anchor="middle">' +
        (tick > 0 ? "+" + tick : tick) + "</text>";
    });
    parts += '<text class="axis-text" x="' + (CHART.left + plotWidth / 2) + '" y="' +
      (CHART.height - 2) + '" text-anchor="middle">Individual feature (standardized)</text>';

    // curve
    const line = values.map(function (value, index) {
      return (index ? "L" : "M") + px(grid[index]).toFixed(2) + " " + py(value).toFixed(2);
    }).join(" ");
    const area = line + " L" + px(3).toFixed(2) + " " + py(0).toFixed(2) +
      " L" + px(-3).toFixed(2) + " " + py(0).toFixed(2) + " Z";
    parts += '<path class="curve-area" d="' + area + '"/>';
    parts += '<path class="curve" d="' + line + '"' +
      (constant ? ' stroke-dasharray="6 5"' : "") + "/>";

    // R = 1 reference
    const oneY = py(1);
    parts += '<line class="threshold" x1="' + CHART.left + '" y1="' + oneY +
      '" x2="' + (CHART.width - CHART.right) + '" y2="' + oneY + '"/>';
    parts += '<text class="threshold-text" x="' + (CHART.width - CHART.right) +
      '" y="' + (oneY - 5) + '" text-anchor="end">R = 1</text>';

    // current position
    const markerX = constant ? CHART.left + plotWidth / 2 : px(state.x);
    const markerY = py(prediction.value);
    if (!constant) {
      parts += '<line class="marker-line" x1="' + markerX + '" y1="' + markerY +
        '" x2="' + markerX + '" y2="' + (CHART.top + plotHeight) + '"/>';
    }
    parts += '<circle class="marker-dot ' + flag.key + '" cx="' + markerX +
      '" cy="' + markerY + '" r="5.5"/>';

    host.innerHTML = svgShell(parts);

    const scenarioLabel = scenario
      ? scenarioById(scenario).label.toLowerCase()
      : "no scenario supplied";
    el("chart-caption").textContent = constant
      ? "Scenario only: the estimate does not vary with the individual feature (" +
        scenarioLabel + ")."
      : "Predicted individual reproduction number across the feature range (" +
        scenarioLabel + ").";
  }

  function svgShell(inner) {
    return '<svg viewBox="0 0 ' + CHART.width + " " + CHART.height +
      '" role="img" aria-label="Predicted individual reproduction number">' +
      inner + "</svg>";
  }

  /** Choose at most five readable y-axis ticks below the plotted maximum. */
  function niceTicks(maxValue) {
    const steps = [0.25, 0.5, 1, 2, 2.5, 5, 10];
    const step = steps.find(function (candidate) {
      return maxValue / candidate <= 5;
    }) || 20;
    const ticks = [];
    for (let value = 0; value <= maxValue; value += step) {
      ticks.push(Math.round(value * 100) / 100);
    }
    return ticks;
  }
})();
