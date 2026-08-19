/*
 * AMPLIFY prototype interface.
 *
 * Loads the exported weights (model.json) and the editable presentation
 * settings (site.json), then re-evaluates the surrogate on every input change.
 *
 * The surrogate takes one standardized individual feature and the organism
 * label. The page presents that feature as a clinically recognisable
 * transmissibility index: site.json describes a continuous driver plus a set of
 * categorical modifiers, and this file collapses them into the single value the
 * model consumes. The mapping is presentation only; nothing here touches the
 * fitted weights. Everything runs locally, with no server component.
 */
(function () {
  "use strict";

  /** The training features are standardized, so the model never saw |x| > 3. */
  const INDEX_LIMIT = 3;

  const state = {
    model: null,
    site: null,
    useProfile: true,
    useScenario: true,
    driver: 0,
    modifiers: {},
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
      buildProfileControls();
      buildScenarioChoice();
      buildMappingTable();
      bindControls();
      bindTabs();
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

    el("profile-legend").textContent = profile().label;
    el("profile-note").textContent = profile().note;
    el("index-label").textContent = profile().index_label;
    el("index-note").textContent = profile().index_note;
    el("organism-legend").textContent = site.organism_label;
    el("organism-note").textContent = site.organism_note;
  }

  function profile() {
    return state.site.profile;
  }

  function driver() {
    return profile().driver;
  }

  function modifiers() {
    return profile().modifiers;
  }

  function scenarios() {
    return state.site.scenarios;
  }

  function scenarioById(id) {
    return scenarios().find(function (item) { return item.id === id; }) || null;
  }

  /* ---------- controls ---------- */

  function buildProfileControls() {
    const spec = driver();
    state.driver = spec.value;

    el("driver-label").textContent = spec.label + " (" + spec.unit + ")";
    el("driver-modality").textContent = spec.modality;
    const slider = el("driver-slider");
    slider.min = String(spec.min);
    slider.max = String(spec.max);
    slider.step = String(spec.step);
    slider.value = String(spec.value);

    const hints = el("driver-hints");
    hints.innerHTML = "";
    (spec.hints || []).forEach(function (hint) {
      const span = document.createElement("span");
      span.textContent = hint;
      hints.appendChild(span);
    });

    const host = el("modifiers");
    host.innerHTML = "";
    modifiers().forEach(function (modifier) {
      state.modifiers[modifier.id] = modifier.value;

      const field = document.createElement("div");
      field.className = "field";

      const head = document.createElement("div");
      head.className = "field-head";
      const label = document.createElement("label");
      label.className = "field-label";
      label.htmlFor = "modifier-" + modifier.id;
      label.textContent = modifier.label;
      const chip = document.createElement("span");
      chip.className = "chip";
      chip.textContent = modifier.modality;
      head.appendChild(label);
      head.appendChild(chip);

      const select = document.createElement("select");
      select.className = "select";
      select.id = "modifier-" + modifier.id;
      modifier.options.forEach(function (option, index) {
        const node = document.createElement("option");
        node.value = String(index);
        node.textContent = option.label;
        select.appendChild(node);
      });
      select.value = String(modifier.value);
      select.addEventListener("change", function (event) {
        state.modifiers[modifier.id] = parseInt(event.target.value, 10);
        update();
      });

      field.appendChild(head);
      field.appendChild(select);
      host.appendChild(field);
    });
  }

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

  /** Document the presentation-only index mapping on the methods tab. */
  function buildMappingTable() {
    const spec = driver();
    const body = el("mapping-table").querySelector("tbody");
    body.innerHTML = "";

    const rows = [{
      label: spec.label + " (" + spec.unit + ")",
      modality: spec.modality,
      contribution:
        formatNumber(spec.weight) + " per unit above " + formatNumber(spec.center, 1) +
        " (range " + formatNumber(spec.min, 1) + "–" + formatNumber(spec.max, 1) + ")"
    }];
    modifiers().forEach(function (modifier) {
      rows.push({
        label: modifier.label,
        modality: modifier.modality,
        contribution: modifier.options.map(function (option) {
          return option.label + " " + signed(option.offset);
        }).join("; ")
      });
    });

    rows.forEach(function (row) {
      const tr = document.createElement("tr");
      const name = document.createElement("th");
      name.scope = "row";
      name.textContent = row.label;
      const modality = document.createElement("td");
      modality.textContent = row.modality;
      const contribution = document.createElement("td");
      contribution.textContent = row.contribution;
      tr.appendChild(name);
      tr.appendChild(modality);
      tr.appendChild(contribution);
      body.appendChild(tr);
    });
  }

  function bindControls() {
    el("use-profile").addEventListener("change", function (event) {
      state.useProfile = event.target.checked;
      update();
    });
    el("use-scenario").addEventListener("change", function (event) {
      state.useScenario = event.target.checked;
      update();
    });
    el("driver-slider").addEventListener("input", function (event) {
      state.driver = parseFloat(event.target.value);
      update();
    });
    // The chart scales with its viewBox, so resizing needs no redraw.
  }

  /* ---------- tabs ---------- */

  function bindTabs() {
    document.querySelectorAll("[role=tab]").forEach(function (tab) {
      tab.addEventListener("click", function () { selectTab(tab.dataset.tab, true); });
    });
    document.querySelectorAll("[data-tab-link]").forEach(function (link) {
      link.addEventListener("click", function (event) {
        event.preventDefault();
        selectTab(link.dataset.tabLink, true);
      });
    });
    window.addEventListener("hashchange", function () {
      selectTab(currentHash(), true);
    });
    // On load, leave the scroll position to the browser.
    selectTab(currentHash(), false);
  }

  function currentHash() {
    return window.location.hash.replace("#", "");
  }

  /** Show one tab panel, falling back to the estimator for an unknown name. */
  function selectTab(name, scrollToTop) {
    const tabs = Array.prototype.slice.call(document.querySelectorAll("[role=tab]"));
    const known = tabs.some(function (tab) { return tab.dataset.tab === name; });
    const target = known ? name : "estimator";
    tabs.forEach(function (tab) {
      const active = tab.dataset.tab === target;
      tab.setAttribute("aria-selected", active ? "true" : "false");
      el("panel-" + tab.dataset.tab).hidden = !active;
    });
    if (known && window.location.hash !== "#" + target) {
      history.replaceState(null, "", "#" + target);
    }
    if (scrollToTop) window.scrollTo({ top: 0, behavior: "auto" });
  }

  /* ---------- the presentation-only index ---------- */

  function clampIndex(value) {
    return Math.max(-INDEX_LIMIT, Math.min(INDEX_LIMIT, value));
  }

  /** Offsets contributed by the categorical modifiers at their current values. */
  function modifierOffset() {
    return modifiers().reduce(function (total, modifier) {
      const option = modifier.options[state.modifiers[modifier.id]];
      return total + (option ? option.offset : 0);
    }, 0);
  }

  /** Collapse the carrier profile into the standardized feature the model takes. */
  function indexFor(driverValue) {
    const spec = driver();
    return clampIndex(spec.weight * (driverValue - spec.center) + modifierOffset());
  }

  /* ---------- formatting ---------- */

  function formatNumber(value, digits) {
    return value.toFixed(typeof digits === "number" ? digits : 2);
  }

  /** Signed to make the direction of an index contribution obvious. */
  function signed(value) {
    if (value === 0) return formatNumber(0);
    return (value > 0 ? "+" : "−") + formatNumber(Math.abs(value));
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
    if (pattern === "both") return "Profile + organism";
    if (pattern === "x_only") return "Profile only";
    return "Organism only";
  }

  /* ---------- update cycle ---------- */

  function update() {
    el("control-profile").classList.toggle("is-off", !state.useProfile);
    el("control-scenario").classList.toggle("is-off", !state.useScenario);

    el("driver-readout").textContent = formatNumber(state.driver, 1);
    el("driver-slider").value = String(state.driver);

    const index = indexFor(state.driver);
    el("index-value").textContent = signed(index);

    document.querySelectorAll("#scenario-choice button").forEach(function (button) {
      const active = button.dataset.scenario === state.scenario;
      button.setAttribute("aria-checked", active ? "true" : "false");
    });
    const active = scenarioById(state.scenario);
    el("scenario-note").textContent = active ? active.note : "";

    const result = el("result");
    result.classList.remove("is-outbreak", "is-borderline", "is-safe");

    if (!state.useProfile && !state.useScenario) {
      el("form-error").hidden = false;
      el("result-value").textContent = "—";
      el("result-flag").textContent = "Waiting for input";
      el("result-meaning").textContent = "";
      el("meta-pattern").textContent = "—";
      el("meta-index").textContent = "—";
      renderChart();
      return;
    }
    el("form-error").hidden = true;

    const scenario = state.useScenario ? state.scenario : null;
    const prediction = AmplifyModel.predict(
      state.model,
      state.useProfile ? index : null,
      scenario
    );
    const flag = band(prediction.value);

    result.classList.add(flag.key);
    el("result-value").textContent = formatNumber(prediction.value);
    el("result-flag").textContent = flag.label;
    el("result-meaning").textContent = flag.meaning;
    el("meta-pattern").textContent = patternLabel(prediction.pattern);
    el("meta-index").textContent = state.useProfile ? signed(index) : "not supplied";

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

  /** Values of the continuous driver at which to evaluate the response curve. */
  function driverGrid() {
    const spec = driver();
    const steps = 120;
    const values = [];
    for (let step = 0; step <= steps; step++) {
      values.push(spec.min + ((spec.max - spec.min) * step) / steps);
    }
    return values;
  }

  function renderChart(prediction, flag) {
    const host = el("chart");
    if (!state.model) return;

    const spec = driver();
    const axisLabel = spec.label + " (" + spec.unit + ")";

    if (!state.useProfile && !state.useScenario) {
      host.innerHTML = svgShell(
        '<text class="empty-text" x="' + (CHART.width / 2) +
        '" y="' + (CHART.height / 2) + '" text-anchor="middle">' +
        "Provide the carrier profile, the organism, or both to see the response.</text>"
      );
      el("chart-caption").textContent =
        "Predicted individual reproduction number across the " + spec.label.toLowerCase() + " range.";
      return;
    }

    const scenario = state.useScenario ? state.scenario : null;
    const grid = driverGrid();
    const constant = !state.useProfile;
    const values = constant
      ? grid.map(function () { return prediction.value; })
      : AmplifyModel.predictCurve(state.model, grid.map(indexFor), scenario);

    const maxValue = Math.max(1.35, Math.max.apply(null, values) * 1.12);
    const plotWidth = CHART.width - CHART.left - CHART.right;
    const plotHeight = CHART.height - CHART.top - CHART.bottom;
    const px = function (value) {
      return CHART.left + ((value - spec.min) / (spec.max - spec.min)) * plotWidth;
    };
    const py = function (value) {
      return CHART.top + plotHeight - (value / maxValue) * plotHeight;
    };

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
    axisTicks(spec.min, spec.max).forEach(function (tick) {
      parts += '<text class="axis-text" x="' + px(tick) + '" y="' +
        (CHART.height - CHART.bottom + 16) + '" text-anchor="middle">' + tick + "</text>";
    });
    parts += '<text class="axis-text" x="' + (CHART.left + plotWidth / 2) + '" y="' +
      (CHART.height - 2) + '" text-anchor="middle">' + axisLabel + "</text>";

    // curve
    const line = values.map(function (value, index) {
      return (index ? "L" : "M") + px(grid[index]).toFixed(2) + " " + py(value).toFixed(2);
    }).join(" ");
    const area = line + " L" + px(spec.max).toFixed(2) + " " + py(0).toFixed(2) +
      " L" + px(spec.min).toFixed(2) + " " + py(0).toFixed(2) + " Z";
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
    const markerX = constant ? CHART.left + plotWidth / 2 : px(state.driver);
    const markerY = py(prediction.value);
    if (!constant) {
      parts += '<line class="marker-line" x1="' + markerX + '" y1="' + markerY +
        '" x2="' + markerX + '" y2="' + (CHART.top + plotHeight) + '"/>';
    }
    parts += '<circle class="marker-dot ' + flag.key + '" cx="' + markerX +
      '" cy="' + markerY + '" r="5.5"/>';

    host.innerHTML = svgShell(parts);

    const scenarioLabel = scenario ? scenarioById(scenario).label : "no organism supplied";
    el("chart-caption").textContent = constant
      ? "Organism only: the estimate does not vary with " + spec.label.toLowerCase() +
        " (" + scenarioLabel + ")."
      : "Predicted individual reproduction number across the " + spec.label.toLowerCase() +
        " range, holding the other profile inputs fixed (" + scenarioLabel + ").";
  }

  function svgShell(inner) {
    return '<svg viewBox="0 0 ' + CHART.width + " " + CHART.height +
      '" role="img" aria-label="Predicted individual reproduction number">' +
      inner + "</svg>";
  }

  /** Whole-number ticks spanning the driver range, at most eight of them. */
  function axisTicks(min, max) {
    const step = Math.max(1, Math.ceil((max - min) / 8));
    const ticks = [];
    for (let value = Math.ceil(min); value <= max; value += step) ticks.push(value);
    return ticks;
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
