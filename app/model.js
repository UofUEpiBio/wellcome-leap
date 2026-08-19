/*
 * Dependency-free inference for the AMPLIFY masked count surrogate.
 *
 * Mirrors the evaluation-mode forward pass of the R `torch` module
 * `masked_count_net` defined in R/torch_model.R. The weights are read from
 * model.json, which is written by scripts/05_export_web_model.R.
 */
const AmplifyModel = (function () {
  "use strict";

  const PATTERNS = {
    both: "both",
    xOnly: "x_only",
    scenarioOnly: "scenario_only"
  };

  function relu(value) {
    return value > 0 ? value : 0;
  }

  function softplus(value) {
    // log(1 + exp(z)), evaluated in the numerically stable form.
    return Math.log1p(Math.exp(-Math.abs(value))) + Math.max(value, 0);
  }

  /** Apply one dense layer stored in row-major order. */
  function dense(values, layer, activation) {
    const output = new Array(layer.units);
    for (let unit = 0; unit < layer.units; unit++) {
      let total = layer.bias[unit];
      const offset = unit * layer.inputs;
      for (let input = 0; input < layer.inputs; input++) {
        total += layer.weight[offset + input] * values[input];
      }
      output[unit] = activation(total);
    }
    return output;
  }

  /** Build the four-column model input from the user-supplied modalities. */
  function buildInput(model, x, scenario) {
    const hasX = typeof x === "number" && Number.isFinite(x);
    const hasScenario =
      typeof scenario === "string" &&
      Object.prototype.hasOwnProperty.call(
        model.preprocessor.scenario_levels,
        scenario
      );
    if (!hasX && !hasScenario) {
      throw new Error("Provide the individual feature, the scenario, or both.");
    }
    const standardized = hasX
      ? (x - model.preprocessor.x_mean) / model.preprocessor.x_sd
      : 0;
    const encoded = hasScenario
      ? model.preprocessor.scenario_levels[scenario]
      : 0;
    return {
      values: [
        standardized,
        encoded,
        hasX ? 1 : 0,
        hasScenario ? 1 : 0
      ],
      pattern: hasX && hasScenario
        ? PATTERNS.both
        : hasX
          ? PATTERNS.xOnly
          : PATTERNS.scenarioOnly
    };
  }

  /**
   * Predict expected secondary infections for one individual.
   *
   * @param {object} model Parsed model.json export.
   * @param {number|null} x Individual feature, or null when unobserved.
   * @param {string|null} scenario "lower", "higher", or null when unobserved.
   * @returns {{value: number, pattern: string}}
   */
  function predict(model, x, scenario) {
    const input = buildInput(model, x, scenario);
    let values = dense(input.values, model.layers.hidden_1, relu);
    values = dense(values, model.layers.hidden_2, relu);
    values = dense(values, model.layers.output, softplus);
    return {
      value: values[0] + model.architecture.output_offset,
      pattern: input.pattern
    };
  }

  /** Predict across a sequence of feature values holding the scenario fixed. */
  function predictCurve(model, xValues, scenario) {
    return xValues.map(function (x) {
      return predict(model, x, scenario).value;
    });
  }

  return {
    PATTERNS: PATTERNS,
    predict: predict,
    predictCurve: predictCurve
  };
})();
