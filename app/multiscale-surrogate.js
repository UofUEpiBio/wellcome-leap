/* Dependency-free inference for the exported multiscale R torch emulator. */
const AmplifyMultiscaleSurrogate = (function () {
  "use strict";

  function relu(value) {
    return value > 0 ? value : 0;
  }

  function dense(values, layer, activation) {
    const output = new Array(layer.units);
    for (let unit = 0; unit < layer.units; unit += 1) {
      let total = layer.bias[unit];
      const offset = unit * layer.inputs;
      for (let input = 0; input < layer.inputs; input += 1) {
        total += layer.weight[offset + input] * values[input];
      }
      output[unit] = activation ? activation(total) : total;
    }
    return output;
  }

  function includes(values, target) {
    return values.indexOf(target) >= 0;
  }

  function patternFor(model, observedBlocks) {
    const names = Object.keys(model.patterns);
    return names.find(function (name) {
      const pattern = model.patterns[name];
      return pattern.length === observedBlocks.length && pattern.every(function (block) {
        return includes(observedBlocks, block);
      });
    }) || null;
  }

  function buildInput(model, rawData, observedBlocks) {
    if (!observedBlocks.length) throw new Error("Provide at least one modality.");
    const prepared = Object.assign({}, rawData);
    ["site_a", "site_b", "site_c", "site_d"].forEach(function (site) {
      prepared[site] = rawData.site_id === site ? 1 : 0;
    });
    const values = [];
    Object.keys(model.preprocessor.blocks).forEach(function (block) {
      const observed = includes(observedBlocks, block);
      model.preprocessor.blocks[block].forEach(function (field) {
        const supplied = Number(prepared[field]);
        const raw = Number.isFinite(supplied) ? supplied : model.preprocessor.center[field];
        values.push(observed
          ? (raw - model.preprocessor.center[field]) / model.preprocessor.scale[field]
          : 0);
      });
    });
    Object.keys(model.preprocessor.blocks).forEach(function (block) {
      values.push(includes(observedBlocks, block) ? 1 : 0);
    });
    return { values: values, pattern: patternFor(model, observedBlocks) };
  }

  function predict(model, rawData, observedBlocks) {
    const input = buildInput(model, rawData, observedBlocks);
    let values = dense(input.values, model.layers.hidden_1, relu);
    values = dense(values, model.layers.hidden_2, relu);
    values = dense(values, model.layers.output, null).map(Math.exp);
    const result = { pattern: input.pattern };
    model.target_names.forEach(function (target, index) {
      result[target] = values[index];
    });
    return result;
  }

  function evaluationFor(model, pattern) {
    return (model.evaluation || []).filter(function (row) {
      return row.evaluation_scope === "profile_holdout" && row.pattern === pattern;
    });
  }

  return {
    buildInput: buildInput,
    predict: predict,
    evaluationFor: evaluationFor
  };
})();
