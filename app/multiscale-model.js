/*
 * Dependency-free browser implementation of the simplified multiscale model.
 * The equations mirror R/within_host.R and R/reproduction.R. Inputs and
 * outputs are synthetic and intended for mechanism exploration, not clinical
 * decision-making.
 */
(function (root, factory) {
  "use strict";
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  root.AmplifyMultiscale = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  function addScaled(left, right, scale) {
    return left.map(function (value, index) {
      return value + scale * right[index];
    });
  }

  function antibioticActivity(time, endTime, config) {
    const start = config.mechanism.antibiotic_start;
    return endTime > start && time >= start && time < endTime ? 1 : 0;
  }

  function initialState(config, resistant, seedSize, donorIndex) {
    const mechanism = config.mechanism;
    const susceptible = config.backgrounds.map(function (background) {
      return mechanism.initial_total * background.initial_share;
    });
    const carrying = config.backgrounds.map(function () { return 0; });
    if (resistant) {
      carrying[donorIndex] = seedSize;
      susceptible[donorIndex] -= seedSize;
    }
    return susceptible.concat(carrying);
  }

  function derivative(time, state, parameters, antibioticEnd, config) {
    const n = config.backgrounds.length;
    const susceptible = state.slice(0, n);
    const resistant = state.slice(n);
    const total = state.reduce(function (sum, value) { return sum + value; }, 0);
    const antibiotic = antibioticActivity(time, antibioticEnd, config);
    const dSusceptible = [];
    const dResistant = [];

    for (let recipient = 0; recipient < n; recipient += 1) {
      const background = config.backgrounds[recipient];
      const birthsS = background.growth_rate * susceptible[recipient];
      const birthsR = background.growth_rate *
        (1 - background.fitness_cost) * resistant[recipient];
      let donorPressure = 0;
      for (let donor = 0; donor < n; donor += 1) {
        donorPressure += config.omega[recipient][donor] * resistant[donor];
      }
      const transfer = susceptible[recipient] / config.mechanism.carrying_capacity *
        parameters.h * donorPressure;
      dSusceptible.push(
        birthsS - birthsS * total / config.mechanism.carrying_capacity -
        parameters.delta * antibiotic * susceptible[recipient] - transfer +
        parameters.gamma * birthsR
      );
      dResistant.push(
        (1 - parameters.gamma) * birthsR -
        birthsR * total / config.mechanism.carrying_capacity + transfer
      );
    }
    return dSusceptible.concat(dResistant);
  }

  function rk4Step(time, state, step, parameters, antibioticEnd, config, endsAtBreak) {
    const endTime = endsAtBreak ? time + step - 1e-10 : time + step;
    const k1 = derivative(time, state, parameters, antibioticEnd, config);
    const k2 = derivative(
      time + step / 2,
      addScaled(state, k1, step / 2),
      parameters,
      antibioticEnd,
      config
    );
    const k3 = derivative(
      time + step / 2,
      addScaled(state, k2, step / 2),
      parameters,
      antibioticEnd,
      config
    );
    const k4 = derivative(
      endTime,
      addScaled(state, k3, step),
      parameters,
      antibioticEnd,
      config
    );
    return state.map(function (value, index) {
      const next = value + step *
        (k1[index] + 2 * k2[index] + 2 * k3[index] + k4[index]) / 6;
      if (!Number.isFinite(next) || next < -1e-8) {
        throw new Error("ODE integration produced an invalid state.");
      }
      return Math.max(0, next);
    });
  }

  function timeGrid(endTime, step, breakpoints) {
    const times = [0];
    for (let time = step; time < endTime - 1e-10; time += step) {
      times.push(Number(time.toFixed(10)));
    }
    (breakpoints || []).forEach(function (time) {
      if (time > 0 && time < endTime) times.push(time);
    });
    times.push(endTime);
    return times.sort(function (left, right) { return left - right; }).filter(
      function (time, index, values) {
        return index === 0 || Math.abs(time - values[index - 1]) > 1e-10;
      }
    );
  }

  function simulate(parameters, antibioticEnd, endTime, step, config, state) {
    const breakpoints = [config.mechanism.antibiotic_start, antibioticEnd];
    const times = timeGrid(endTime, step, breakpoints);
    const states = [state.slice()];
    let current = state.slice();
    for (let index = 0; index < times.length - 1; index += 1) {
      const dt = times[index + 1] - times[index];
      const endsAtBreak = breakpoints.some(function (point) {
        return Math.abs(times[index + 1] - point) < 1e-10;
      });
      current = rk4Step(
        times[index],
        current,
        dt,
        parameters,
        antibioticEnd,
        config,
        endsAtBreak
      );
      states.push(current);
    }
    return { times: times, states: states };
  }

  function trapz(times, values) {
    let total = 0;
    for (let index = 0; index < times.length - 1; index += 1) {
      total += (times[index + 1] - times[index]) *
        (values[index] + values[index + 1]) / 2;
    }
    return total;
  }

  function spectralRadius(matrix) {
    const n = matrix.length;
    let vector = matrix.map(function () { return 1 / n; });
    let value = 0;
    for (let iteration = 0; iteration < 500; iteration += 1) {
      const product = matrix.map(function (row) {
        return row.reduce(function (sum, entry, index) {
          return sum + entry * vector[index];
        }, 0);
      });
      const nextValue = product.reduce(function (sum, entry) { return sum + entry; }, 0);
      if (nextValue === 0) return 0;
      vector = product.map(function (entry) { return entry / nextValue; });
      if (Math.abs(nextValue - value) < 1e-12) return nextValue;
      value = nextValue;
    }
    return value;
  }

  function nextGeneration(parameters, antibioticEnd, config) {
    const n = config.backgrounds.length;
    const mechanism = config.mechanism;
    const susceptibleOnly = initialState(config, false, 0, 0);
    const reference = simulate(
      parameters,
      antibioticEnd,
      mechanism.ngm_horizon,
      mechanism.ngm_step,
      config,
      susceptibleOnly
    );
    const susceptible = reference.states.map(function (state) { return state.slice(0, n); });
    const totalReference = susceptible.map(function (state) {
      return state.reduce(function (sum, value) { return sum + value; }, 0);
    });
    const lambda = config.backgrounds.map(function () {
      return config.backgrounds.map(function () { return 0; });
    });

    for (let donor = 0; donor < n; donor += 1) {
      const background = config.backgrounds[donor];
      const growth = background.growth_rate * (1 - background.fitness_cost);
      const lineage = [mechanism.ngm_seed];
      for (let index = 0; index < reference.times.length - 1; index += 1) {
        const dt = reference.times[index + 1] - reference.times[index];
        const rate1 = growth *
          (1 - parameters.gamma - totalReference[index] / mechanism.carrying_capacity);
        const rate2 = growth * (1 - parameters.gamma -
          (totalReference[index] + totalReference[index + 1]) /
          (2 * mechanism.carrying_capacity));
        const rate4 = growth *
          (1 - parameters.gamma - totalReference[index + 1] / mechanism.carrying_capacity);
        const y = lineage[index];
        const k1 = rate1 * y;
        const k2 = rate2 * (y + dt * k1 / 2);
        const k3 = rate2 * (y + dt * k2 / 2);
        const k4 = rate4 * (y + dt * k3);
        lineage.push(Math.max(0, y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6));
      }
      for (let recipient = 0; recipient < n; recipient += 1) {
        if (recipient === donor) continue;
        const hazard = reference.times.map(function (_, index) {
          return parameters.h * config.omega[recipient][donor] *
            susceptible[index][recipient] * lineage[index] /
            mechanism.carrying_capacity;
        });
        lambda[recipient][donor] = trapz(reference.times, hazard) /
          mechanism.ngm_seed;
      }
    }

    const matrix = lambda.map(function (row, recipient) {
      return row.map(function (exposure, donor) {
        if (recipient === donor) return 0;
        return 1 - Math.exp(-exposure * config.backgrounds[recipient].establishment);
      });
    });
    return {
      exposure: lambda,
      matrix: matrix,
      reproduction: spectralRadius(matrix)
    };
  }

  function betweenHost(trajectory, contactRate, susceptibleFraction, config) {
    const n = config.backgrounds.length;
    const resistant = trajectory.states.map(function (state) {
      return state.slice(n).reduce(function (sum, value) { return sum + value; }, 0);
    });
    const total = trajectory.states.map(function (state) {
      return state.reduce(function (sum, value) { return sum + value; }, 0);
    });
    const fraction = resistant.map(function (value, index) {
      return value / Math.max(total[index], 1e-12);
    });
    let last = -1;
    resistant.forEach(function (value, index) {
      if (value >= config.mechanism.carriage_detection) last = index;
    });
    if (last < 0) {
      return { reproduction: 0, productComparator: 0, duration: 0, integratedRisk: 0 };
    }
    const times = trajectory.times.slice(0, last + 1);
    const fractions = fraction.slice(0, last + 1);
    const risks = fractions.map(function (value) {
      return 1 - Math.exp(-config.mechanism.inoculum_kappa * value);
    });
    const duration = times[times.length - 1];
    const integratedRisk = trapz(times, risks);
    const averageFraction = duration > 0 ? trapz(times, fractions) / duration : 0;
    const multiplier = contactRate * susceptibleFraction *
      config.mechanism.between_establishment;
    return {
      reproduction: multiplier * integratedRisk,
      productComparator: multiplier * duration *
        (1 - Math.exp(-config.mechanism.inoculum_kappa * averageFraction)),
      duration: duration,
      integratedRisk: integratedRisk
    };
  }

  function targets(parameters, conditions, config) {
    const mechanism = config.mechanism;
    const donor = config.backgrounds.findIndex(function (background) {
      return background.id === mechanism.initial_donor;
    });
    const state = initialState(
      config,
      true,
      mechanism.initial_resistant,
      donor < 0 ? 0 : donor
    );
    const referenceEnd = mechanism.antibiotic_start;
    const referenceTrajectory = simulate(
      parameters,
      referenceEnd,
      mechanism.carriage_horizon,
      mechanism.ode_step,
      config,
      state
    );
    const currentTrajectory = simulate(
      parameters,
      conditions.antibioticDays,
      mechanism.carriage_horizon,
      mechanism.ode_step,
      config,
      state
    );
    const withinReference = nextGeneration(parameters, referenceEnd, config);
    const withinCurrent = nextGeneration(parameters, conditions.antibioticDays, config);
    const betweenReference = betweenHost(
      referenceTrajectory,
      conditions.contactRate,
      1,
      config
    );
    const betweenCurrent = betweenHost(
      currentTrajectory,
      conditions.contactRate,
      conditions.susceptibleFraction,
      config
    );
    return {
      r0Within: withinReference.reproduction,
      reWithin: withinCurrent.reproduction,
      r0Between: betweenReference.reproduction,
      reBetween: betweenCurrent.reproduction,
      productBetween: betweenCurrent.productComparator,
      carriageDuration: betweenCurrent.duration,
      matrix: withinCurrent.matrix,
      exposure: withinCurrent.exposure
    };
  }

  return {
    targets: targets,
    nextGeneration: nextGeneration,
    spectralRadius: spectralRadius,
    simulate: simulate
  };
});
