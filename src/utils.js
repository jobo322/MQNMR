function getDistFromJ(jCoupling) {
    let center = -jCoupling.reduce((a, b) => a + b, 0) / 2; //@TODO check if always it is true;
    console.log(center)
    let dist = [center];
    for (let i = 0; i < jCoupling.length; i++) {
      center += jCoupling[i];
      dist.push(center);
    }
    return dist;
  }
  
  function getCandidates(peaks, jcp, pattern, candidates, options = {}) {
    let debug = false;
    if (candidates.length === 0) return null;
    if (
      pattern.length === 0 ||
      candidates.some(e => e.indexs.length === pattern.length)
    ) {
      let { delta, range, nH } = options;
      return candidates.map(cand => {
        let indexs = cand.indexs;
        let toExport = {
          peaks: indexs.map(index => {
            return peaks[index];
          }),
          score: cand.score,
          delta,
          range,
          nH
        };
        return toExport;
      });
    }
    let len = peaks.length;
    let newCandidates = [];
    for (let i = 0; i < candidates.length; i++) {
      let { indexs, score } = candidates[i];
      let index = indexs[indexs.length - 1];
      let iPattern = indexs.length - 1;
      if (debug) console.log("ipattern ", iPattern);
      let maxDiff = jcp[iPattern] + 0.001;
      for (let j = index + 1; j < len; j++) {
        let c = Math.abs(peaks[index].x - peaks[j].x);
        let diff = Math.abs(c - jcp[iPattern]) / jcp[iPattern];
  
        if (debug) console.log(`index ${index} j: ${j} length: ${len}`);
        if (debug) console.log("ref / exp / diff ---> ", jcp[iPattern], c, diff);
  
        if (c > maxDiff) {
          if (debug) console.log(c, maxDiff);
          if (debug) console.log("entra break per diff J");
          break;
        }
        if (diff < 0.05) {
          if (debug) console.log("pasa diff J");
          let RIP = pattern[iPattern] / pattern[iPattern + 1];
          let RIC = peaks[index].y / peaks[j].y;
  
          if (debug) {
            console.log(
              "pattern,",
              pattern[iPattern],
              pattern[iPattern + 1],
              RIP
            );
          }
          if (debug) {
            console.log("experimental ", peaks[index].y, peaks[j].y, RIC);
          }
  
          let diffRI = Math.abs(RIP - RIC) / RIP;
  
          if (debug) console.log(`diffRI -> ${diffRI} max: 0.30`);
  
          if (diffRI < 0.1) {
            score += 1 - diffRI;
            newCandidates.push({ indexs: indexs.concat([j]), score });
            if (debug) console.log("candidate added");
          }
        } else {
          if (debug) console.log("Descartado diff > max");
        }
      }
    }
    return getCandidates(peaks, jcp, pattern, newCandidates, options);
  }
  
  function getCandidatesByJ(peaks, signal, options) {
    console.log('signal')
    console.log('sinal')
    let { field } = options;
    let delta = signal.delta;
    let pattern = signal.pattern;
    let jCoupling = signal.jCoupling.map(j => j / field);
    let candidates = [];
    let candOptions = { delta, range: signal.range, nH: signal.integral };
    peaks.forEach((_peak, pp, arr) => {
      let cand = getCandidateByJ(
        arr,
        jCoupling,
        pattern,
        [{ indexs: [pp], score: 0 }],
        candOptions
      ); // generate combinations from J
      if (cand !== null) candidates = candidates.concat(cand);
    });
    return candidates;
  }
  function getCandidateByJ(peaks, jcp, pattern, candidates, options = {}) {
    let debug = false;
    if (candidates.length === 0) return null;
    if (
      pattern.length === 0 ||
      candidates.some(e => e.indexs.length === pattern.length)
    ) {
      let { delta, range, nH } = options;
      return candidates.map(cand => {
        let indexs = cand.indexs;
        let toExport = {
          peaks: indexs.map(index => {
            return peaks[index];
          }),
          score: cand.score,
          delta,
          range,
          nH
        };
        return toExport;
      });
    }
    let len = peaks.length;
    let newCandidates = [];
    candidates.forEach((candidate, i, arr) => {
      let { indexs, score } = candidate;
      let index = indexs[indexs.length - 1];
      let iPattern = indexs.length - 1;
      if (debug) console.log("ipattern ", iPattern);
      let maxDiff = jcp[iPattern] + 0.001;
      for (let j = index + 1; j < len; j++) {
        let c = Math.abs(peaks[index].x - peaks[j].x);
        let diff = Math.abs(c - jcp[iPattern]) / jcp[iPattern];
  
        if (debug) console.log(`index ${index} j: ${j} length: ${len}`);
        if (debug) console.log("ref / exp / diff ---> ", jcp[iPattern], c, diff);
  
        if (c > maxDiff) {
          if (debug) console.log(c, maxDiff);
          if (debug) console.log("entra break per diff J");
          break;
        }
        if (diff < 0.05) {
          newCandidates.push({ indexs: indexs.concat([j]), score });
        } else {
          if (debug) console.log("Descartado diff > max");
        }
      }
    });
    return getCandidateByJ(peaks, jcp, pattern, newCandidates, options);
  }
  function getNoiseLevel(y) {
    let mean = 0;
  
    let stddev = 0;
    let length = y.length;
    for (let i = 0; i < length; ++i) {
      mean += y[i];
    }
    mean /= length;
    let averageDeviations = new Array(length);
    for (let i = 0; i < length; ++i) {
      averageDeviations[i] = Math.abs(y[i] - mean);
    }
    averageDeviations.sort((a, b) => a - b);
    if (length % 2 === 1) {
      stddev = averageDeviations[(length - 1) / 2] / 0.6745;
    } else {
      stddev =
        (0.5 *
          (averageDeviations[length / 2] + averageDeviations[length / 2 - 1])) /
        0.6745;
    }
  
    return stddev;
  }
  function getCombinations(arrayOfArrays, options = {}) {
    if (Object.prototype.toString.call(arrayOfArrays) !== "[object Array]") {
      throw new Error("combinations method was passed a non-array argument");
    }
  
    let combinations = [];
    let numOfCombos = arrayOfArrays.length ? 1 : 0;
    let arrayOfArraysLength = arrayOfArrays.length;
    for (var n = 0; n < arrayOfArraysLength; ++n) {
      if (Object.prototype.toString.call(arrayOfArrays[n]) !== "[object Array]") {
        throw new Error("combinations method was passed a non-array argument");
      }
      numOfCombos = numOfCombos * arrayOfArrays[n].length;
    }
  
    for (var x = 0; x < numOfCombos; ++x) {
      let carry = x;
      let comboKeys = [];
      let combo = [];
  
      for (var i = 0; i < arrayOfArraysLength; ++i) {
        comboKeys[i] = carry % arrayOfArrays[i].length;
        carry = Math.floor(carry / arrayOfArrays[i].length);
      }
      for (var i = 0; i < comboKeys.length; ++i) {
        combo.push(arrayOfArrays[i][comboKeys[i]]);
      }
      combinations.push(combo);
    }
    return combinations;
  }
  
  function getCombinationsScored(arrayOfArrays, options = {}) {
    let { eretic, sqrtPI, intPattern, parentPort } = options;
      parentPort.postMessage('intPattern ' + JSON.stringify(intPattern));
    if (Object.prototype.toString.call(arrayOfArrays) !== "[object Array]") {
      throw new Error("combinations method was passed a non-array argument");
    }
  
    let combinations = [];
    let numOfCombos = arrayOfArrays.length ? 1 : 0;
    let arrayOfArraysLength = arrayOfArrays.length;
    for (var n = 0; n < arrayOfArraysLength; ++n) {
      if (Object.prototype.toString.call(arrayOfArrays[n]) !== "[object Array]") {
        throw new Error("combinations method was passed a non-array argument");
      }
      numOfCombos = numOfCombos * arrayOfArrays[n].length;
    }
  
    for (var x = 0; x < numOfCombos; ++x) {
      let carry = x;
      let comboKeys = [];
      let combo = [];
  
      for (var i = 0; i < arrayOfArraysLength; ++i) {
        comboKeys[i] = carry % arrayOfArrays[i].length;
        carry = Math.floor(carry / arrayOfArrays[i].length);
      }
      for (var i = 0; i < comboKeys.length; ++i) {
        combo.push(arrayOfArrays[i][comboKeys[i]]);
      }
      let integrals = [];
      let mean = 0;
      let similarityPatternScore = 0;
      let deltaScore = 0;
      combo.forEach((c, i, arr) => {
        let peaks = c.peaks;
        similarityPatternScore += c.score;
        parentPort.postMessage('score ' + c.score)
        if (c.deltaScore) {
          deltaScore += c.deltaScore;
        }
        let integral =
          (peaks.reduce((a, b) => {
            let peak = b;
            parentPort.postMessage(peak.x + ' ' + peak.y + ' ' + peak.width)
            return (
              peak.y * peak.width * sqrtPI * (1 - peak.xL + peak.xL * sqrtPI) + a
            );
          }, 0) /
            intPattern[i] /
            eretic) *
          10;
        arr[i].integral = integral;
        mean += integral;
        parentPort.postMessage(integral)
  
        integrals.push(integral);
      });
      similarityPatternScore /= combo.length;
      if (deltaScore !== 0) deltaScore /= combo.length;
      // console.log('combo', combo)
      // combo.forEach(c => {
      //     console.log(c.peaks.map(ee => ([ee.y,ee.x])))
      // })
      // console.log(JSON.stringify(integrals))
      // console.log('\n');
      // if (integrals.some(e => e < 0.05)) continue;
      mean /= combo.length;
      let std = integrals.reduce((a, b) => Math.pow(b - mean, 2) + a, 0);
      std = Math.sqrt(std / combo.length);
      parentPort.postMessage("mean " + mean);
      parentPort.postMessage("std " + std);
      let less = std / mean || 0;
      let score = 1 - less;
      // if (std / mean > 0.05) continue;
      let result = {
        signals: combo,
        meanIntegral: mean,
        similarityPatternScore,
        deltaScore,
        IntegralScore: score,
        score: score * 0.1 + similarityPatternScore * 0.1 + deltaScore * 0.8
      };
      // combo.meanIntegral = mean;
      // combo.similarityPatternScore = similarityPatternScore;
      // combo.IntegralScore = score;
      // combo.deltaScore = deltaScore;
      // combo.score = ;
      // console.log('pasa')
      // console.log('\n mean ', std)
      // console.log(integrals);
      // console.log(combo.map(e => e.score))
      combinations.push(result);
    }
    return combinations;
  }
  
  function checkIntegrals(expPeaks, toAssignPeaks, toCombine, options) {
    let combinations = getCombinationsScored(toCombine, options);
    return combinations;
  }
  function getDelta(peaks) {
    let delta = peaks.reduce((a, b) => a + b.x, 0) / peaks.length;
    return delta;
  }

  
  
  module.exports = {
    getDistFromJ,
    getCandidates,
    getCandidatesByJ,
    getNoiseLevel,
    getCombinations,
    getCombinationsScored,
    checkIntegrals,
    getDelta
  };
  