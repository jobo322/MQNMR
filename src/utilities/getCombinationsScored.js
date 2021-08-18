'use strict';

function getCombinationsScored(arrayOfArrays, options = {}) {
  let { eretic, sqrtPI, intPattern, parentPort } = options;
  // parentPort.postMessage('intPattern ' + JSON.stringify(intPattern));
  if (Object.prototype.toString.call(arrayOfArrays) !== '[object Array]') {
    throw new Error('combinations method was passed a non-array argument');
  }

  let combinations = [];
  let numOfCombos = arrayOfArrays.length ? 1 : 0;
  let arrayOfArraysLength = arrayOfArrays.length;
  for (let n = 0; n < arrayOfArraysLength; ++n) {
    if (Object.prototype.toString.call(arrayOfArrays[n]) !== '[object Array]') {
      throw new Error('combinations method was passed a non-array argument');
    }
    numOfCombos = numOfCombos * arrayOfArrays[n].length;
  }

  for (let x = 0; x < numOfCombos; ++x) {
    let carry = x;
    let comboKeys = [];
    let combo = [];

    for (let i = 0; i < arrayOfArraysLength; ++i) {
      comboKeys[i] = carry % arrayOfArrays[i].length;
      carry = Math.floor(carry / arrayOfArrays[i].length);
    }
    for (let i = 0; i < comboKeys.length; ++i) {
      combo.push(arrayOfArrays[i][comboKeys[i]]);
    }
    let integrals = [];
    let mean = 0;
    let similarityPatternScore = 0;
    let deltaScore = 0;
    combo.forEach((c, i, arr) => {
      let peaks = c.peaks;
      similarityPatternScore += c.score;
      // parentPort.postMessage('score ' + c.score)
      if (c.deltaScore) {
        deltaScore += c.deltaScore;
      }
      // parentPort.postMessage('peak length ' +  peaks.length)
      let integral =
        (peaks.reduce((a, b) => {
          let peak = b;
          // parentPort.postMessage(peak.x + ' ' + peak.y + ' ' + peak.width + ' ' + peak.xL)
          return (
            peak.y * peak.width * sqrtPI * (1 - peak.xL + peak.xL * sqrtPI) + a
          );
        }, 0) /
          intPattern[i] /
          eretic) *
        10;
      arr[i].integral = integral;
      mean += integral;
      integrals.push(integral);
    });
    similarityPatternScore /= combo.length;
    if (deltaScore !== 0) deltaScore /= combo.length;
    // if (integrals.some(e => e < 0.05)) continue;
    mean /= combo.length;
    let std = integrals.reduce((a, b) => Math.pow(b - mean, 2) + a, 0);
    std = Math.sqrt(std / combo.length);
    // parentPort.postMessage("mean " + mean);
    // parentPort.postMessage("std " + std);
    let less = std / mean || 0;
    let score = 1 - less;
    // if (std / mean > 0.05) continue;
    let result = {
      signals: combo,
      meanIntegral: mean,
      similarityPatternScore,
      deltaScore,
      IntegralScore: score,
      score: score * 0.1 + similarityPatternScore * 0.1 + deltaScore * 0.8,
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

module.exports = getCombinationsScored;
