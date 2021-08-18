'use strict';

const optimizePeaks = require('./optimizePeaks');

function runOptimization(xy, peaks, candidates, optOptions) {
  for (let i = 0; i < candidates.length; i++) {
    let candPeaks = candidates[i].peaks;
    let candPeakIndexs = candPeaks.map((e) => e.index);
    //parentPort.postMessage(`candidatesPeakIndexs ${JSON.stringify(candPeakIndexs)}`);
    let first = candPeaks[0];
    let last = candPeaks[candPeaks.length - 1];
    let from = first.x - first.width * 4;
    let to = last.x + last.width * 4;
    let filteredPeaks = peaks.filter((peak) => {
      let w3 = peak.width * 3;
      return peak.x + w3 >= from && peak.x - w3 <= to;
    });
    let optPeaks = optimizePeaks(filteredPeaks, xy.x, xy.y, optOptions);

    candidates[i].peaks = optPeaks.filter((e) =>
      candPeakIndexs.includes(e.index),
    );
    candidates[i].optPeaks = optPeaks;
  }

  return candidates;
};

module.exports = runOptimization;