const fs = require('fs');
const path = require('path');

const pathToData = path.resolve('.');
let listSamples = fs.readdirSync(pathToData);
let samples = listSamples.filter((s) => s.match(/newByWorkers_new.*/));
let result = '[';
samples.forEach((s) => {
  let peaks = fs.readFileSync(path.join(pathToData, s), 'utf8');
  result = result.concat(peaks);
});

result = result.slice(0, result.length - 1).concat(']');
fs.writeFileSync('newByWorkers_new.json', result);

function runOptimization(xy, candidates, optOptions) {
  for (let i = 0; i < candidates.length; i++) {
    let candPeaks = candidates[i].peaks;
    let first = candPeaks[0];
    let last = candPeaks[candPeaks.length - 1];
    from = first.x - first.width * 6;
    to = last.x + last.width * 6;
    let filteredPeaks = peaks.filter((peak) => {
      let w3 = peak.width * 3;
      return peak.x + w3 >= from && peak.x - w3 <= to;
    });
    let optPeaks = optimizePeaks(filteredPeaks, xy.x, xy.y, optOptions);
    let peakIndex = candidates[i].peaks[0].index;
    candidates[i].peaks = [optPeaks.find((e) => e.index === peakIndex)];
    candidates[i].optPeaks = optPeaks;
  }
  return candidates;
}
