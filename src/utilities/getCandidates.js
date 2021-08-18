'use strict';

function getCandidates(peaks, jcp, pattern, candidates, options = {}) {
  let debug = false;
  if (candidates.length === 0) return null;
  if (
    pattern.length === 0 ||
    candidates.some((e) => e.indexs.length === pattern.length)
  ) {
    let { delta, range, nH } = options;
    return candidates.map((cand) => {
      let indexs = cand.indexs;
      let toExport = {
        peaks: indexs.map((index) => {
          return peaks[index];
        }),
        score: cand.score,
        delta,
        range,
        nH,
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
    if (debug) console.log('ipattern ', iPattern);
    let maxDiff = jcp[iPattern] + 0.001;
    for (let j = index + 1; j < len; j++) {
      let c = Math.abs(peaks[index].x - peaks[j].x);
      let diff = Math.abs(c - jcp[iPattern]) / jcp[iPattern];

      if (debug) console.log(`index ${index} j: ${j} length: ${len}`);
      if (debug) console.log('ref / exp / diff ---> ', jcp[iPattern], c, diff);

      if (c > maxDiff) {
        if (debug) console.log(c, maxDiff);
        if (debug) console.log('entra break per diff J');
        break;
      }
      if (diff < 0.05) {
        if (debug) console.log('pasa diff J');
        let RIP = pattern[iPattern] / pattern[iPattern + 1];
        let RIC = peaks[index].y / peaks[j].y;

        if (debug) {
          console.log(
            'pattern,',
            pattern[iPattern],
            pattern[iPattern + 1],
            RIP,
          );
        }
        if (debug) {
          console.log('experimental ', peaks[index].y, peaks[j].y, RIC);
        }

        let diffRI = Math.abs(RIP - RIC) / RIP;

        if (debug) console.log(`diffRI -> ${diffRI} max: 0.30`);

        if (diffRI < 0.1) {
          score += 1 - diffRI;
          newCandidates.push({ indexs: indexs.concat([j]), score });
          if (debug) console.log('candidate added');
        }
      } else {
        if (debug) console.log('Descartado diff > max');
      }
    }
  }
  return getCandidates(peaks, jcp, pattern, newCandidates, options);
}

module.exports = getCandidates;
