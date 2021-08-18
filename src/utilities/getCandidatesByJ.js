function getCandidatesByJ(peaks, signal, options) {
  let { field } = options;
  let delta = signal.delta;
  let pattern = signal.pattern;
  let jCoupling = signal.jCoupling.map((j) => j / field);
  let candidates = [];
  let candOptions = { delta, range: signal.range, nH: signal.integral };
  peaks.forEach((_peak, pp, arr) => {
    let cand = getCandidateByJ(
      arr,
      jCoupling,
      pattern,
      [{ indexs: [pp], score: 0 }],
      candOptions,
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
  candidates.forEach((candidate) => {
    let { indexs, score } = candidate;
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
        newCandidates.push({ indexs: indexs.concat([j]), score });
      } else {
        if (debug) console.log('Descartado diff > max');
      }
    }
  });
  return getCandidateByJ(peaks, jcp, pattern, newCandidates, options);
}

module.exports = {
  getCandidatesByJ,
  getCandidateByJ,
};
