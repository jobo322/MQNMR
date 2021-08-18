const { getPeaks, optimizePeaks } = require('../../utilities/utils');

module.exports = function (ps, xy, options) {
  let { sqrtPI, defaultOptions } = options;

  let cluster = ps.peaks[0];
  let delta = cluster.delta;
  let signal = ps.toSearch.find((e) => e.delta === delta);
  if (!signal) return;
  let shift = { delta, nH: signal.integral, selected: [], integral: -0.1 };

  let { peakList, data } = getPeaks(xy, cluster, options);

  let optOptions = Object.assign({}, defaultOptions, cluster.gsdOptions);
  let optPeaks = optimizePeaks(peakList, data.x, data.y, optOptions);
  let selectedPeak = getMaxPeak(optPeaks);
  integral =
    (selectedPeak.y *
      selectedPeak.width *
      sqrtPI *
      (1 - selectedPeak.xL + selectedPeak.xL * sqrtPI)) /
    ps.peaks[0].integral;

  shift = Object.assign({}, shift, {
    integral,
    selected: [selectedPeak],
    optPeaks,
  });
  return { metabolite: { signals: [shift] } };
};

function getMaxPeak(peaks) {
  let maxPeak = peaks[0];
  for (let i = 1; i < peaks.length; i++) {
    if (maxPeak.y < peaks[i].y) {
      maxPeak = peaks[i];
    }
  }
  return maxPeak;
}
