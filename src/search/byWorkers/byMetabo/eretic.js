const { getPeaks, optimizePeaks } = require('../../../utilities/utils');

module.exports = function(ps, xy, options) {
  let { sqrtPI } = options;

  let cluster = ps.peaks[0];
  let delta = cluster.delta;
  let signal = ps.toSearch.find((e) => e.delta === delta);
  if (!signal) return;
  let shift = { delta, nH: signal.integral, selected: [], integral: -0.1 };

  let { peakList, data } = getPeaks(cluster);
  let optPeaks = optimizePeaks(peakList, data.x, data.y, options);
  let selectedPeaks = ps.callback(optPeaks);
  integral =
    selectedPeaks.reduce((a, b) => {
      let peak = b;
      return (
        peak.y * peak.width * sqrtPI * (1 - peak.xL + peak.xL * sqrtPI) + a
      );
    }, 0) / ps.peaks[0].integral;

  return Object.assign({}, shift, {
    integral,
    selected: selectedPeaks,
    optPeaks,
  });
};
