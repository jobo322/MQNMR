const { gsd } = require('ml-gsd');
const spectraProcessing = require('ml-spectra-processing');

function getPeaks(xy, cluster, options) {
  const { defaultOptions } = options;
  let { from, to } = cluster.range || cluster;
  if (from > to) [from, to] = [to, from];
  let reduceOptions = { from, to, nbPoints: Number.MAX_SAFE_INTEGER };
  let data = spectraProcessing.XY.reduce(xy, reduceOptions);
  let gsdOptions = cluster.gsdOptions || {};
  gsdOptions = Object.assign({}, defaultOptions, gsdOptions);

  let peakList = gsd(data.x, data.y, gsdOptions); //should reduce number of peaks by noise level threshold

  return {
    peakList,
    data,
  };
}

module.exports = getPeaks;
