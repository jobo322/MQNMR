const spectraProcessing = require('ml-spectra-processing');
const { gsd } = require('ml-gsd');

module.exports.getPeaks = function(cluster) {
  let { from, to } = cluster.range || cluster;
  if (from > to) [from, to] = [to, from];
  let reduceOptions = { from, to };
  let { x, y } = xy;
  data = spectraProcessing.XY.reduce(x, y, reduceOptions);
  // parentPort.postMessage("pasa reduce");
  let gsdOptions = cluster.gsdOptions || {};
  gsdOptions = Object.assign({}, defaultOptions, gsdOptions);

  let peakList = gsd(data.x, data.y, gsdOptions); //should reduce number of peaks by noise level threshold

  return {
    peakList,
    data,
  };
};
