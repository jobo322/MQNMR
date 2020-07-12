'use strict';

const { parentPort } = require('worker_threads');

const { xyzAutoPeaksPicking } = require('nmr-processing');

module.exports = function getPeaks2D(zMatrix, cluster, options = {}) {
  parentPort.postMessage({ mix: zMatrix.minX, maxx: zMatrix.maxX });
  let data = reduceByX(zMatrix, cluster);
  let peakList2D = xyzAutoPeaksPicking(data, options);
  peakList2D.sort((a, b) => a.shiftX - b.shiftX);
  // parentPort.postMessage(peakList2D.map(e => e.shiftX))
  // parentPort.postMessage(`peak len ${peakList2D.length}`);
  let groups = [];
  let lastIndex = 0;
  const tolerance = 0.005;
  while (lastIndex < peakList2D.length) {
    let result = getGroup(peakList2D, tolerance, lastIndex);
    groups.push(result.group);
    lastIndex = result.lastIndex;
  }
  // parentPort.postMessage(`groups len ${groups.length}, `)
  let clusters = new Array(groups.length);
  for (let i = 0; i < groups.length; i++) {
    let group = groups[i];
    clusters[i] = {
      signals: group.map((s) => ({
        intensity: Math.max(...s.peaks.map((p) => p.z)),
        xx: s.shiftX,
        x: s.shiftY,
        width: 0.005,
      })),
      deltaX: group.reduce((a, b) => a + b.shiftX, 0) / group.length,
    };
  }
  return { jRes: clusters, data2D: data };
};

function getGroup(peaks, tolerance, startIndex) {
  let lastIndex = startIndex + 1;
  let group = [peaks[startIndex]];
  for (let i = lastIndex; i < peaks.length; i++) {
    if (Math.abs(peaks[startIndex].shiftX - peaks[i].shiftX) <= tolerance) {
      group.push(peaks[i]);
      lastIndex = i + 1;
    } else {
      lastIndex = i;
      break;
    }
  }
  return {
    group,
    lastIndex,
  };
}

function reduceByX(zMatrix, cluster) {
  let { minX, maxX, z } = zMatrix;
  let { from, to } = cluster.range || cluster;
  if (from > to) [from, to] = [to, from];

  let dx = (maxX - minX) / (z[0].length - 1);
  let range = getIndexs([from, to], minX, dx, z[0].length);

  return Object.assign({}, zMatrix, {
    z: z.map((row) => row.slice(range.index[0], range.index[1] + 1)),
    minX: range.value[0],
    maxX: range.value[1],
  });
}

function getIndexs(fromTo, minX, dx, length) {
  let result = [];
  for (let value of fromTo) {
    let index = Math.round((value - minX) / dx);
    result.push(index);
  }
  let n = result[1] - result[0] + 1;
  if (n === 0 || (n & (n - 1)) !== 0) {
    let residual = Math.pow(2, Math.round(Math.log2(n))) - n;
    let newFrom = (residual / 2) << (0 + (residual % 2));
    let newTo = residual - newFrom;
    if (result[0] - newFrom < 0) {
      if (result[1] + residual < length) {
        result[1] += residual;
      } else {
        throw new Error('not possible to get range 2^M with M integer');
      }
    } else if (result[1] + newTo >= length) {
      if (result[0] - residual >= 0) {
        result[0] -= residual;
      } else {
        throw new Error('not possible to get range 2^M with M integer');
      }
    } else {
      result[0] -= newFrom;
      result[1] += newTo;
    }
  }

  return {
    index: result,
    value: result.map((index) => index * dx + minX),
  };
}
