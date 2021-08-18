'use strict';

function getDistFromJ(jCoupling) {
  let center = -jCoupling.reduce((a, b) => a + b, 0) / 2; //@TODO check if always it is true;
  let dist = [center];
  for (let i = 0; i < jCoupling.length; i++) {
    center += jCoupling[i];
    dist.push(center);
  }
  return dist;
}

function getDelta(peaks) {
  let delta = peaks.reduce((a, b) => a + b.x, 0) / peaks.length;
  return delta;
}

module.exports = {
  getDistFromJ,
  getDelta,
};
