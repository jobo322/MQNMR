'use strict';

module.exports = function chemicalShift(peaks, mask) {
  let sum = 0;
  let cs = 0;
  let area;
  if (mask) {
    for (let i = 0; i < peaks.length; i++) {
      if (mask[i] === true) {
        area = getArea(peaks[i]);
        sum += area;
        cs += area * peaks[i].x;
      }
    }
  } else {
    for (let i = 0; i < peaks.length; i++) {
      area = getArea(peaks[i]);
      sum += area;
      cs += area * peaks[i].x;
    }
  }
  return cs / sum;
};

function getArea(peak) {
  return Math.abs(peak.y * peak.width * 1.57); // 1.772453851);
}
