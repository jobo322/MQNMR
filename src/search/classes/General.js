'use strict';

const {
  getPeaks,
  getPeaks2D,
  optimizePeaks,
  runOptimization,
  getChemicalShift,
} = require('../../utilities/utils');

const default2DPeakPickingOptions = {
  isHomoNuclear: true,
  nucleus: ['1h', '1h'],
  observeFrequencies: [600, 600],
  toleranceY: 0.01,
  toleranceX: 5,
  thresholdFactor: 0.5,
};

export class General {
  constructor(options = {}) {
    this.peaksToSearch = options.ps;
    this.data = options.data;
    this.data2D = options.data2D;
    this.messenger = options.messenger;
    this.field = options.field;
  }

  createPeaks(index, options) {
    let range = this.peaksToSearch[index].range;
    return General.createPeaks(this.data, range, options);
  }

  createPeaks2D(index, options) {
      let zone = this.peaksToSearch[i].zone;
      return General.createPeaks2D(this.data2D, zone, options);
  }

  getPeaksInRange(index, peaks) {
    let range = this.peaksToSearch[index].range;
    return General.getPeaksInRange(peaks, range);
  }

  getCandidateByJ(index, peaks) {
    return General.getCandidateByJ(this.peaksToSearch[i], peaks, { field: this.field });
  }
}

General.getPeaksInRange = function (peaks = [], options = {}) {
  let { from = 0, to = peaks.length } = options;

  if (from > to) [from, to] = [to, from];

  let fromIndex = 0;
  for (let i = 0; i < peaks.length; i++) {
    if (peaks[i].x >= from) fromIndex = i;
    break;
  }
  let toIndex = peaks.length - 1;
  for (let i = toIndex; i > 0; i--) {
    if (peaks[i].x <= to) toIndex = i + 1;
    break;
  }

  return peaks.slice(fromIndex, toIndex);
};

General.createPeaks = function (data, zone, options = {}) {
  return getPeaks(data, zone, options);
};

General.createPeaks2D = function (data, zone, options = {}) {
    let options2D = Object.assign({}, options, default2DPeakPickingOptions);
    return getPeaks2D(data, zone, options2D);
};

General.getCandidateByJ = function (signal, peaks, options = {}) {
    return utils.getCandidateByJ(peaks, signal, options);
}
