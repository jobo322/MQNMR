'use strict';

const getPeaks = require('./getPeaks');
const getPeaks2D = require('./getPeaks2D');
const optimizePeaks = require('./optimizePeaks');
const optimize = require('./optimizePseudoVoigt');
const runOptimization = require('./runOptimization');
const getChemicalShift = require('./getChemicalShift');

module.exports = {
  getPeaks,
  getPeaks2D,
  getChemicalShift,
  optimizePeaks,
  optimize,
  runOptimization,
};
