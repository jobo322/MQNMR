'use strict';

const getCombinationsScored = require('./getCombinationsScored');

function checkIntegrals(expPeaks, toAssignPeaks, toCombine, options) {
  let combinations = getCombinationsScored(toCombine, options);
  return combinations;
}

module.exports = checkIntegrals;
