'use strict';

function getCombinations(arrayOfArrays, options = {}) {
  if (Object.prototype.toString.call(arrayOfArrays) !== '[object Array]') {
    throw new Error('combinations method was passed a non-array argument');
  }

  let combinations = [];
  let numOfCombos = arrayOfArrays.length ? 1 : 0;
  let arrayOfArraysLength = arrayOfArrays.length;
  for (let n = 0; n < arrayOfArraysLength; ++n) {
    if (Object.prototype.toString.call(arrayOfArrays[n]) !== '[object Array]') {
      throw new Error('combinations method was passed a non-array argument');
    }
    numOfCombos = numOfCombos * arrayOfArrays[n].length;
  }

  for (let x = 0; x < numOfCombos; ++x) {
    let carry = x;
    let comboKeys = [];
    let combo = [];

    for (let i = 0; i < arrayOfArraysLength; ++i) {
      comboKeys[i] = carry % arrayOfArrays[i].length;
      carry = Math.floor(carry / arrayOfArrays[i].length);
    }
    for (let i = 0; i < comboKeys.length; ++i) {
      combo.push(arrayOfArrays[i][comboKeys[i]]);
    }
    combinations.push(combo);
  }
  return combinations;
}

module.exports = getCombinations;
