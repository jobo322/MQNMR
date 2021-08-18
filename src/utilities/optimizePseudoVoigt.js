'use strict';

const LM = require('ml-curve-fitting');
const Matrix = require('ml-matrix');

function optimizePseudoVoigtSum(xy, group, options = {}) {
  let {
    percentage = 0,
    LMOptions = [3, 100, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1],
    noiseLevel = 0,
  } = options;

  let xy2 = parseData(xy, percentage);
  if (xy2 === null || xy2[0].rows < 3) {
    return null; //Cannot run an optimization with less than 3 points
  }

  let t = xy2[0];
  let yData = xy2[1];
  let maxY = xy2[2];
  let nbPoints = t.rows;
  let weight = [nbPoints / Math.sqrt(yData.dot(yData))];
  let consts = [[noiseLevel / maxY]]; // optional vector of constants
  let nL = group.length;
  let pInit = new Matrix(nL * 4, 1);
  let pMin = new Matrix(nL * 4, 1);
  let pMax = new Matrix(nL * 4, 1);
  let dx = new Matrix(nL * 4, 1);
  let dt = Math.abs(t[0][0] - t[1][0]);

  for (let i = 0; i < nL; i++) {
    pInit[i][0] = group[i].x;
    pInit[i + nL][0] = group[i].y / maxY;
    pInit[i + 2 * nL][0] = group[i].width;
    pInit[i + 3 * nL][0] = 0.5;

    pMin[i][0] = group[i].x - dt;
    pMin[i + nL][0] = 0;
    pMin[i + 2 * nL][0] = group[i].width / 4;
    pMin[i + 3 * nL][0] = 0;

    pMax[i][0] = group[i].x + dt;
    pMax[i + nL][0] = 1.25;
    pMax[i + 2 * nL][0] = group[i].width * 2;
    pMax[i + 3 * nL][0] = 1;

    dx[i][0] = -dt / 1000;
    dx[i + nL][0] = -1e-3;
    dx[i + 2 * nL][0] = -dt / 1000;
    dx[i + 3 * nL][0] = -0.001;
  }
  dx = -Math.abs(t[0][0] - t[1][0]) / 10000;

  let pFit = LM.optimize(
    sumOfPseudoVoigt,
    pInit,
    t,
    yData,
    weight,
    dx,
    pMin,
    pMax,
    consts,
    LMOptions,
  );
  pFit = pFit.p;

  //Put back the result in the correct format
  let result = new Array(nL);
  for (let i = 0; i < nL; i++) {
    result[i] = [
      pFit[i],
      [pFit[i + nL][0] * maxY],
      pFit[i + 2 * nL],
      pFit[i + 3 * nL],
    ];
  }
  return result;
}

module.exports = {
  optimizePseudoVoigtSum,
};

function sumOfPseudoVoigt(t, p, c) {
  let nL = p.length / 4;
  let factorG;
  let factorL;
  let cols = t.rows;
  let p2;
  let result = Matrix.eye(t.length, 1, c[0][0]);
  // console.log('the first is %s', result[0][0])
  for (let i = 0; i < nL; i++) {
    let xL = p[i + nL * 3][0];
    let xG = 1 - xL;
    p2 = Math.pow(p[i + nL * 2][0], 2);
    factorL = xL * p[i + nL][0] * p2;
    factorG = xG * p[i + nL][0];
    for (let j = 0; j < cols; j++) {
      result[j][0] +=
        factorG * Math.exp(-Math.pow(t[j][0] - p[i][0], 2) / p2) +
        factorL / (Math.pow(t[j][0] - p[i][0], 2) + p2);
    }
  }
  return result;
}
/**
 *
 * Converts the given input to the required x, y column matrices. y data is normalized to max(y)=1
 * @param xy
 * @returns {*[]}
 */
function parseData(xy, threshold) {
  let nbSeries = xy.length;
  let t = null;
  let yData = null;
  let x;
  let y;
  let maxY = 0;
  let nbPoints;

  if (nbSeries === 2) {
    //Looks like row wise matrix [x,y]
    nbPoints = xy[0].length;
    //if(nbPoints<3)
    //    throw new Exception(nbPoints);
    //else{
    t = new Array(nbPoints); //new Matrix(nbPoints,1);
    yData = new Array(nbPoints); //new Matrix(nbPoints,1);
    x = xy[0];
    y = xy[1];
    if (typeof x[0] === 'number') {
      for (let i = 0; i < nbPoints; i++) {
        t[i] = x[i];
        yData[i] = y[i];
        if (y[i] > maxY) maxY = y[i];
      }
    } else {
      //It is a colum matrix
      if (typeof x[0] === 'object') {
        for (let i = 0; i < nbPoints; i++) {
          t[i] = x[i][0];
          yData[i] = y[i][0];
          if (y[i][0] > maxY) maxY = y[i][0];
        }
      }
    }

    //}
  } else {
    //Looks like a column wise matrix [[x],[y]]
    nbPoints = nbSeries;
    //if(nbPoints<3)
    //    throw new SizeException(nbPoints);
    //else {
    t = new Array(nbPoints); //new Matrix(nbPoints, 1);
    yData = new Array(nbPoints); //new Matrix(nbPoints, 1);
    for (let i = 0; i < nbPoints; i++) {
      t[i] = xy[i][0];
      yData[i] = xy[i][1];
      if (yData[i] > maxY) maxY = yData[i];
    }
    //}
  }

  for (let i = 0; i < nbPoints; i++) {
    yData[i] /= maxY;
  }

  if (threshold) {
    for (let i = nbPoints - 1; i >= 0; i--) {
      if (yData[i] < threshold) {
        yData.splice(i, 1);
        t.splice(i, 1);
      }
    }
  }

  if (t.length > 0) {
    return [new Matrix([t]).transpose(), new Matrix([yData]).transpose(), maxY];
  }
  return null;
}
