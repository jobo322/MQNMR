const Matrix = require("ml-matrix");
const LM = require("ml-curve-fitting");
// const { direct } = require('ml-direct-optimization');

function optimizePseudoVoigtSum(xy, group, options = {}) {
  var {
    percentage = 0,
    LMOptions = [3, 100, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1],
    noiseLevel = 0,
  } = options;
  
  var xy2 = parseData(xy, percentage);
  if (xy2 === null || xy2[0].rows < 3) {
    return null; //Cannot run an optimization with less than 3 points
  }
  
  var t = xy2[0];
  var yData = xy2[1];
  var maxY = xy2[2];
  var nbPoints = t.rows;
  
  var weight = [nbPoints / Math.sqrt(yData.dot(yData))];
  var consts = [[noiseLevel / maxY]]; // optional vector of constants
  var nL = group.length;
  var pInit = new Matrix(nL * 4, 1);
  var pMin = new Matrix(nL * 4, 1);
  var pMax = new Matrix(nL * 4, 1);
  var dx = new Matrix(nL * 4, 1);
  var dt = Math.abs(t[0][0] - t[1][0]);
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
  var dx = -Math.abs(t[0][0] - t[1][0]) / 10000;
  
  var pFit = LM.optimize(
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
  // var pFit = direct( //still does not work well
  //   sumOfPseudoVoigtArray,
  //   pInit,
  //   t,
  //   yData,
  //   weight,
  //   dx,
  //   pMin,
  //   pMax,
  //   consts,
  //   LMOptions
  // );
  pFit = pFit.p;
  //Put back the result in the correct format
  var result = new Array(nL);
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
}

function sumOfPseudoVoigt(t, p, c) {
    var nL = p.length / 4,
      factorG,
      factorL,
      cols = t.rows,
      p2;
    var result = Matrix.eye(t.length, 1, c[0][0]);
    // console.log('the first is %s', result[0][0])
    for (let i = 0; i < nL; i++) {
      var xL = p[i + nL * 3][0];
      var xG = 1 - xL;
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
  var nbSeries = xy.length;
  var t = null;
  var y_data = null,
    x,
    y;
  var maxY = 0,
    i,
    j;
    
  if (nbSeries == 2) {
    //Looks like row wise matrix [x,y]
    var nbPoints = xy[0].length;
    //if(nbPoints<3)
    //    throw new Exception(nbPoints);
    //else{
    t = new Array(nbPoints); //new Matrix(nbPoints,1);
    y_data = new Array(nbPoints); //new Matrix(nbPoints,1);
    x = xy[0];
    y = xy[1];
    
    if (typeof x[0] === 'number') {
      for (i = 0; i < nbPoints; i++) {
        t[i] = x[i];
        y_data[i] = y[i];
        if (y[i] > maxY) maxY = y[i];
      }
    } else {
      //It is a colum matrix
      if (typeof x[0] === 'object') {
        for (i = 0; i < nbPoints; i++) {
          t[i] = x[i][0];
          y_data[i] = y[i][0];
          if (y[i][0] > maxY) maxY = y[i][0];
        }
      }
    }
    
    //}
  } else {
    //Looks like a column wise matrix [[x],[y]]
    var nbPoints = nbSeries;
    //if(nbPoints<3)
    //    throw new SizeException(nbPoints);
    //else {
    t = new Array(nbPoints); //new Matrix(nbPoints, 1);
    y_data = new Array(nbPoints); //new Matrix(nbPoints, 1);
    for (i = 0; i < nbPoints; i++) {
      t[i] = xy[i][0];
      y_data[i] = xy[i][1];
      if (y_data[i] > maxY) maxY = y_data[i];
    }
    //}
  }
  
  for (i = 0; i < nbPoints; i++) {
    y_data[i] /= maxY;
  }
  
  if (threshold) {
    for (i = nbPoints - 1; i >= 0; i--) {
      if (y_data[i] < threshold) {
        y_data.splice(i, 1);
        t.splice(i, 1);
      }
    }
  }
  
  if (t.length > 0) {
    return [
      new Matrix([t]).transpose(),
      new Matrix([y_data]).transpose(),
      maxY,
    ];
  }
  return null;
}


