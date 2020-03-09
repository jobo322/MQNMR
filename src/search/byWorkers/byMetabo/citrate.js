const utils = require("../utils");
const spectraProcessing = require("ml-spectra-processing");
const { gsd } = require("ml-gsd");
const optimize = require('ml-optimize-lorentzian');
const debug = false;
module.exports = function(ps, xy, options) {
  let { field, delta, peaksToSearch, parentPort, defaultOptions, sqrtPI } = options; 
  let toCombine = [];
  let multiplicities = [];
  let index = 0;
  let optPeaks, data;
  ps.peaks.forEach(cluster => {
    let delta = cluster.delta;
    let signal = ps.toSearch.find(e => e.delta === delta);
    if (!signal) return;
    // intPattern.push(signal.integral); // part of checkIntegral filter
    let shift = { delta, nH: signal.integral, selected: [], integral: -0.1 };

    let { from, to } = cluster.range || cluster;
    if (from > to) [from, to] = [to, from];
    let reduceOptions = { from, to };
    var { x, y } = xy;
    data = spectraProcessing.XY.reduce(x, y, reduceOptions);
    // parentPort.postMessage("pasa reduce");
    let gsdOptions = cluster.gsdOptions || {};
    let options = Object.assign({}, defaultOptions, gsdOptions);
    if (debug)
      parentPort.postMessage(
        `data length ${JSON.stringify(data.y.length % 2)}`
      );
    var peakList = gsd(data.x, data.y, options);
    peakList.forEach((p, pi, arr) => {
      arr[pi].index = index++;
    });
    // parentPort.postMessage("peakList length1 " + peakList.length);
    if (peakList.length === 0) {
      parentPort.postMessage(
        `entry: ${entry} range: ${JSON.stringify({ from, to })}`
      );
      return;
    }
    optPeaks = peakList; //optimizePeaks(peakList, data.x, data.y, options);
    from = signal.range.from || signal.from; //@TODO: check it
    to = signal.range.to || signal.to; //@TODO: check it
    if (from > to) [from, to] = [to, from];
    let peaks = optPeaks.filter(e => e.x < to && e.x > from);

    if (signal.getNoise) {
      let noiseLevel = utils.getNoiseLevel(peaks.map(e => e.y));
      peaks = peaks.filter(e => e.y > noiseLevel);
    }
    // parentPort.postMessage("peakList length2 " + peaks.length);
    // parentPort.postMessage(field);
    let candidates = utils.getCandidatesByJ(peaks, signal, { field });
    let interfCand = {};
    interfCand[ps.name] = candidates.slice();
    if (candidates.length > 0) {
      if (signal.interferences) {
        signal.interferences.forEach(interference => {
          let dataTemp = peaksToSearch.find(
            pst => pst.name === interference.name
          );
          let signalTemp = dataTemp.toSearch.find(
            e => e.delta === interference.delta
          );
          if (signalTemp.jCoupling.length < 1) return; // filtro para no generar candidates con singuletes
          multiplicities.push(signalTemp.jCoupling.length);

          // parentPort.postMessage(JSON.stringify(jCoupling));
          candidates = utils.getCandidatesByJ(peaks, signalTemp, { field });
          if (candidates.length > 0)
            interfCand[interference.name] = candidates.slice();
        });
      }
      let candMatrix = Object.keys(interfCand).map(e => {
        interfCand[e].forEach((_c, ic, cArr) => {
          cArr[ic].name = e;
        });
        return interfCand[e];
      });
      let combinations = utils.getCombinations(candMatrix);

      // parentPort.postMessage("finalComb " + combinations.length);
      // parentPort.postMessage("finalComb1 " + combinations[0].length);

      let range = 1.2 / field;
      let filteredCombinations = combinations.filter(combination => {
        // parentPort.postMessage(JSON.stringify(combination))
        let limit = combination.length;
        let main = combination[0].peaks;
        let includeIt = limit === 1 ? true : false;
        for (let ic = 1; ic < limit; ic++) {
          includeIt = main.some(mp => {
            return combination[ic].peaks.some(cp => {
              let diff = Math.abs(cp.x - mp.x);
              return diff <= range;
            });
          });
        }
        // parentPort.postMessage(includeIt)
        return includeIt;
      });

      if (filteredCombinations.length > 0) combinations = filteredCombinations;
      toCombine.push(combinations);
    }
    // parentPort.postMessage(
    //   `toCombine ${JSON.stringify(utils.getCombinations(toCombine))}`
    // );
  });
  toCombine = utils.getCombinations(toCombine);

  candidates = [];
  // parentPort.postMessage(multiplicities);
  if (multiplicities.some(e => e > 0)) {
    toCombine.forEach((combi) => {
      let tempArr = [];
      let guest = [];
      let constants = [[0]];
      let toRescue = [];
      // parentPort.postMessage("entra");
      let candidatesToOptimize = {};
      combi.forEach(combination => {
        // parentPort.postMessage("pasa");
        combination.forEach(cand => {
          let minY = Number.MAX_SAFE_INTEGER;
          let maxW = Number.MIN_SAFE_INTEGER;
          // parentPort.postMessage(cand.name);
          if (!candidatesToOptimize[cand.name])
            candidatesToOptimize[cand.name] = [];
          let dataCandTemp = peaksToSearch.find(pst => pst.name === cand.name);
          let dataSignalTemp = dataCandTemp.toSearch.find(
            e => e.delta === cand.delta
          );
          // parentPort.postMessage(
          //   JSON.stringify(dataSignalTemp.jCoupling.map(j => j / field))
          // );
          let candJcoupling = [];
          let integral = 0;
          cand.peaks.forEach((ee, i, arrCand) => {
            // parentPort.postMessage(JSON.stringify(ee))
            if (i > 0) {
              candJcoupling.push(ee.x - arrCand[i - 1].x);
            }
            if (minY > ee.y) minY = ee.y;
            if (maxW < ee.width) maxW = ee.width;
            integral += ee.y * ee.width * sqrtPI //* (1 - ee.xL + ee.xL * sqrtPI)
            tempArr.push(ee.index);
            // parentPort.postMessage(` ${ee.index} ${ee.x}`);
          });
         
          let distPeaks = utils.getDistFromJ(candJcoupling);
          let pattern = dataSignalTemp.pattern;
          parentPort.postMessage(pattern)
          toRescue.push({ name: cand.name, index: guest.length });
          candidatesToOptimize[cand.name].push({
            delta: cand.delta,
            guest: {
              x: utils.getDelta(cand.peaks),
              integral: integral,
              y: minY,
              width: maxW
            },
            constants: { x: distPeaks, y: pattern, relation: cand.nH }
          });
          // guest.push({ x: utils.getDelta(cand.peaks), y: minY, width: maxW });
          // //parentPort.postMessage(`${JSON.stringify(cand.peaks)} delta: ${utils.getDelta(cand.peaks)}`)
          // parentPort.postMessage(
          //   "cand Jcoupling2 " + JSON.stringify(utils.getDistFromJ(candJcoupling))
          // );
          // constants.push({ x: distPeaks, y: pattern, relation: cand.nH });
          // parentPort.postMessage("sale 1");
        });
      });

      // parentPort.postMessage("tme " + JSON.stringify(tempArr));
      // // se tiene todos los candidates de cada combinacion final agrupados en metabolitos, de esta manera se
      // // debe encontrar la manera de reducir el numero de parametros a optimizar.
      
      // let filteredPeaks = optPeaks.filter(e => !tempArr.includes(e.index)); // @ICHANGE
      // parentPort.postMessage("filteredPeaks length " + filteredPeaks.length);

      // filteredPeaks.forEach((fp, i) => {
      //   candidatesToOptimize['singlet' + i] = [{
      //     guest: fp,
      //     constants: { x: [0], y: [1], relation: 1 }
      //   }]
      // });

      // parentPort.postMessage('candidatesToOptimize \n' + JSON.stringify(candidatesToOptimize))
      // parentPort.postMessage(`rescue ${JSON.stringify(toRescue)}`);
      // parentPort.postMessage("antes de optimize");
      // parentPort.postMessage("data length " + data.x.length);
      // parentPort.postMessage(JSON.stringify(guest.length));
      // parentPort.postMessage(JSON.stringify(constants.length));
      // let optimizedPeaks = [];
      
      try {
        optimizedPeaks = optimize.optimizeClusters(
          [data.x, data.y],
          guest,
          {
            candidates: candidatesToOptimize,
            LMOptions: [3, 1000, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1],
          }
        );
      } catch (err) {
        parentPort.postMessage(err);
        parentPort.postMessage(JSON.stringify(err));
      }
      parentPort.postMessage('optimizePeaks')
      parentPort.postMessage(JSON.stringify(candidatesToOptimize) + '\n')
      parentPort.postMessage(JSON.stringify(optimizedPeaks.candidates) + '\n');
      parentPort.postMessage(JSON.stringify(optimizedPeaks.pFit));

      // return
      
      // parentPort.postMessage("pasa optimize");
      // // parentPort.postMessage(JSON.stringify(optimizedPeaks));
      // let result = [];
      // let nL = optimizedPeaks.result[0].length;
      // let keys = ["x", "y", "width", "xL"];
      // for (var j = 0; j < optimizedPeaks.result.length; j++) {
      //   let toPush = {};
      //   for (let k = 0; k < nL; k++) {
      //     toPush[keys[k]] = optimizedPeaks.result[j][k][0];
      //   }
      //   result.push(toPush);
      // }
      // let finalOptimizePeaks = [];
      // parentPort.postMessage("constantes length " + constants.length);
      // for (let ii = 1; ii < constants.length; ii++) {
      //   let c = constants[ii];
      //   let p = result[ii - 1];
      //   // parentPort.postMessage(JSON.stringify(p));
      //   // parentPort.postMessage(JSON.stringify(c));
      //   for (let j = 0; j < c.x.length; j++) {
      //     finalOptimizePeaks.push({
      //       x: p.x + c.x[j],
      //       y: p.y * c.y[j],
      //       width: p.width,
      //       xL: p.xL
      //     });
      //   }
      // }
      // let c = constants[1];
      // let p = result[0];
      // let candPeaks = [];
      // for (let j = 0; j < c.x.length; j++) {
      //   candPeaks.push({
      //     x: p.x + c.x[j],
      //     y: p.y * c.y[j],
      //     width: p.width,
      //     xL: p.xL
      //   });
      // }
      // // parentPort.postMessage('chi ' + JSON.stringify(optimizedPeaks.x2))
      // candidates.push({
      //   peaks: candPeaks,
      //   score: 1 / optimizedPeaks.x2,
      //   nH: signal.nH,
      //   range: signal.range,
      //   delta,
      //   optPeaks: finalOptimizePeaks
      // });
      // // peaks = finalOptimizePeaks.slice();
      // parentPort.postMessage("finalOptimizePeaks " + peaks.length);
      // // parentPort.postMessage("candidates" + JSON.stringify(candidates));
    });
  } else {
    //Assume that there is just singlets
    optPeaks = optimizePeaks(optPeaks, data.x, data.y, options);
    peaks = optPeaks.filter(e => e.x < to && e.x > from);
    peaks.forEach((_peak, pp, arr) => {
      let cand = getCandidates(
        arr,
        jCoupling,
        pattern,
        [{ indexs: [pp], score: 0 }],
        candOptions
      ); // generate combinations from J
      if (cand !== null) candidates = candidates.concat(cand);
    });
    candidates.forEach((cand, ic, cArr) => {
      cArr[ic].optPeaks = optPeaks.slice();
    });
  }
};
