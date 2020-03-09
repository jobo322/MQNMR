"use strict";

const { parentPort, workerData } = require("worker_threads");
const path = require("path");
const fs = require("fs");
let pathToPredictor = path.resolve(
  path.join("src/prediction", "predictor.js")
);
const byMetabo = require('./byMetabo/index');
const { singletPredictor } = require(pathToPredictor);
const converter = require("jcampconverter");
const spectraProcessing = require("ml-spectra-processing");
const { gsd } = require("ml-gsd");
const utils = require('../utils');
const {
  getDistFromJ,
  getCandidates,
  getCandidatesByJ,
  getNoiseLevel,
  getCombinations,
  getCombinationsScored,
  checkIntegrals,
  getDelta
} = utils;
const optimizePeaks = require("../optimizePeaks");
const optimize = require("ml-optimize-lorentzian");
const defaultOptions = {
  thresholdFactor: 1,
  minMaxRatio: 0.0001,
  broadRatio: 0.0025,
  smoothY: false,
  widthFactor: 4,
  realTop: false,
  fnType: "voigt",
  broadWidth: 0.25,
  sgOptions: { windowSize: 15, polynomial: 4 }
};

var {
  index,
  samples,
  pathToData,
  pathInfo,
  toSearch,
  rangeToOpt,
  sqrtPI,
  field,
  subFix
} = workerData;
const msg = parentPort.postMessage;
var info = require(pathInfo);
const { peaksToSearch } = info;
let filteredPeaksToSearch = peaksToSearch.filter(e => {
  return toSearch.includes(e.name);
});

let debug = false;
let first = true;
let filename = "searchNuevoParadigma2.json";
parentPort.postMessage(`sample length ${samples.length}`);
for (let i = 0; i < samples.length; i++) {
  parentPort.postMessage(`--------------- ${String(index)} - ${i}`);
  let sample = samples[i];
  let entry = sample.replace(/\.[a-z]*/g, "");
  parentPort.postMessage(entry)
  let pathToJcamp = path.join(pathToData, sample);
  var jcamp = fs.readFileSync(pathToJcamp, "utf8");
  var spectrum = converter.convert(jcamp, { xy: true });
  var xy = spectrum.spectra[0].data[0];
  if (xy.x[0] > xy.x[1]) {
    xy.x = xy.x.reverse();
    xy.y = xy.y.reverse();
  }
  parentPort.postMessage(xy.x.length);
  let toExport = { sampleid: entry };
  filteredPeaksToSearch.forEach(ps => {
    parentPort.postMessage(ps.name);
    let toCombine = [];
    let intPattern = [];
    let selectedPeaks = [];
    let integral = 0;
    let metabolite = { signals: [] };
    if (ps.name === 'citrate' || ps.name === 'lactate') {
      parentPort.postMessage('ohalo');
      parentPort.postMessage(JSON.stringify(byMetabo));
      byMetabo.citrate(ps, xy, { field, peaksToSearch, parentPort, defaultOptions, sqrtPI });
      parentPort.postMessage('pasa')
    }
    return
    ps.peaks.forEach(cluster => {
      let delta = cluster.delta;
      let signal = ps.toSearch.find(e => e.delta === delta);
      if (!signal) return
      intPattern.push(signal.integral); // part of checkIntegral filter
      let shift = { delta, nH: signal.integral, selected: [], integral: -0.1 };
      
      let { from, to } = cluster.range || cluster;
      if (from > to) [from, to] = [to, from];
      let reduceOptions = { from, to };
      var { x, y } = xy;
      var data = spectraProcessing.XY.reduce(x, y, reduceOptions);
      parentPort.postMessage('pasa reduce')
      let gsdOptions = cluster.gsdOptions || {};
      let options = Object.assign({}, defaultOptions, gsdOptions);
      if (debug)
        parentPort.postMessage(
          `data length ${JSON.stringify(data.y.length % 2)}`
        );
      var peakList = gsd(data.x, data.y, options);
      parentPort.postMessage("peakList length1 " + peakList.length);
      if (peakList.length === 0) {
        parentPort.postMessage(
          `entry: ${entry} range: ${JSON.stringify({ from, to })}`
        );
        return;
      }
      let optPeaks = peakList;//optimizePeaks(peakList, data.x, data.y, options);
      from = signal.range.from || signal.from; //@TODO: check it
      to = signal.range.to || signal.to; //@TODO: check it
      if (from > to) [ from, to ] = [ to, from ];
      let peaks = optPeaks.filter(e => e.x < to && e.x > from);
      if (ps.name !== "eretic") {
        console.log('---------------------------', ps.name);
        // if (ps.name == 'citrate') {
        //   parentPort.postMessage('entra _______ citratet')
        //   byMetabo.citrate(peaks, signal, toExport, parentPort, {field, delta});
        // }
        let multiplicities = [];
        if (signal.getNoise) {
          let noiseLevel = getNoiseLevel(peaks.map(e => e.y));
          peaks = peaks.filter(e => e.y > noiseLevel);
        }
        multiplicities.push(signal.jCoupling.length);
        parentPort.postMessage("peakList length2 " + peaks.length);
        let candidates = getCandidatesByJ(peaks, signal, { field });
        // parentPort.postMessage('lactate' + JSON.stringify(candidates));
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
              multiplicities.push(signalTemp.jCoupling.length);
              if (signalTemp.jCoupling.length < 1) return; // filtro para no generar candidates con singuletes
              candidates = getCandidatesByJ(peaks, signalTemp, { field })
              if (candidates.length > 0)
                interfCand[interference.name] = candidates.slice();
            });
            // parentPort.postMessage('threonine ' + JSON.stringify(interfCand['threonine']));
          }
          let candMatrix = Object.keys(interfCand).map(e => {
            interfCand[e].forEach((_c, ic, cArr) => {
              cArr[ic].name = e;
            });
            return interfCand[e];
          });
          let combinations = getCombinations(candMatrix);
          let range = 1.2 / field;
          parentPort.postMessage("finalComb " + combinations.length);
          parentPort.postMessage("finalComb1 " + combinations[0].length);

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

          if (filteredCombinations.length > 0)
            combinations = filteredCombinations;
          parentPort.postMessage("finalComb " + combinations.length);
          candidates = [];
          parentPort.postMessage(multiplicities)
          if (multiplicities.some(e => e > 0)) {
            combinations.forEach((combination, icomb) => {
              // if (icomb > 0) return;
              let tempArr = [];
              let guest = [];
              let constants = [[0]];
              combination.forEach(cand => {
                let maxY = Number.MIN_SAFE_INTEGER;
                let maxW = Number.MIN_SAFE_INTEGER;
                let dataCandTemp = peaksToSearch.find(
                  pst => pst.name === cand.name
                );
                let dataSignalTemp = dataCandTemp.toSearch.find(
                  e => e.delta === cand.delta
                );
                parentPort.postMessage(JSON.stringify(dataSignalTemp.jCoupling.map(j => j / field)))
                let candJcoupling = [];
                cand.peaks.forEach((ee,i, arrCand) => {
                   if (i > 0) {
                     candJcoupling.push(ee.x - arrCand[i-1].x)
                   }
                  if (maxY < ee.y) maxY = ee.y;
                  if (maxW < ee.width) maxW = ee.width;
                  tempArr.push(ee.index);
                });
                parentPort.postMessage('cand Jcoupling ' + JSON.stringify(candJcoupling));
                guest.push({ x: getDelta(cand.peaks), y: maxY, width: maxW });
                parentPort.postMessage('cand Jcoupling2 ' + JSON.stringify(getDistFromJ(candJcoupling)));
  
                let distPeaks = getDistFromJ(candJcoupling);
                let pattern = dataSignalTemp.pattern;
                constants.push({ x: distPeaks, y: pattern });
              });
              parentPort.postMessage(peaks.length);
              parentPort.postMessage(peaks[0]);
              let filteredPeaks = optPeaks.filter(e => !tempArr.includes(e.index)); // @ICHANGE
              parentPort.postMessage(
                "filteredPeaks length " + filteredPeaks.length
              );
              filteredPeaks.forEach(fp => {
                guest.push(fp);
                constants.push({ x: [0], y: [1] });
              });
              parentPort.postMessage("antes de optimize");
              parentPort.postMessage('data length ' + data.x.length);
              parentPort.postMessage(JSON.stringify(guest.length));
              parentPort.postMessage(JSON.stringify(constants.length));
              let optimizedPeaks = [];

              try {
                optimizedPeaks = optimize.optimizeClustersPVoigtSum(
                  [data.x, data.y],
                  guest,
                  {
                    consts: constants,
                    LMOptions: [3, 100, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1]
                  }
                );
              } catch (err) {
                parentPort.postMessage(err);
                parentPort.postMessage(JSON.stringify(err));
              }
              parentPort.postMessage("pasa optimize");
              // parentPort.postMessage(JSON.stringify(optimizedPeaks));
              let result = [];
              let nL = optimizedPeaks.result[0].length;
              let keys = ["x", "y", "width", "xL"];
              for (var j = 0; j < optimizedPeaks.result.length; j++) {
                let toPush = {};
                for (let k = 0; k < nL; k++) {
                  toPush[keys[k]] = optimizedPeaks.result[j][k][0];
                }
                result.push(toPush);
              }
              let finalOptimizePeaks = [];
              parentPort.postMessage("constantes length " + constants.length);
              for (let ii = 1; ii < constants.length; ii++) {
                let c = constants[ii];
                let p = result[ii - 1];
                // parentPort.postMessage(JSON.stringify(p));
                // parentPort.postMessage(JSON.stringify(c));
                for (let j = 0; j < c.x.length; j++) {
                  finalOptimizePeaks.push({
                    x: p.x + c.x[j],
                    y: p.y * c.y[j],
                    width: p.width,
                    xL: p.xL
                  });
                }
              }
              let c = constants[1];
              let p = result[0];
              let candPeaks = [];
              for (let j = 0; j < c.x.length; j++) {
                candPeaks.push({
                  x: p.x + c.x[j],
                  y: p.y * c.y[j],
                  width: p.width,
                  xL: p.xL
                });
              }
              // parentPort.postMessage('chi ' + JSON.stringify(optimizedPeaks.x2))
              candidates.push({
                peaks: candPeaks,
                score: 1/optimizedPeaks.x2,
                nH: signal.nH,
                range: signal.range,
                delta,
                optPeaks: finalOptimizePeaks
              });
              // peaks = finalOptimizePeaks.slice();
              parentPort.postMessage("finalOptimizePeaks " + peaks.length);
              // parentPort.postMessage("candidates" + JSON.stringify(candidates));
            });
          } else { //Assume that there is just singlets
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
            })
          }
          parentPort.postMessage('candidates length ' + candidates.length)
          if (candidates.length > 0) {
            if (
              ps.name === "glycine" ||
              ps.name === "formate" ||
              ps.name === "alanine" ||
              ps.name === "trigonelline" ||
              ps.name === "tartrate" ||
              ps.name === "creatine" ||
              ps.name === "succinate"
            ) {
              parentPort.postMessage('---------\n----------\n');
              parentPort.postMessage('eentra');
              let prediction = singletPredictor(toExport, ps.name);
              // console.log(prediction);
              // console.log(candidates.map(e => getDelta(e.peaks)));
              if (prediction !== null) {
                candidates.forEach((e, ii, arr) => {
                  let cs = getDelta(e.peaks);
                  let score = (1 - Math.abs(cs - prediction[0])) * 10;
                  arr[ii].deltaScore = score;
                });
              }
            }
            toCombine.push(candidates);
          }
          parentPort.postMessage("pasa");
        }
        shift.optPeaks = peaks;
      } else {
        optPeaks = optimizePeaks(peakList, data.x, data.y, options);
        selectedPeaks = ps.callback(optPeaks);
        integral =
          selectedPeaks.reduce((a, b) => {
            let peak = b;
            return (
              peak.y * peak.width * sqrtPI * (1 - peak.xL + peak.xL * sqrtPI) +
              a
            );
          }, 0) / ps.peaks[0].integral;
        shift = Object.assign({}, shift, {
          integral, 
          selected: selectedPeaks,
          optPeaks: peaks
        }); // @TODO
      } // es la parte donde se escoje a eretic
      metabolite.signals.push(shift);
    });

    if (ps.name !== "eretic") {
      // if (ps.name === 'creatine') {
      // console.log('toCombine', toCombine);
      // }

      parentPort.postMessage(
        "toCombine length " + JSON.stringify(toCombine.length)
      );
      parentPort.postMessage(
        "toSearch length " + JSON.stringify(ps.toSearch.length)
      );
      if (toCombine.length !== ps.toSearch.length) toCombine = [];
      let eretic = toExport.eretic;
      let finalCandidates = getCombinationsScored(toCombine, {
        sqrtPI,
        eretic: eretic.meanIntegral,
        intPattern,
        parentPort
      });
      if (finalCandidates.length === 0) {
        parentPort.postMessage(`sin candidatos ${sample}  ${i}`);
      } else {
        selectedPeaks = [];
        finalCandidates.sort((a, b) => b.score - a.score);
        parentPort.postMessage('finalCandidates');
        parentPort.postMessage(JSON.stringify(finalCandidates.map(a => a.score)))
        parentPort.postMessage(JSON.stringify(finalCandidates.map(a => a.IntegralScore)))
        parentPort.postMessage(JSON.stringify(finalCandidates.map(a => a.similarityPatternScore)))
        let index = 0;
        index = index >= finalCandidates.length ? 0 : index;
        finalCandidates[index].signals.forEach((c, i, arr) => {
          let peaks = c.peaks;
          let delta = c.delta;
          let shift = metabolite.signals.find(e => e.delta === delta);
          shift.selected = peaks;
          shift.integral = c.integral;
          shift.range = c.range;
          shift.nH = c.nH;
          shift.optPeaks = c.optPeaks;
        });
        metabolite.meanIntegral = finalCandidates[0].meanIntegral;
      }
    } else {
      metabolite.meanIntegral = metabolite.signals[0].integral;
    }
    toExport[ps.name] = metabolite;
    parentPort.postMessage("metaboites");
    // parentPort.postMessage(JSON.stringify(toExport));
  });
  //   fs.appendFileSync(`${subFix + String(index)}.json`, `${JSON.stringify(toExport)},`);
  // if (Object.keys(toExport).length < 3) {
  //   if (i === samples.length - 1) {
  //     fs.appendFileSync(filename, "]");
  //   }
  // }

  // if (first) {
  //   first = false;
  //   fs.appendFileSync(filename, `[${JSON.stringify(toExport)},`);
  // } else if (i === samples.length - 1) {
  //   fs.appendFileSync(filename, `${JSON.stringify(toExport)}]`);
  // } else {
  //   fs.appendFileSync(filename, `${JSON.stringify(toExport)},`);
  // }
  parentPort.postMessage("sale");
}
process.exit();

