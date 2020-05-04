'use strict';

const path = require('path');
const { getPeaks, optimizePeaks, runOptimization } = require('../../../utilities/utils')
const utils = require('../../../utils.js');
const debug = false;
module.exports = function(ps, xy, options) {
  let {
    field,
    toExport,
    delta,
    peaksToSearch,
    parentPort,
    defaultOptions,
    sqrtPI,
  } = options;
  let toCombine = [];
  let intPattern = [];
  let metabolite = { signals: [] };
  let index = 0;
  let optPeaks;
  
  ps.peaks.forEach((cluster) => {
    let delta = cluster.delta;
    let signal = ps.toSearch.find((e) => e.delta === delta);
    if (!signal) return;
    intPattern.push(signal.integral); // part of checkIntegral filter
    let shift = { delta, nH: signal.integral, selected: [], integral: -0.1 };
    let { peakList, data } = getPeaks(xy, cluster, options);

    peakList.forEach((p, pi, arr) => {
      arr[pi].index = index++;
    });
    
    if (peakList.length === 0) {
      parentPort.postMessage(
        `entry: ${entry} range: `,
      );
      return;
    }
    optPeaks = peakList; //optimizePeaks(peakList, data.x, data.y, options);
    let pattern = signal.pattern;
    let from = signal.range.from || signal.from; //@TODO: check it
    let to = signal.range.to || signal.to; //@TODO: check it
    
    if (from > to) [from, to] = [to, from];
    let peaks = optPeaks.filter((e) => e.x < to && e.x > from);

    if (signal.getNoise) {
      let noiseLevel = utils.getNoiseLevel(peaks.map((e) => e.y));
      peaks = peaks.filter((e) => e.y > noiseLevel);
    }
    
    let candidates = utils.getCandidatesByJ(peaks, signal, { field });

    // pasar de un rango amplio o muchos rangos pequeñós para realizar la optimizacion de parametros
    let optOptions = Object.assign({}, defaultOptions, cluster.gsdOptions);
        
    if (pattern.length > 1) {
      candidates.forEach((_e, i, arr) => {
        arr[i].score /= pattern.length - 1;
      });
    }
    
    if (candidates.length > 0) {
      if (
        // false
        ps.name === 'glycine' ||
        ps.name === 'formate' ||
        ps.name === 'alanine' ||
        ps.name === 'trigonelline' ||
        ps.name === 'tartrate' ||
        ps.name === 'creatine' ||
        ps.name === 'succinate'
      ) {
        let pathToPredictor = path.resolve(
          path.join('src/search/prediction', 'predictor'),
        );
        
        let { singletPredictor } = require(pathToPredictor);
        
        // console.log('---------\n----------\n');
        let prediction = singletPredictor(toExport, ps.name);
        //parentPort.postMessage(pathToPredictor)
        // console.log(prediction);
        //parentPort.postMessage(`prediction ${prediction}`)
        if (prediction !== null) {
          candidates.forEach((e, i, arr) => {
            let delta = utils.getDelta(e.peaks);
            let score = (1 - Math.abs(delta - prediction[0])) * 10;
            arr[i].deltaScore = score;
          });
        }
        candidates = candidates.filter(candidate => {
          let score = candidate.deltaScore;
          //parentPort.postMessage(`score ${Object.keys(candidate)}`)
          return score >= 9.8
        })
      }
    }

    if (pattern.length > 0) {
      candidates = runOptimization(xy, peaks, candidates, optOptions)
    }
    if (pattern.length === 0) {
      if (
        // false &&
        ps.name === 'glycine' ||
        ps.name === 'formate' ||
        ps.name === 'alanine' ||
        ps.name === 'trigonelline' ||
        ps.name === 'tartrate' ||
        ps.name === 'creatine' ||
        ps.name === 'succinate'
      ) {
        candidates = runOptimization(xy, peaks, candidates, optOptions)
        // for(let i = 0; i < candidates.length; i++) {
        //   let candPeaks = candidates[i].peaks;
        //   let first = candPeaks[0];
        //   let last = candPeaks[candPeaks.length - 1];
        //   from = first.x - first.width * 4;
        //   to = last.x + last.width * 4;
        //   let filteredPeaks = peaks.filter(peak => {
        //     let w3 = peak.width * 3;
        //     return (peak.x + w3) >= from && (peak.x - w3) <= to;
        //   })
        //   let optPeaks = optimizePeaks(filteredPeaks, xy.x, xy.y, optOptions);
        //   let peakIndex = candidates[i].peaks[0].index;
        //   candidates[i].peaks = [optPeaks.find(e => e.index === peakIndex)];
        //   candidates[i].optPeaks = optPeaks;
        // }
        //parentPort.postMessage(`candidate filtered ${candidates.length}`)
      } else {
        //parentPort.postMessage(`optimize ${JSON.stringify(candidates.map(e=>e.peaks[0].index))}`)
        let optPeaks = optimizePeaks(peaks, data.x, data.y, optOptions);
        for(let i = 0; i < candidates.length; i++) {
          let peakIndex = candidates[i].peaks[0].index;
          //parentPort.postMessage(`peakIndex ${peakIndex}`)
          candidates[i].peaks = [optPeaks.find(e => e.index === peakIndex)];
          candidates[i].optPeaks = optPeaks;
          //parentPort.postMessage(`peakIndexs ${candidates[i].peaks.length}`)
        }
        //parentPort.postMessage(`\n optimize2 ${JSON.stringify(optPeaks.map(e=>e.index))}`);
      }
    }
    if (candidates.length > 0) toCombine.push(candidates);
    metabolite.signals.push(shift);
  });
  return { toCombine, metabolite, intPattern };
};

