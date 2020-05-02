'use strict';

const path = require('path');
const { getPeaks } = require('../../../utilities/utils')
const utils = require('../../../utils.js');
const debug = false;
module.exports = function(ps, xy, options) {
  let {
    field,
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
    
    if (pattern.length > 1) {
      candidates.forEach((_e, i, arr) => {
        arr[i].score /= pattern.length - 1;
      });
    }
    
    if (candidates.length > 0) {
      if (
        ps.name === 'glycine' ||
        ps.name === 'formate' ||
        ps.name === 'alanine' ||
        ps.name === 'trigonelline' ||
        ps.name === 'tartrate' ||
        ps.name === 'creatine' ||
        ps.name === 'succinate'
      ) {
        let pathToPredictor = path.resolve(
          path.join('src/search/prediction', 'predictor.js'),
        );
        let { singletPredictor } = require(pathToPredictor);
        // console.log('---------\n----------\n');
        let prediction = singletPredictor(toExport, ps.name);
        // console.log(prediction);
        // console.log(candidates.map((e) => getDelta(e.peaks)));
        if (prediction !== null) {
          candidates.forEach((e, i, arr) => {
            let delta = getDelta(e.peaks);
            let score = (1 - Math.abs(delta - prediction[0])) * 10;
            arr[i].deltaScore = score;
          });
        }
      }
      toCombine.push(candidates);
    }
    
    shift.optPeaks = peaks; //It is not actually optimum values
    metabolite.signals.push(shift);
  });
  return { toCombine, metabolite, intPattern };
};
