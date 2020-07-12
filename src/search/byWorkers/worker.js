'use strict';

const { parentPort, workerData } = require('worker_threads');
const path = require('path');
const fs = require('fs');
const { writeFileSync } = require('fs');
const { join } = require('path');

const converter = require('jcampconverter');

const utils = require('../../utils');

const byMetabo = require('./byMetabo/index');

const { getCombinationsScored } = utils;

const defaultOptions = {
  thresholdFactor: 1,
  minMaxRatio: 0.0001,
  broadRatio: 0.0025,
  smoothY: false,
  widthFactor: 4,
  realTop: false,
  fnType: 'voigt',
  broadWidth: 0.25,
  sgOptions: { windowSize: 15, polynomial: 4 },
};

let {
  index,
  samples,
  pathToData,
  pathInfo,
  toSearch,
  sqrtPI,
  field,
  subFix,
  use2D,
} = workerData;

let info = require(pathInfo);
const { peaksToSearch } = info;
let filteredPeaksToSearch = peaksToSearch.filter((e) => {
  return toSearch.includes(e.name);
});

let filename = `${subFix + index}.json`;
parentPort.postMessage(`sample length ${samples.length}`);
parentPort.postMessage(`name ${filename}`);

//loop over partial list of samples assigned to the worker
for (let i = 0; i < samples.length; i++) {
  parentPort.postMessage(`--------------- ${String(index)} - ${i}`);

  let sample = samples[i];
  let entry = sample.replace(/\.[a-z]*/g, '');
  parentPort.postMessage(entry);

  let pathToJcamp = path.join(pathToData, sample);
  let jcamp = fs.readFileSync(pathToJcamp, 'utf8');
  let spectrum = converter.convert(jcamp, { xy: true });
  let xy = spectrum.spectra[0].data[0];
  if (xy.x[0] > xy.x[1]) {
    xy.x = xy.x.reverse();
    xy.y = xy.y.reverse();
  }

  let zMatrix;
  if (use2D) {
    let number =
      parseInt(sample.replace(/\w+\-\w+\-(\w+)\-\w+/g, '$1'), 10) + 1;
    pathToJcamp = path.join(
      pathToData,
      sample.replace(/(\w+\-\w+\-)(\w+)(\-\w+)/g, `$1${number}$3`),
    );
    let jcamp2D = fs.readFileSync(pathToJcamp, 'utf8');
    spectrum = converter.convert(jcamp2D, {
      noContour: true,
      xy: true,
      keepRecordsRegExp: /.*/,
    });
    zMatrix = spectrum.minMax;

    // writeFileSync(
    //   join(sample + '.json'),
    //   JSON.stringify(zMatrix, null, 2),
    //   'utf8',
    // );
  }
  // parentPort.postMessage(zMatrix);

  let toExport = { sampleid: entry };
  for (let ii = 0; ii < filteredPeaksToSearch.length; ii++) {
    let ps = filteredPeaksToSearch[ii];
    parentPort.postMessage(ps.name);
    let { metabolite, toCombine, intPattern } =
      ps.name !== 'eretic'
        ? byMetabo.general(ps, xy, {
            entry,
            field,
            toExport,
            use2D,
            peaksToSearch,
            parentPort,
            defaultOptions,
            sqrtPI,
            zMatrix,
          })
        : byMetabo.eretic(ps, xy, {
            field,
            peaksToSearch,
            parentPort,
            defaultOptions,
            sqrtPI,
          });
    // parentPort.postMessage(`toCombine ${JSON.stringify(toCombine)}`)
    if (ps.name !== 'eretic') {
      // parentPort.postMessage(
      //   'toCombine length ' + JSON.stringify(toCombine.length),
      // );
      // parentPort.postMessage(
      //   'toSearch length ' + JSON.stringify(ps.toSearch.length),
      // );
      if (toCombine.length !== ps.toSearch.length) toCombine = [];
      let eretic = toExport.eretic;
      // parentPort.postMessage(eretic)
      let finalCandidates = getCombinationsScored(toCombine, {
        sqrtPI,
        eretic: eretic.meanIntegral,
        intPattern,
        parentPort,
      });
      if (finalCandidates.length === 0) {
        parentPort.postMessage(
          `sin candidatos ${sample}:${ps.name} worker:${index}`,
        );
      } else {
        finalCandidates.sort((a, b) => b.score - a.score);
        // parentPort.postMessage('finalCandidates');
        // parentPort.postMessage(
        //   JSON.stringify(finalCandidates.map((a) => a.score)),
        // );
        // parentPort.postMessage(
        //   JSON.stringify(finalCandidates.map((a) => a.IntegralScore)),
        // );
        // parentPort.postMessage(
        //   JSON.stringify(finalCandidates.map((a) => a.similarityPatternScore)),
        // );
        let indexCand = 0;
        indexCand = indexCand >= finalCandidates.length ? 0 : indexCand;
        finalCandidates[indexCand].signals.forEach((c) => {
          let peaks = c.peaks;
          let delta = c.delta;
          let shift = metabolite.signals.find((e) => e.delta === delta);
          shift.selected = peaks;
          shift.integral = c.integral;
          shift.range = c.range;
          shift.nH = c.nH;
          shift.optPeaks = c.optPeaks;
        });
        metabolite.meanIntegral = finalCandidates[indexCand].meanIntegral;
      }
    } else {
      metabolite.meanIntegral = metabolite.signals[0].integral;
    }
    toExport[ps.name] = metabolite;
    // parentPort.postMessage('metabolites');
    // parentPort.postMessage(JSON.stringify(toExport));
  }
  fs.appendFileSync(filename, `${JSON.stringify(toExport)},`);
  // if (Object.keys(toExport).length < 3) {
  //   if (i === samples.length - 1) {
  //     fs.appendFileSync(filename, ']');
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
}
parentPort.postMessage('sale');
process.exit();
