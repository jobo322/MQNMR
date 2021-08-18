'use strict';

const path = require('path');
const fs = require('fs');

const Random = require('random-js');
const Pool = require('worker-threads-pool');

const pathInfo = path.resolve(path.join('data', 'infoLocalHMDB.json'));

//cosas importadas en worker
const converter = require('jcampconverter');

const utils = require('../../utils');

const byMetabo = require('../byMetabo/index');

const { getCombinationsScored } = utils;

const messenger = console.log;

//termina lo importado en worker

const defaultOptions = {
  thresholdFactor: 1,
  minMaxRatio: 0.0001,
  broadRatio: 0.0025,
  smoothY: false,
  widthFactor: 4,
  realTop: false,
  functionName: 'pseudovoigt',
  broadWidth: 0.25,
  sgOptions: { windowSize: 15, polynomial: 4 },
  optimizationOptions: { maxIterations: 500, shape: { kind: 'pseudovoigt' } },
};

let info = JSON.parse(fs.readFileSync(pathInfo, 'utf-8'));

const maxThreads = 1;
const maxWaiting = 4000;
const field = 600.89;
let subFix = 'newByWorkers3';

let sqrtPI = Math.sqrt(Math.PI);

let { peaksToSearch, pathToData = '', selectIt = [], existListAll = [] } = info;

let listSamples = fs.readdirSync(pathToData);

let toGet = [
  //   'airwave_13_470_1',
  //   'airwave_13_430_1',
  //   'airwave_22_330_1',
  'airwave_22_270_1',
  // 'airwave_37_370_1',
  // 'airwave_32_780_1',
  // 'airwave_33_40_1',
  // 'airwave_30_210_1',
];
console.log('listSamples.length before filter', listSamples.length);
listSamples = excludeIt(
  [
    { include: true, list: toGet },
    // { include: true, list: selectIt },
    // { include: true, list: existListAll },
  ],
  listSamples,
);

console.log('listSamples.length after filter', listSamples.length);

let samples = [];
let indexes = [];
let nbSamples = listSamples.length;
let r = new Random.Random(Random.MersenneTwister19937.seed(0x12345678));
for (let i = 0; i < nbSamples; ) {
  let index = r.integer(0, nbSamples - 1);
  if (indexes.includes(index)) continue;
  samples.push(listSamples[index]);
  indexes.push(index);
  i++;
}

let batchSize = Math.floor(samples.length / maxThreads);

let list = new Array(maxThreads).fill(0);
let diff = samples.length - batchSize * maxThreads;
list.forEach((_e, i, arr) => (arr[i] = samples.splice(0, batchSize)));
for (let i = 0; i < diff; i++) {
  list[i] = list[i].concat(samples.splice(0, 1));
}
messenger('first list', list[0].length);

let toSearch = [
  'eretic',
  // 'creatinine',
  'citrate',
];

for (let k = 0; k < maxThreads; k++) {
    let index = k;
    let samples = list[k];
    let use2D = false;
  let filteredPeaksToSearch = peaksToSearch.filter((e) => {
    return toSearch.includes(e.name);
  });

  let filename = `${subFix + index}.json`;

  messenger(`sample length ${samples.length}`);
  messenger(`name ${filename}`);

  //loop over partial list of samples assigned to the worker
  for (let i = 0; i < samples.length; i++) {
    messenger(`--------------- ${String(index)} - ${i}`);

    let sample = samples[i];
    let entry = sample.replace(/\.[a-z]*/g, '');
    messenger(entry);

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

    let toExport = { sampleid: entry };
    for (let ii = 0; ii < filteredPeaksToSearch.length; ii++) {
      let ps = filteredPeaksToSearch[ii];
      messenger(`ps.name ${ps.name}`);
      let { metabolite, toCombine, intPattern } =
        ps.name !== 'eretic'
          ? byMetabo.general(ps, xy, {
              entry,
              field,
              toExport,
              use2D,
              peaksToSearch,
              messenger,
              defaultOptions,
              sqrtPI,
              zMatrix,
            })
          : byMetabo.eretic(ps, xy, {
              field,
              peaksToSearch,
              messenger,
              defaultOptions,
              sqrtPI,
            });
      // messenger(`toCombine ${JSON.stringify(toCombine)}`)
      if (ps.name !== 'eretic') {
        // messenger(
        //   'toCombine length ' + JSON.stringify(toCombine.length),
        // );
        // messenger(
        //   'toSearch length ' + JSON.stringify(ps.toSearch.length),
        // );
        if (toCombine.length !== ps.toSearch.length) toCombine = [];
        let eretic = toExport.eretic;
        // messenger(eretic)
        let finalCandidates = getCombinationsScored(toCombine, {
          sqrtPI,
          eretic: eretic.meanIntegral,
          intPattern,
          messenger,
        });
        if (finalCandidates.length === 0) {
          messenger(`sin candidatos ${sample}:${ps.name} worker:${index}`);
        } else {
          finalCandidates.sort((a, b) => b.score - a.score);
          // messenger('finalCandidates');
          // messenger(
          //   JSON.stringify(finalCandidates.map((a) => a.score)),
          // );
          // messenger(
          //   JSON.stringify(finalCandidates.map((a) => a.IntegralScore)),
          // );
          // messenger(
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
      messenger(JSON.stringify(toExport));
    }
    // fs.appendFileSync(filename, `${JSON.stringify(toExport)},`);
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
}

function excludeIt(listOfList, peaks) {
  peaks = peaks.filter((sample, i, arr) => {
    let name = sample.sampleID ? sample.sampleID : sample;
    let entry = name.toLowerCase().replace(/\.[a-z]*/g, '');
    entry = entry.replace(/-/g, '_');
    return !listOfList.some((e) => {
      let list = e.list;
      return e.include ? !list.includes(entry) : list.includes(entry);
    });
  });
  return peaks;
}
