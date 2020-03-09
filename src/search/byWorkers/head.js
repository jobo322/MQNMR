"use strict";

const path = require("path");
const fs = require("fs");
const { workerData } = require("worker_threads");

const Random = require("random-js");
const Pool = require("worker-threads-pool");

const pathInfo = path.resolve(path.join("data", "infoLocalHMDB.js"));

var info = require(pathInfo);

const maxThreads = 1;
const maxWaiting = 4000;
const field = 600.89;
var subFix = "optimizedPeaksByMetab";
var worker = path.resolve(
  path.join("src", "searchFromOptPeaks", "worker.js")
);

//This is temporal
// const pathToOptData = path.resolve('optimizeSet1160Updated.json');//'optimizeSet131.json')//'optimize131TrainingFiles.json'); //'optimizeAllUpdated.json');//
// var optimizedPeaks = JSON.parse(fs.readFileSync(pathToData, 'utf8'));
// optimizedPeaks = excludeIt([
//     { include: true, list: selectIt },
//   //   { include: true, list: existList3 },
//   ], optimizedPeaks);

// Finish temporal data

var pool = new Pool({ max: maxThreads, maxWaiting: maxWaiting });
var sqrtPI = Math.sqrt(Math.PI);

var { peaksToSearch, pathToData = '', selectIt = [], existListAll = [] } = info;
var listSamples = fs.readdirSync(pathToData);


let toGet = [
//   'airwave_13_470_1',
//   'airwave_13_430_1',
//   'airwave_22_330_1',
  'airwave_22_270_1',
  // 'airwave_37_370_1',
  // 'airwave_32_780_1',
  // 'airwave_33_40_1',
  'airwave_30_210_1'
];
console.log('listSamples.length before filter', listSamples.length);
listSamples = excludeIt([
  { include: true, list: toGet },
  // { include: true, list: selectIt },
  // { include: true, list: existListAll },
], listSamples);

console.log('listSamples.length after filter', listSamples.length);

var samples = [];
var indexes = [];
var nbSamples = listSamples.length;
var r = new Random.Random(Random.MersenneTwister19937.seed(0x12345678));
for (let i = 0; i < nbSamples;) {
  let index = r.integer(0, nbSamples - 1);
  if (indexes.includes(index)) continue;
  samples.push(listSamples[index]);
  indexes.push(index);
  i++;
}
var batchSize = Math.floor(samples.length / maxThreads);

var list = new Array(maxThreads).fill(0);
var diff = samples.length - batchSize * maxThreads;
list.forEach((_e, i, arr) => (arr[i] = samples.splice(0, batchSize)));
for (let i = 0; i < diff; i++) {
  list[i] = list[i].concat(samples.splice(0, 1));
}
console.log('samples.length',list[0].length);
var toSearch = [
  'eretic',
  //   'creatinine',
    // 'citrate',
  //   'glycine',
  //   'dimethylamine',
  //   'formate',
  // 'hippurate',
  // 'trigonelline',
  // 'tartrate',
  // 'succinate',
  // "alanine",
  // "taurine",
  "lactate",
  // 'acetate',
  // 'ethanol'
];
let peakList = [];
let rangeToOpt = peakList;

for (let i = 0; i < maxThreads; i++) {
    pool.acquire(
        worker,
        {
          workerData: {
            index: i,
            pathToData,
            sqrtPI,
            field,
            pathInfo,
            toSearch, 
            samples: list[i],
            rangeToOpt,
            subFix
          }
        },
        function(err, worker) {
          if (err) throw err;
          worker.on("message", data => console.log(data));
          console.log(`started worker ${i} (pool size: ${pool.size})`);
          worker.on("exit", function() {
            console.log(`worker ${i} exited (pool size: ${pool.size})`);
            // if (pool.size === 0) {
            //   const pathToData = path.resolve(".");
            //   var listSamples = fs.readdirSync(pathToData);
            //   var samples = listSamples.filter((s) =>
            //     s.match(/optimizedPeaksByMetab.*/)
            //   );
            //   var result = "[";
            //   samples.forEach((s) => {
            //     var peaks = fs.readFileSync(path.join(pathToData, s), "utf8");
            //     result = result.concat(peaks);
            //   });
    
            //   result = result.slice(0, result.length - 1).concat("]");
            //   fs.writeFileSync('optimizedPeaksByMetab.json', result);
            // }
          });
        }
      );
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