'use strict';

const fs = require('fs');
const path = require('path');

const sqrtPI = Math.sqrt(Math.PI);

// ----------- TO WRITE RESULT ------------//
const filename = 'test-creatine.json'; //'searchAll_oldVersion131.json';
const pathInfo = path.resolve(path.join('data', 'infoLocalHMDB.js'));
const pathToData = path.resolve('optimizedData/optimizeAllUpdated.json'); //optimize131TrainingFiles.json'); //optimizeSet1160Updated.json');//'optimizeSet131.json')//'optimizeAllUpdated.json');//

let info = require(pathInfo);
let optimizedPeaks = JSON.parse(fs.readFileSync(pathToData, 'utf8'));

let {
  peaksToSearch,
  existList3 = [],
  existListAll = [],
  selectIt = [],
  exclusionListEretic = [],
  exclusionListOrganic = [],
} = info;

let toGet = [
  //   'airwave_13_470_1',
  //   'airwave_13_430_1',
  //   'airwave_22_330_1',
  'airwave_22_770_1',
  // 'airwave_37_370_1',
  // 'airwave_32_780_1',
  // 'airwave_33_40_1',
  'airwave_30_210_1',
];
optimizedPeaks = excludeIt(
  [
    // { include: true, list: ['airwave_27_510_1'] },
    // { include: true, list: selectIt },
    // { include: true, list: toGet },
    // { include: true, list: existList3 },
    { include: true, list: existListAll },
    // { include: false, list: exclusionListEretic },
    // { include: false, list: exclusionListOrganic }
  ],
  optimizedPeaks,
);
// console.log(optimizedPeaks.length);
// return
// fs.writeFileSync('optimizeSet131.json', JSON.stringify(optimizedPeaks));

let toSearch = [
  'eretic',
  'creatinine',
  'citrate',
  'dimethylamine',
  'glycine',
  'creatine',
  'tartrate',
  'formate',
  'alanine',
  'sarcosine',
  'hippurate',
  'ethanol',
  'paracetamol',
  'guanidinoacetate',
  'glucose',
  'betaine',
  'lactate',
  'taurine',
  'dimethylglycine',
  'glycolate',
  'succinate',
  'trigonelline',
  'acetate',
  'trimethylamine',
  'pantothenic_acid',
];

toSearch = [
  'eretic',
  'creatinine',
  'citrate',
  'dimethylamine',
  'glycine',
  'formate',
  'creatine',
  'succinate',
  'alanine',
];

peaksToSearch = peaksToSearch.filter((e) => {
  return toSearch.includes(e.name);
});

let first = true;
optimizedPeaks.forEach((op, i) => {
  let toExport = { sampleid: op.sampleID };
  console.log(op.sampleID);
  const optData = op.peaks;
  const field = 600; // it should be inside of optimized peaks
  let assignedPeaks = [];
  peaksToSearch.forEach((ps, j) => {
    let toCombine = [];
    console.log(ps.name);
    let intPattern = [];
    let selectedPeaks = [];
    let integral = 0;
    let metabolite = { signals: [] };
    ps.toSearch.forEach((signal, k) => {
      let delta = signal.delta;
      let shift = { delta, nH: signal.integral, selected: [], integral: -0.1 }; // initialize partial result
      let { from, to } = signal.range;
      if (from > to) [from, to] = [to, from];

      intPattern.push(signal.integral); // part of checkIntegral filter

      let optPeaks = optData[ps.name].find((p) => {
        return p.delta.some((e) => e === delta);
      });

      if (!optPeaks) {
        console.log('no tiene optPeaks -------- ', op.sampleID);
        // console.log('delta', delta, from, to )
        return;
      }

      let peaks = optPeaks.peaks.filter((e) => e.x < to && e.x > from);
      if (!peaks.length) {
        console.log(peaks);
        console.log('sin picos ', toExport.sampleid);
        console.log('delta', delta, from, to);
      }

      if (ps.name !== 'eretic') {
        if (signal.getNoise) {
          let noiseLevel = getNoiseLevel(peaks.map((e) => e.y));
          peaks = peaks.filter((e) => e.y > noiseLevel);
        }
        let pattern = signal.pattern;
        let jCoupling = signal.jCoupling.map((j) => j / field);
        let candidates = [];
        let options = { delta, range: signal.range, nH: signal.integral };
        peaks.forEach((peak, pp, arr) => {
          let cand = getCandidates(
            arr,
            jCoupling,
            pattern,
            [{ indexs: [pp], score: 0 }],
            options,
          ); // generate combinations from J
          if (cand !== null) candidates = candidates.concat(cand);
        });
        console.log('candidates', candidates);
        // mover dentro de la funcion getCandidates
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
            console.log('---------\n----------\n');
            let prediction = singletPredictor(toExport, ps.name);
            console.log(prediction);
            console.log(candidates.map((e) => getDelta(e.peaks)));
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
        shift.optPeaks = peaks;
      } else {
        // if eretic peak is assigned
        // console.log(optPeaks.peaks)
        selectedPeaks = ps.callback(peaks);
        integral =
          selectedPeaks.reduce((a, b) => {
            let peak = b;
            return (
              peak.y * peak.width * sqrtPI * (1 - peak.xL + peak.xL * sqrtPI) +
              a
            );
          }, 0) / ps.peaks[0].integral;
        // let deltaS =
        //   selectedPeaks.reduce((a, b) => a + b.x, 0) / selectedPeaks.length;
        shift = Object.assign({}, shift, {
          integral,
          selected: selectedPeaks,
          optPeaks: peaks,
        }); // @TODO
      }
      metabolite.signals.push(shift);
    });
    // return;
    if (ps.name !== 'eretic') {
      if (toCombine.length !== ps.toSearch.length) toCombine = [];
      let eretic = toExport.eretic;
      let finalCandidates = checkIntegrals(optData, ps.peaks, toCombine, {
        sqrtPI,
        eretic: eretic.meanIntegral,
        intPattern,
      });
      if (finalCandidates.length === 0) {
        console.log(`sin candidatos ${op.sampleID}  ${i}`);
      } else {
        selectedPeaks = [];
        finalCandidates.sort((a, b) => b.score - a.score);

        finalCandidates[0].signals.forEach((c, i, arr) => {
          let peaks = c.peaks;
          let delta = c.delta;
          let shift = metabolite.signals.find((e) => e.delta === delta);
          shift.selected = peaks;
          shift.integral = c.integral;
          shift.range = c.range;
          shift.nH = c.nH;
        });
        metabolite.meanIntegral = finalCandidates[0].meanIntegral;
      }
    } else {
      metabolite.meanIntegral = metabolite.signals[0].integral;
    }

    toExport[ps.name] = metabolite;
  });

  if (Object.keys(toExport).length < 3) {
    if (i === optimizedPeaks.length - 1) {
      fs.appendFileSync(filename, ']');
    }
  }
  if (first) {
    first = false;
    fs.appendFileSync(filename, `[${JSON.stringify(toExport)},`);
  } else if (i === optimizedPeaks.length - 1) {
    fs.appendFileSync(filename, `${JSON.stringify(toExport)}]`);
  } else {
    fs.appendFileSync(filename, `${JSON.stringify(toExport)},`);
  }
});

// ----- functions -------//

function getCandidates(peaks, jcp, pattern, candidates, options = {}) {
  let debug = false;
  if (candidates.length === 0) return null;
  if (
    pattern.length === 0 ||
    candidates.some((e) => e.indexs.length === pattern.length)
  ) {
    let { delta, range, nH } = options;
    return candidates.map((cand) => {
      let indexs = cand.indexs;
      let toExport = {
        peaks: indexs.map((index) => {
          return peaks[index];
        }),
        score: cand.score,
        delta,
        range,
        nH,
      };
      return toExport;
    });
  }
  let len = peaks.length;
  let newCandidates = [];
  candidates.forEach((candidate, i, arr) => {
    let { indexs, score } = candidate;
    let index = indexs[indexs.length - 1];
    let iPattern = indexs.length - 1;
    if (debug) console.log('ipattern ', iPattern);
    let maxDiff = jcp[iPattern] + 0.4;
    for (let j = index + 1; j < len; j++) {
      let c = Math.abs(peaks[index].x - peaks[j].x);
      let diff = Math.abs(c - jcp[iPattern]) / jcp[iPattern];

      if (debug) console.log(`index ${index} j: ${j} length: ${len}`);
      if (debug) console.log('ref / exp / diff ---> ', jcp[iPattern], c, diff);

      if (c > maxDiff) {
        if (debug) console.log(c, maxDiff);
        if (debug) console.log('entra break per diff J');
        break;
      }
      if (diff < 0.1) {
        if (debug) console.log('pasa diff J');
        let RIP = pattern[iPattern] / pattern[iPattern + 1];
        let RIC = peaks[index].y / peaks[j].y;

        if (debug) {
          console.log(
            'pattern,',
            pattern[iPattern],
            pattern[iPattern + 1],
            RIP,
          );
        }
        if (debug) {
          console.log('experimental ', peaks[index].y, peaks[j].y, RIC);
        }

        let diffRI = Math.abs(RIP - RIC) / RIP;

        if (debug) console.log(`diffRI -> ${diffRI} max: 0.60`);

        if (diffRI < 0.6) {
          score += 1 - diffRI;
          newCandidates.push({ indexs: indexs.concat([j]), score });
          if (debug) console.log('candidate added');
        }
      } else {
        if (debug) console.log('Descartado diff > max');
      }
    }
  });
  return getCandidates(peaks, jcp, pattern, newCandidates, options);
}

function getCombinations(arrayOfArrays, options = {}) {
  let { eretic, sqrtPI, intPattern } = options;

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
    let integrals = [];
    let mean = 0;
    let similarityPatternScore = 0;
    let deltaScore = 0;
    combo.forEach((c, i, arr) => {
      let peaks = c.peaks;
      similarityPatternScore += c.score;
      if (c.deltaScore) {
        deltaScore += c.deltaScore;
      }
      let integral =
        (peaks.reduce((a, b) => {
          let peak = b;
          return (
            peak.y * peak.width * sqrtPI * (1 - peak.xL + peak.xL * sqrtPI) + a
          );
        }, 0) /
          intPattern[i] /
          eretic) *
        10;
      arr[i].integral = integral;
      mean += integral;
      integrals.push(integral);
    });
    similarityPatternScore /= combo.length;
    if (deltaScore !== 0) deltaScore /= combo.length;
    // console.log('combo', combo)
    // combo.forEach(c => {
    //     console.log(c.peaks.map(ee => ([ee.y,ee.x])))
    // })
    // console.log(JSON.stringify(integrals))
    // console.log('\n');
    // if (integrals.some(e => e < 0.05)) continue;
    mean /= combo.length;
    let std = integrals.reduce((a, b) => Math.pow(b - mean, 2) + a, 0);
    std = Math.sqrt(std / combo.length);
    let score = 1 - std / mean;
    // console.log(integrals, mean, std, std / mean);
    if (std / mean > 0.05) continue;
    let result = {
      signals: combo,
      meanIntegral: mean,
      similarityPatternScore,
      deltaScore,
      IntegralScore: score,
      score: score * 0.1 + similarityPatternScore * 0.1 + deltaScore * 0.8,
    };
    // combo.meanIntegral = mean;
    // combo.similarityPatternScore = similarityPatternScore;
    // combo.IntegralScore = score;
    // combo.deltaScore = deltaScore;
    // combo.score = ;
    // console.log('pasa')
    // console.log('\n mean ', std)
    // console.log(integrals);
    // console.log(combo.map(e => e.score))
    combinations.push(result);
  }
  return combinations;
}

function excludeIt(listOfList, peaks) {
  peaks = peaks.filter((sample) => {
    let entry = sample.sampleID.toLowerCase().replace(/\.[a-z]*/g, '');
    entry = entry.replace(/-/g, '_');
    return !listOfList.some((e) => {
      let list = e.list;
      return e.include ? !list.includes(entry) : list.includes(entry);
    });
  });
  return peaks;
}

function getDelta(peaks) {
  let delta = peaks.reduce((a, b) => a + b.x, 0) / peaks.length;
  return delta;
}

function checkIntegrals(expPeaks, toAssignPeaks, toCombine, options) {
  let combinations = getCombinations(toCombine, options);
  return combinations;
}

function getNoiseLevel(y) {
  let mean = 0;

  let stddev = 0;
  let length = y.length;
  for (let i = 0; i < length; ++i) {
    mean += y[i];
  }
  mean /= length;
  let averageDeviations = new Array(length);
  for (let i = 0; i < length; ++i) {
    averageDeviations[i] = Math.abs(y[i] - mean);
  }
  averageDeviations.sort((a, b) => a - b);
  if (length % 2 === 1) {
    stddev = averageDeviations[(length - 1) / 2] / 0.6745;
  } else {
    stddev =
      (0.5 *
        (averageDeviations[length / 2] + averageDeviations[length / 2 - 1])) /
      0.6745;
  }

  return stddev;
}
