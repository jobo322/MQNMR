'use strict';

const MLR = require('ml-regression-multivariate-linear');
const { parentPort, workerData } = require('worker_threads');
// El orden de entrada
const models = {
  creatine: {
    name: 'multivariateLinearRegression',
    weights: [
      [0.28485266864299774],
      [0.04659555200487375],
      [0.0232930404599756],
      [-0.06136702001094818],
      [-0.29100191593170166],
      [0.6449560225009918],
      [-2.5835208892822266],
    ],
    inputs: 6,
    outputs: 1,
    intercept: true,
    summary: {
      regressionStatistics: {
        standardError: 0.0002592204235578109,
        observations: 1,
      },
      variables: [
        {
          label: 'X Variable 1',
          coefficients: [0.28485266864299774],
          standardError: 0.041030063112751515,
          tStat: 6.942535473567671,
        },
        {
          label: 'X Variable 2',
          coefficients: [0.04659555200487375],
          standardError: 0.03610487182155064,
          tStat: 1.29056134682249,
        },
        {
          label: 'X Variable 3',
          coefficients: [0.0232930404599756],
          standardError: 0.012375782208051122,
          tStat: 1.8821469276359928,
        },
        {
          label: 'X Variable 4',
          coefficients: [-0.06136702001094818],
          standardError: 0.019223226618653794,
          tStat: -3.19233712572472,
        },
        {
          label: 'X Variable 5',
          coefficients: [-0.29100191593170166],
          standardError: 0.09008674841954985,
          tStat: -3.2302410846982124,
        },
        {
          label: 'X Variable 6',
          coefficients: [0.6449560225009918],
          standardError: 0.12266857293653724,
          tStat: 5.257711955568772,
        },
        {
          label: 'Intercept',
          coefficients: [-2.5835208892822266],
          standardError: 0.8698697763212632,
          tStat: -2.9700087985676515,
        },
      ],
    },
  },
  glycine: {
    name: 'multivariateLinearRegression',
    weights: [
      [-0.15635141637176275],
      [-0.03511440381407738],
      [0.03422107093501836],
      [-0.006793254520744085],
      [0.9111678726039827],
      [1.7561778724193573],
    ],
    inputs: 5,
    outputs: 1,
    intercept: true,
    summary: {
      regressionStatistics: {
        standardError: 0.0004777815558630025,
        observations: 1,
      },
      variables: [
        {
          label: 'X Variable 1',
          coefficients: [-0.15635141637176275],
          standardError: 0.11022465687215362,
          tStat: -1.4184795018514798,
        },
        {
          label: 'X Variable 2',
          coefficients: [-0.03511440381407738],
          standardError: 0.17689745934390158,
          tStat: -0.19850145923131893,
        },
        {
          label: 'X Variable 3',
          coefficients: [0.03422107093501836],
          standardError: 0.021094560207705424,
          tStat: 1.6222699405943568,
        },
        {
          label: 'X Variable 4',
          coefficients: [-0.006793254520744085],
          standardError: 0.03786995207614271,
          tStat: -0.1793837633352511,
        },
        {
          label: 'X Variable 5',
          coefficients: [0.9111678726039827],
          standardError: 0.05811636998412149,
          tStat: 15.678334225846028,
        },
        {
          label: 'Intercept',
          coefficients: [1.7561778724193573],
          standardError: 0.19792628418232816,
          tStat: 8.8728886093854,
        },
      ],
    },
  },
  formate2: {
    name: 'multivariateLinearRegression',
    weights: [
      [-0.43388133496046066],
      [1.2181273754686117],
      [0.031324637588113546],
      [0.012073548510670662],
      [0.9743838608264923],
    ],
    inputs: 4,
    outputs: 1,
    intercept: true,
    summary: {
      regressionStatistics: {
        standardError: 0.0007934667591700934,
        observations: 1,
      },
      variables: [
        {
          label: 'X Variable 1',
          coefficients: [-0.43388133496046066],
          standardError: 0.10639040461909835,
          tStat: -4.078199876331458,
        },
        {
          label: 'X Variable 2',
          coefficients: [1.2181273754686117],
          standardError: 0.17517927786791976,
          tStat: 6.953604274970499,
        },
        {
          label: 'X Variable 3',
          coefficients: [0.031324637588113546],
          standardError: 0.017422741438173676,
          tStat: 1.797916688327846,
        },
        {
          label: 'X Variable 4',
          coefficients: [0.012073548510670662],
          standardError: 0.038568835047522,
          tStat: 0.31303897293746163,
        },
        {
          label: 'Intercept',
          coefficients: [0.9743838608264923],
          standardError: 0.2026132162903806,
          tStat: 4.809083428348661,
        },
      ],
    },
  },
  formate3: {
    name: 'multivariateLinearRegression',
    weights: [
      [-0.5055788224563003],
      [0.5382648333907127],
      [-0.01448969915509224],
      [0.04933098331093788],
      [-0.30990108847618103],
      [0.8846853673458099],
      [6.473986655473709],
    ],
    inputs: 6,
    outputs: 1,
    intercept: true,
    summary: {
      regressionStatistics: {
        standardError: 0.0002946838329429175,
        observations: 1,
      },
      variables: [
        {
          label: 'X Variable 1',
          coefficients: [-0.5055788224563003],
          standardError: 0.04823222570834511,
          tStat: -10.482178979537021,
        },
        {
          label: 'X Variable 2',
          coefficients: [0.5382648333907127],
          standardError: 0.08032238245659507,
          tStat: 6.701305625260586,
        },
        {
          label: 'X Variable 3',
          coefficients: [-0.01448969915509224],
          standardError: 0.00937084444691468,
          tStat: -1.5462533005617145,
        },
        {
          label: 'X Variable 4',
          coefficients: [0.04933098331093788],
          standardError: 0.015519936610339903,
          tStat: 3.1785557215531357,
        },
        {
          label: 'X Variable 5',
          coefficients: [-0.30990108847618103],
          standardError: 0.06026896585912651,
          tStat: -5.141967910990012,
        },
        {
          label: 'X Variable 6',
          coefficients: [0.8846853673458099],
          standardError: 0.05900160450684264,
          tStat: 14.994259473794642,
        },
        {
          label: 'Intercept',
          coefficients: [6.473986655473709],
          standardError: 0.13727179274328571,
          tStat: 47.161813261817166,
        },
      ],
    },
  },
  formate: {
    name: 'multivariateLinearRegression',
    weights: [
      [-0.43405030854046345],
      [0.424820576707134],
      [-0.0063019925728440285],
      [0.03086139727383852],
      [-0.0838621212169528],
      [0.682109571993351],
      [6.662017375230789],
    ],
    inputs: 6,
    outputs: 1,
    intercept: true,
    summary: {
      regressionStatistics: {
        standardError: 0.0004097112087502757,
        observations: 1,
      },
      variables: [
        {
          label: 'X Variable 1',
          coefficients: [-0.43405030854046345],
          standardError: 0.020859260638700066,
          tStat: -20.808518387040642,
        },
        {
          label: 'X Variable 2',
          coefficients: [0.424820576707134],
          standardError: 0.03046862309157686,
          tStat: 13.942887259141582,
        },
        {
          label: 'X Variable 3',
          coefficients: [-0.0063019925728440285],
          standardError: 0.0035115720058622834,
          tStat: -1.794635725060846,
        },
        {
          label: 'X Variable 4',
          coefficients: [0.03086139727383852],
          standardError: 0.006632403260958711,
          tStat: 4.653124374312776,
        },
        {
          label: 'X Variable 5',
          coefficients: [-0.0838621212169528],
          standardError: 0.01391309655682363,
          tStat: -6.027566967169714,
        },
        {
          label: 'X Variable 6',
          coefficients: [0.682109571993351],
          standardError: 0.016091108331971103,
          tStat: 42.39046546210127,
        },
        {
          label: 'Intercept',
          coefficients: [6.662017375230789],
          standardError: 0.0470132139650028,
          tStat: 141.70521037319583,
        },
      ],
    },
  },
  alanine: {
    name: 'multivariateLinearRegression',
    weights: [
      [-0.06235690042376518],
      [0.019673950038850307],
      [-0.013574167620390654],
      [0.024713167433219496],
      [-0.1635812595486641],
      [0.6436380073428154],
      [0.10917951166629791],
      [-1.1233350038528442],
    ],
    inputs: 7,
    outputs: 1,
    intercept: true,
    summary: {
      regressionStatistics: {
        standardError: 0.0001286776893145295,
        observations: 1,
      },
      variables: [
        {
          label: 'X Variable 1',
          coefficients: [-0.06235690042376518],
          standardError: 0.012085711866845372,
          tStat: -5.159555441233737,
        },
        {
          label: 'X Variable 2',
          coefficients: [0.019673950038850307],
          standardError: 0.013511036943706674,
          tStat: 1.4561391639162282,
        },
        {
          label: 'X Variable 3',
          coefficients: [-0.013574167620390654],
          standardError: 0.003713897025092564,
          tStat: -3.654966071670319,
        },
        {
          label: 'X Variable 4',
          coefficients: [0.024713167433219496],
          standardError: 0.006033638686573774,
          tStat: 4.095897801804398,
        },
        {
          label: 'X Variable 5',
          coefficients: [-0.1635812595486641],
          standardError: 0.026366585531533707,
          tStat: -6.204112373709726,
        },
        {
          label: 'X Variable 6',
          coefficients: [0.6436380073428154],
          standardError: 0.03928017222752305,
          tStat: 16.38582447181399,
        },
        {
          label: 'X Variable 7',
          coefficients: [0.10917951166629791],
          standardError: 0.03306277953164503,
          tStat: 3.3021879349798793,
        },
        {
          label: 'Intercept',
          coefficients: [-1.1233350038528442],
          standardError: 0.22337072584498227,
          tStat: -5.029016222261958,
        },
      ],
    },
  },
  trigonelline: {
    name: 'multivariateLinearRegression',
    weights: [
      [-0.7660428490489721],
      [1.9493951611220837],
      [-0.02106030867435038],
      [0.023250930476933718],
      [1.608876246958971],
    ],
    inputs: 4,
    outputs: 1,
    intercept: true,
    summary: {
      regressionStatistics: {
        standardError: 0.00035424106176396627,
        observations: 1,
      },
      variables: [
        {
          label: 'X Variable 1',
          coefficients: [-0.7660428490489721],
          standardError: 0.04433374514005019,
          tStat: -17.279001506167475,
        },
        {
          label: 'X Variable 2',
          coefficients: [1.9493951611220837],
          standardError: 0.06784605243837384,
          tStat: 28.732624685758466,
        },
        {
          label: 'X Variable 3',
          coefficients: [-0.02106030867435038],
          standardError: 0.007868738398166222,
          tStat: -2.676453023175659,
        },
        {
          label: 'X Variable 4',
          coefficients: [0.023250930476933718],
          standardError: 0.017178853292836208,
          tStat: 1.3534623109349004,
        },
        {
          label: 'Intercept',
          coefficients: [1.608876246958971],
          standardError: 0.07090377613572704,
          tStat: 22.690981138708203,
        },
      ],
    },
  },
  tartrate: {
    name: 'multivariateLinearRegression',
    weights: [
      [-0.45710545778274536],
      [0.5758174192160368],
      [0.25951624661684036],
      [1.0107949003577232],
      [-0.4825252704322338],
      [1.454115405678749],
    ],
    inputs: 5,
    outputs: 1,
    intercept: true,
    summary: {
      regressionStatistics: {
        standardError: 0.0006331215469030368,
        observations: 1,
      },
      variables: [
        {
          label: 'X Variable 1',
          coefficients: [-0.45710545778274536],
          standardError: 0.07838915641678443,
          tStat: -5.831233281199993,
        },
        {
          label: 'X Variable 2',
          coefficients: [0.5758174192160368],
          standardError: 0.10740513943171888,
          tStat: 5.361171935185687,
        },
        {
          label: 'X Variable 3',
          coefficients: [0.25951624661684036],
          standardError: 0.010378819747616128,
          tStat: 25.00440829762437,
        },
        {
          label: 'X Variable 4',
          coefficients: [1.0107949003577232],
          standardError: 0.09692981319033565,
          tStat: 10.428111507580045,
        },
        {
          label: 'X Variable 5',
          coefficients: [-0.4825252704322338],
          standardError: 0.07947148088710126,
          tStat: -6.071678356135311,
        },
        {
          label: 'Intercept',
          coefficients: [1.454115405678749],
          standardError: 0.19851822736443525,
          tStat: 7.324845808789724,
        },
      ],
    },
  },
  succinate: {
    name: 'multivariateLinearRegression',
    weights: [
      [0.6374383233487606],
      [-0.7709003407508135],
      [0.046518510673195124],
      [0.026243122294545174],
      [0.44427351653575897],
      [0.18457704782485962],
      [-0.6133911684155464],
      [3.5202205181121826],
    ],
    inputs: 7,
    outputs: 1,
    intercept: true,
    summary: {
      regressionStatistics: {
        standardError: 0.0005589507409700816,
        observations: 1,
      },
      variables: [
        {
          label: 'X Variable 1',
          coefficients: [0.6374383233487606],
          standardError: 0.09429625532767763,
          tStat: 6.759953734468829,
        },
        {
          label: 'X Variable 2',
          coefficients: [-0.7709003407508135],
          standardError: 0.12492895531384636,
          tStat: -6.170709895189299,
        },
        {
          label: 'X Variable 3',
          coefficients: [0.046518510673195124],
          standardError: 0.016382418408975516,
          tStat: 2.8395386756639542,
        },
        {
          label: 'X Variable 4',
          coefficients: [0.026243122294545174],
          standardError: 0.030012443711784808,
          tStat: 0.8744080470941605,
        },
        {
          label: 'X Variable 5',
          coefficients: [0.44427351653575897],
          standardError: 0.07795091927616976,
          tStat: 5.699400605678002,
        },
        {
          label: 'X Variable 6',
          coefficients: [0.18457704782485962],
          standardError: 0.09662270328912506,
          tStat: 1.9102865221287373,
        },
        {
          label: 'X Variable 7',
          coefficients: [-0.6133911684155464],
          standardError: 0.0695694566715041,
          tStat: -8.816960743446392,
        },
        {
          label: 'Intercept',
          coefficients: [3.5202205181121826],
          standardError: 0.8889582681495087,
          tStat: 3.9599390030310593,
        },
      ],
    },
  },
  succinateOld: {
    name: 'multivariateLinearRegression',
    weights: [
      [-0.21366441110149026],
      [0.08969045709818602],
      [0.08570078830234706],
      [0.03173782583326101],
      [0.21205786429345608],
      [1.9364686533808708],
    ],
    inputs: 5,
    outputs: 1,
    intercept: true,
    summary: {
      regressionStatistics: {
        standardError: 0.001034301971577298,
        observations: 1,
      },
      variables: [
        {
          label: 'X Variable 1',
          coefficients: [-0.21366441110149026],
          standardError: 0.10332584276301733,
          tStat: -2.067870006069436,
        },
        {
          label: 'X Variable 2',
          coefficients: [0.08969045709818602],
          standardError: 0.15805765152144166,
          tStat: 0.5674540665057196,
        },
        {
          label: 'X Variable 3',
          coefficients: [0.08570078830234706],
          standardError: 0.022115709383515245,
          tStat: 3.875109173130449,
        },
        {
          label: 'X Variable 4',
          coefficients: [0.03173782583326101],
          standardError: 0.041621047900190514,
          tStat: 0.762542690164121,
        },
        {
          label: 'X Variable 5',
          coefficients: [0.21205786429345608],
          standardError: 0.05327808739092012,
          tStat: 3.980207899309762,
        },
        {
          label: 'Intercept',
          coefficients: [1.9364686533808708],
          standardError: 0.2366612923042114,
          tStat: 8.182447727411532,
        },
      ],
    },
  },
};

const signalsReporter = {
  glycine: [
    { name: 'creatinine', delta: 4.05 },
    { name: 'creatinine', delta: 3.03 },
    { name: 'citrate', delta: 2.66 },
    { name: 'citrate', delta: 2.54 },
    { name: 'dimethylamine', delta: 2.7 },
  ],
  formate: [
    { name: 'creatinine', delta: 4.05 },
    { name: 'creatinine', delta: 3.03 },
    { name: 'citrate', delta: 2.66 },
    { name: 'citrate', delta: 2.54 },
    { name: 'dimethylamine', delta: 2.7 },
    { name: 'glycine', delta: 3.54 },
  ],
  alanine: [
    { name: 'creatinine', delta: 4.05 },
    { name: 'creatinine', delta: 3.03 },
    { name: 'citrate', delta: 2.66 },
    { name: 'citrate', delta: 2.54 },
    { name: 'dimethylamine', delta: 2.7 },
    { name: 'glycine', delta: 3.54 },
    { name: 'formate', delta: 8.44 },
  ],
  trigonelline: [
    { name: 'creatinine', delta: 4.05 },
    { name: 'creatinine', delta: 3.03 },
    { name: 'citrate', delta: 2.66 },
    { name: 'citrate', delta: 2.54 },
  ],
  tartrate: [
    { name: 'creatinine', delta: 4.05 },
    { name: 'creatinine', delta: 3.03 },
    { name: 'citrate', delta: 2.66 },
    { name: 'glycine', delta: 3.54 },
    { name: 'dimethylamine', delta: 2.7 },
  ],
  succinate: [
    { name: 'creatinine', delta: 4.05 },
    { name: 'creatinine', delta: 3.03 },
    { name: 'citrate', delta: 2.66 },
    { name: 'citrate', delta: 2.54 },
    { name: 'glycine', delta: 3.54 },
    { name: 'formate', delta: 8.44 },
    { name: 'hippurate', delta: 7.62 },
  ],
  succinateOld: [
    { name: 'creatinine', delta: 4.05 },
    { name: 'creatinine', delta: 3.03 },
    { name: 'citrate', delta: 2.66 },
    { name: 'citrate', delta: 2.54 },
    { name: 'glycine', delta: 3.54 },
  ],
  creatine: [
    { name: 'creatinine', delta: 4.05 },
    { name: 'creatinine', delta: 3.03 },
    { name: 'citrate', delta: 2.66 },
    { name: 'citrate', delta: 2.54 },
    { name: 'glycine', delta: 3.54 },
    { name: 'formate', delta: 8.44 },
  ],
};
module.exports.signalsReporter = signalsReporter;
module.exports.singletPredictor = function(state, name) {
  
  let reporters = signalsReporter[name];
  let model = models[name];
  
  let xValues = reporters.map((r) => {
    return getDelta(state, r);
  });
 
  if (xValues.some((e) => e === null)) {
    return null;
  }
  let mlrModel = new MLR(true, model);
  let predicted = mlrModel._predict(xValues);
  return predicted;
};

function getDelta(state, reporter) {
  let metabolite = state[reporter.name];
  if (!metabolite) return null;
  let signal = metabolite.signals.find((me) => me.delta === reporter.delta);
  
  if (!signal) return null;
  // parentPort.postMessage(`reporter ${reporter.name}`)
  // parentPort.postMessage(`reporter ${signal.selected.length}`)
  if (signal.selected.length === 0) return null;
  let delta =
    signal.selected.reduce((a, b) => a + b.x, 0) / signal.selected.length;
  return [delta];
}
