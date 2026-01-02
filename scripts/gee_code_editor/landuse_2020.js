/*
========================================================
北海道紋別市・湧別町 土地利用分類 2020年版
Google Earth Engine Code Editor (JavaScript)

使い方:
1. https://code.earthengine.google.com/ を開く
2. このコードをコピー&ペースト
3. 「Run」ボタンをクリック
========================================================
*/

// ============================================
// 1. 設定パラメータ
// ============================================

// 対象地域（紋別市・湧別町）
var STUDY_AREA = ee.Geometry.Rectangle([143.0, 43.8, 144.5, 44.6]);

// 分析対象年
var TARGET_YEAR = 2020;
var START_DATE = TARGET_YEAR + '-05-01';
var END_DATE = TARGET_YEAR + '-09-30';

// 雲量閾値
var CLOUD_COVER_MAX = 20;

// Random Forest パラメータ
var RF_TREES = 100;

// 分類クラス
var CLASS_NAMES = [
  '牧草地',           // 1
  'トウモロコシ畑',   // 2
  '畑作農地',         // 3
  '森林',             // 4
  '市街地',           // 5
  '水域',             // 6
  '裸地'              // 7
];

var CLASS_PALETTE = [
  '90EE90',  // 1: 牧草地 - ライトグリーン
  'FFD700',  // 2: トウモロコシ - ゴールド
  'FFA500',  // 3: 畑作 - オレンジ
  '006400',  // 4: 森林 - ダークグリーン
  'FF0000',  // 5: 市街地 - 赤
  '0000FF',  // 6: 水域 - 青
  '808080'   // 7: 裸地 - グレー
];

// ============================================
// 2. Landsat 8 データ取得・前処理
// ============================================

// 雲マスキング関数
function maskL8Clouds(image) {
  var qa = image.select('QA_PIXEL');
  var cloudBitMask = 1 << 3;
  var cloudShadowBitMask = 1 << 4;

  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0));

  return image.updateMask(mask);
}

// スケールファクター適用
function applyScaleFactors(image) {
  var opticalBands = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'])
    .multiply(0.0000275).add(-0.2);

  return image.addBands(opticalBands, null, true);
}

// Landsat 8 コレクション取得
var collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(STUDY_AREA)
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_COVER_MAX))
  .map(maskL8Clouds)
  .map(applyScaleFactors);

print('取得画像数:', collection.size());

// ============================================
// 3. 植生指数の計算
// ============================================

function addIndices(image) {
  // NDVI
  var ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

  // EVI
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('SR_B5'),
      'RED': image.select('SR_B4'),
      'BLUE': image.select('SR_B2')
    }).rename('EVI');

  // NDWI
  var ndwi = image.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');

  // LSWI
  var lswi = image.normalizedDifference(['SR_B5', 'SR_B6']).rename('LSWI');

  // NDBI
  var ndbi = image.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI');

  return image.addBands([ndvi, evi, ndwi, lswi, ndbi]);
}

// ============================================
// 4. 月別コンポジット作成
// ============================================

function createMonthlyComposite(month) {
  var start = ee.Date.fromYMD(TARGET_YEAR, month, 1);
  var end = start.advance(1, 'month');

  var monthlyCol = collection.filterDate(start, end);
  var composite = monthlyCol.median();
  composite = addIndices(composite);

  var bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7',
               'NDVI', 'EVI', 'NDWI', 'LSWI', 'NDBI'];
  var prefix = 'M' + month + '_';
  var renamedBands = bands.map(function(b) { return prefix + b; });

  return composite.select(bands).rename(renamedBands);
}

// 5-9月のコンポジット
var m5 = createMonthlyComposite(5);
var m6 = createMonthlyComposite(6);
var m7 = createMonthlyComposite(7);
var m8 = createMonthlyComposite(8);
var m9 = createMonthlyComposite(9);

var monthlyStack = m5.addBands(m6).addBands(m7).addBands(m8).addBands(m9);

// ============================================
// 5. 時系列統計量
// ============================================

var collectionWithIndices = collection.map(addIndices);

var ndviMax = collectionWithIndices.select('NDVI').max().rename('NDVI_max');
var ndviMin = collectionWithIndices.select('NDVI').min().rename('NDVI_min');
var ndviMean = collectionWithIndices.select('NDVI').mean().rename('NDVI_mean');
var ndviStd = collectionWithIndices.select('NDVI').reduce(ee.Reducer.stdDev()).rename('NDVI_std');

var eviMax = collectionWithIndices.select('EVI').max().rename('EVI_max');
var eviMean = collectionWithIndices.select('EVI').mean().rename('EVI_mean');

var ndwiMean = collectionWithIndices.select('NDWI').mean().rename('NDWI_mean');
var lswiMean = collectionWithIndices.select('LSWI').mean().rename('LSWI_mean');

var temporalStats = ndviMax.addBands(ndviMin).addBands(ndviMean).addBands(ndviStd)
  .addBands(eviMax).addBands(eviMean).addBands(ndwiMean).addBands(lswiMean);

// ============================================
// 6. フェノロジー特徴量（牧草地 vs トウモロコシ識別用）
// ============================================

// 7月、8月、9月のNDVI
var julyCol = collection.filterDate(TARGET_YEAR + '-07-01', TARGET_YEAR + '-07-31').map(addIndices);
var augCol = collection.filterDate(TARGET_YEAR + '-08-01', TARGET_YEAR + '-08-31').map(addIndices);
var septCol = collection.filterDate(TARGET_YEAR + '-09-01', TARGET_YEAR + '-09-30').map(addIndices);

var ndviJuly = julyCol.select('NDVI').median().rename('NDVI_july');
var ndviAug = augCol.select('NDVI').median().rename('NDVI_aug');
var ndviSept = septCol.select('NDVI').median().rename('NDVI_sept');

// NDVI差分（トウモロコシ識別に重要）
var ndviDiffAugJuly = ndviAug.subtract(ndviJuly).rename('NDVI_diff_aug_july');
var ndviDiffSeptAug = ndviSept.subtract(ndviAug).rename('NDVI_diff_sept_aug');

// 夏季最大NDVI
var summerMaxNdvi = julyCol.merge(augCol).select('NDVI').max().rename('NDVI_summer_max');

var phenology = ndviJuly.addBands(ndviAug).addBands(ndviSept)
  .addBands(ndviDiffAugJuly).addBands(ndviDiffSeptAug).addBands(summerMaxNdvi);

// ============================================
// 7. 地形データ
// ============================================

var dem = ee.Image('USGS/SRTMGL1_003');
var elevation = dem.select('elevation').rename('elevation');
var slope = ee.Terrain.slope(dem).rename('slope');
var aspect = ee.Terrain.aspect(dem).rename('aspect');
var aspectSin = aspect.multiply(Math.PI / 180).sin().rename('aspect_sin');
var aspectCos = aspect.multiply(Math.PI / 180).cos().rename('aspect_cos');

var terrain = elevation.addBands(slope).addBands(aspectSin).addBands(aspectCos);

// ============================================
// 8. 特徴量スタック
// ============================================

var featureStack = monthlyStack
  .addBands(temporalStats)
  .addBands(phenology)
  .addBands(terrain)
  .clip(STUDY_AREA);

print('特徴量バンド:', featureStack.bandNames());

// ============================================
// 9. 教師データの定義
// ============================================
/*
★★★ 教師データをここに追加 ★★★

方法1: 手動でポイントを追加（下のサンプルを編集）
方法2: Geometry Toolsでポイントを作成
  - 左側パネルの「Geometry Imports」で各クラスのポイントを作成
  - 作成したジオメトリを training_data にマージ

クラス:
  1 = 牧草地
  2 = トウモロコシ畑
  3 = 畑作農地
  4 = 森林
  5 = 市街地
  6 = 水域
  7 = 裸地
*/

// サンプル教師データ（実際の座標に置き換えてください）
var grassland = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([143.45, 44.18]), {class: 1}),
  ee.Feature(ee.Geometry.Point([143.50, 44.22]), {class: 1}),
  ee.Feature(ee.Geometry.Point([143.55, 44.25]), {class: 1}),
  ee.Feature(ee.Geometry.Point([143.48, 44.15]), {class: 1}),
  ee.Feature(ee.Geometry.Point([143.52, 44.20]), {class: 1}),
]);

var corn = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([143.65, 44.12]), {class: 2}),
  ee.Feature(ee.Geometry.Point([143.70, 44.15]), {class: 2}),
  ee.Feature(ee.Geometry.Point([143.68, 44.10]), {class: 2}),
  ee.Feature(ee.Geometry.Point([143.72, 44.18]), {class: 2}),
  ee.Feature(ee.Geometry.Point([143.67, 44.14]), {class: 2}),
]);

var otherCrop = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([143.60, 44.08]), {class: 3}),
  ee.Feature(ee.Geometry.Point([143.62, 44.12]), {class: 3}),
  ee.Feature(ee.Geometry.Point([143.58, 44.10]), {class: 3}),
]);

var forest = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([143.35, 44.40]), {class: 4}),
  ee.Feature(ee.Geometry.Point([143.40, 44.45]), {class: 4}),
  ee.Feature(ee.Geometry.Point([143.38, 44.42]), {class: 4}),
  ee.Feature(ee.Geometry.Point([143.42, 44.48]), {class: 4}),
  ee.Feature(ee.Geometry.Point([143.45, 44.50]), {class: 4}),
]);

var urban = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([143.35, 44.35]), {class: 5}),
  ee.Feature(ee.Geometry.Point([143.36, 44.34]), {class: 5}),
  ee.Feature(ee.Geometry.Point([143.34, 44.36]), {class: 5}),
]);

var water = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([143.90, 44.00]), {class: 6}),
  ee.Feature(ee.Geometry.Point([143.88, 44.02]), {class: 6}),
  ee.Feature(ee.Geometry.Point([143.92, 43.98]), {class: 6}),
]);

var bare = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([143.80, 44.05]), {class: 7}),
  ee.Feature(ee.Geometry.Point([143.82, 44.03]), {class: 7}),
]);

// 全教師データを統合
var trainingData = grassland.merge(corn).merge(otherCrop)
  .merge(forest).merge(urban).merge(water).merge(bare);

print('教師データ数:', trainingData.size());

// ============================================
// 10. Random Forest 分類
// ============================================

// 教師データにバンド値を抽出
var trainingSamples = featureStack.sampleRegions({
  collection: trainingData,
  properties: ['class'],
  scale: 30,
  tileScale: 8
});

print('訓練サンプル数:', trainingSamples.size());

// 分類器を訓練
var classifier = ee.Classifier.smileRandomForest(RF_TREES)
  .train({
    features: trainingSamples,
    classProperty: 'class',
    inputProperties: featureStack.bandNames()
  });

// 分類実行
var classified = featureStack.classify(classifier);

// ============================================
// 11. 精度評価
// ============================================

// データ分割
var withRandom = trainingData.randomColumn('random');
var trainSet = withRandom.filter(ee.Filter.lt('random', 0.7));
var testSet = withRandom.filter(ee.Filter.gte('random', 0.7));

// 訓練
var trainSamples = featureStack.sampleRegions({
  collection: trainSet,
  properties: ['class'],
  scale: 30,
  tileScale: 8
});

var classifierEval = ee.Classifier.smileRandomForest(RF_TREES)
  .train({
    features: trainSamples,
    classProperty: 'class',
    inputProperties: featureStack.bandNames()
  });

// テスト
var testSamples = featureStack.sampleRegions({
  collection: testSet,
  properties: ['class'],
  scale: 30,
  tileScale: 8
});

var validated = testSamples.classify(classifierEval);

// 混同行列
var confusionMatrix = validated.errorMatrix('class', 'classification');
print('混同行列:', confusionMatrix);
print('全体精度:', confusionMatrix.accuracy());
print('カッパ係数:', confusionMatrix.kappa());

// ============================================
// 12. 特徴量重要度
// ============================================

var importance = classifier.explain();
print('特徴量重要度:', importance);

// ============================================
// 13. 地図表示
// ============================================

// 地図の中心を設定
Map.setCenter(143.7, 44.2, 10);

// True Color
var composite = collection.median();
var visTC = {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0, max: 0.3};
Map.addLayer(composite.clip(STUDY_AREA), visTC, 'True Color');

// NDVI（8月）
var visNDVI = {min: 0.2, max: 0.9, palette: ['red', 'yellow', 'green', 'darkgreen']};
Map.addLayer(ndviAug.clip(STUDY_AREA), visNDVI, '8月 NDVI');

// 分類結果
var visClass = {min: 1, max: 7, palette: CLASS_PALETTE};
Map.addLayer(classified.clip(STUDY_AREA), visClass, '土地利用分類 2020');

// 教師データ
Map.addLayer(trainingData, {color: 'white'}, '教師データ');

// 対象地域の境界
Map.addLayer(ee.Image().paint(STUDY_AREA, 0, 2), {palette: 'blue'}, '対象地域');

// ============================================
// 14. 凡例パネル
// ============================================

var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

var legendTitle = ui.Label({
  value: '土地利用分類',
  style: {fontWeight: 'bold', fontSize: '16px', margin: '0 0 4px 0'}
});
legend.add(legendTitle);

for (var i = 0; i < CLASS_NAMES.length; i++) {
  var colorBox = ui.Label({
    style: {
      backgroundColor: '#' + CLASS_PALETTE[i],
      padding: '8px',
      margin: '0 0 4px 0'
    }
  });

  var description = ui.Label({
    value: (i + 1) + ': ' + CLASS_NAMES[i],
    style: {margin: '0 0 4px 6px'}
  });

  var row = ui.Panel({
    widgets: [colorBox, description],
    layout: ui.Panel.Layout.Flow('horizontal')
  });

  legend.add(row);
}

Map.add(legend);

// ============================================
// 15. エクスポート
// ============================================

Export.image.toDrive({
  image: classified.toInt8(),
  description: 'Mombetsu_LandUse_2020',
  folder: 'LandUse_Mombetsu',
  region: STUDY_AREA,
  scale: 30,
  maxPixels: 1e13,
  crs: 'EPSG:32654'
});

print('');
print('========================================');
print('分類完了！');
print('右側の Tasks タブからエクスポートを実行してください');
print('========================================');
