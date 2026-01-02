/*
========================================================
北海道紋別市・湧別町 土地利用分類 2000年版
Google Earth Engine Code Editor (JavaScript)

Landsat 5 TM を使用
擬似教師データ（スペクトル閾値）による分類

使い方:
1. https://code.earthengine.google.com/ を開く
2. このコードをコピー&ペースト
3. 「Run」ボタンをクリック

⚠️ 精度の限界:
- 2000年は直接的な教師データがないため精度に限界がある
- 結果は「参考値」として扱い、変化検出等に使用推奨
========================================================
*/

// ============================================
// 1. 設定パラメータ
// ============================================

var STUDY_AREA = ee.Geometry.Rectangle([143.0, 43.8, 144.5, 44.6]);

var TARGET_YEAR = 2000;
var START_DATE = TARGET_YEAR + '-05-01';
var END_DATE = TARGET_YEAR + '-09-30';

var CLOUD_COVER_MAX = 20;
var RF_TREES = 100;

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
  '90EE90',  // 1: 牧草地
  'FFD700',  // 2: トウモロコシ
  'FFA500',  // 3: 畑作
  '006400',  // 4: 森林
  'FF0000',  // 5: 市街地
  '0000FF',  // 6: 水域
  '808080'   // 7: 裸地
];

// ============================================
// 2. Landsat 5 TM データ取得・前処理
// ============================================

/*
Landsat 5 TM vs Landsat 8 OLI バンド対応:
  Blue:  L5=SR_B1, L8=SR_B2
  Green: L5=SR_B2, L8=SR_B3
  Red:   L5=SR_B3, L8=SR_B4
  NIR:   L5=SR_B4, L8=SR_B5
  SWIR1: L5=SR_B5, L8=SR_B6
  SWIR2: L5=SR_B7, L8=SR_B7
*/

function maskL5Clouds(image) {
  var qa = image.select('QA_PIXEL');
  var cloudBitMask = 1 << 3;
  var cloudShadowBitMask = 1 << 4;

  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0));

  return image.updateMask(mask);
}

function applyScaleFactorsL5(image) {
  var opticalBands = image.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'])
    .multiply(0.0000275).add(-0.2);

  return image.addBands(opticalBands, null, true);
}

var collection = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(STUDY_AREA)
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_COVER_MAX))
  .map(maskL5Clouds)
  .map(applyScaleFactorsL5);

print('取得画像数:', collection.size());

// ============================================
// 3. 植生指数の計算（Landsat 5用）
// ============================================

function addIndicesL5(image) {
  // NDVI (NIR=B4, Red=B3)
  var ndvi = image.normalizedDifference(['SR_B4', 'SR_B3']).rename('NDVI');

  // EVI
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('SR_B4'),
      'RED': image.select('SR_B3'),
      'BLUE': image.select('SR_B1')
    }).rename('EVI');

  // NDWI (Green=B2, NIR=B4)
  var ndwi = image.normalizedDifference(['SR_B2', 'SR_B4']).rename('NDWI');

  // LSWI (NIR=B4, SWIR1=B5)
  var lswi = image.normalizedDifference(['SR_B4', 'SR_B5']).rename('LSWI');

  // NDBI (SWIR1=B5, NIR=B4)
  var ndbi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDBI');

  return image.addBands([ndvi, evi, ndwi, lswi, ndbi]);
}

// ============================================
// 4. 月別コンポジット作成
// ============================================

function createMonthlyCompositeL5(month) {
  var start = ee.Date.fromYMD(TARGET_YEAR, month, 1);
  var end = start.advance(1, 'month');

  var monthlyCol = collection.filterDate(start, end);
  var composite = monthlyCol.median();
  composite = addIndicesL5(composite);

  var bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7',
               'NDVI', 'EVI', 'NDWI', 'LSWI', 'NDBI'];
  var prefix = 'M' + month + '_';
  var renamedBands = bands.map(function(b) { return prefix + b; });

  return composite.select(bands).rename(renamedBands);
}

var m5 = createMonthlyCompositeL5(5);
var m6 = createMonthlyCompositeL5(6);
var m7 = createMonthlyCompositeL5(7);
var m8 = createMonthlyCompositeL5(8);
var m9 = createMonthlyCompositeL5(9);

var monthlyStack = m5.addBands(m6).addBands(m7).addBands(m8).addBands(m9);

// ============================================
// 5. 時系列統計量
// ============================================

var collectionWithIndices = collection.map(addIndicesL5);

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
// 6. フェノロジー特徴量
// ============================================

var julyCol = collection.filterDate(TARGET_YEAR + '-07-01', TARGET_YEAR + '-07-31').map(addIndicesL5);
var augCol = collection.filterDate(TARGET_YEAR + '-08-01', TARGET_YEAR + '-08-31').map(addIndicesL5);
var septCol = collection.filterDate(TARGET_YEAR + '-09-01', TARGET_YEAR + '-09-30').map(addIndicesL5);

var ndviJuly = julyCol.select('NDVI').median().rename('NDVI_july');
var ndviAug = augCol.select('NDVI').median().rename('NDVI_aug');
var ndviSept = septCol.select('NDVI').median().rename('NDVI_sept');

var ndviDiffAugJuly = ndviAug.subtract(ndviJuly).rename('NDVI_diff_aug_july');
var ndviDiffSeptAug = ndviSept.subtract(ndviAug).rename('NDVI_diff_sept_aug');
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
// 9. 擬似教師データの作成（スペクトル閾値）
// ============================================

print('');
print('========================================');
print('擬似教師データ作成中...');
print('※ スペクトル閾値に基づく自動分類');
print('========================================');

// 閾値ベースの擬似分類
// 水域: NDWI > 0
var waterMask = ndwiMean.gt(0);

// 森林: NDVI_max > 0.7 かつ elevation > 100m
var forestMask = ndviMax.gt(0.7).and(elevation.gt(100));

// 市街地: NDVI_max < 0.3
var urbanMask = ndviMax.lt(0.3).and(ndwiMean.lt(0));

// 農地（高NDVI平野部）
var croplandMask = ndviMax.gt(0.5).and(elevation.lt(100));

// トウモロコシ: 8月NDVIが高い（> 0.75）
var cornMask = croplandMask.and(ndviAug.gt(0.75));

// 牧草地: 8月NDVIが中程度（0.5-0.75）
var grassMask = croplandMask.and(ndviAug.gte(0.5)).and(ndviAug.lte(0.75));

// その他農地
var otherCropMask = croplandMask.and(ndviAug.lt(0.5));

// 裸地
var bareMask = ndviMax.lt(0.2).and(ndwiMean.lt(0));

// 擬似クラス画像を作成
var pseudoClass = ee.Image(0)
  .where(bareMask, 7)
  .where(otherCropMask, 3)
  .where(grassMask, 1)
  .where(cornMask, 2)
  .where(urbanMask, 5)
  .where(forestMask, 4)
  .where(waterMask, 6)
  .rename('class')
  .clip(STUDY_AREA);

// ============================================
// 10. 擬似教師データからサンプリング
// ============================================

var trainingSamples = pseudoClass.stratifiedSample({
  numPoints: 300,
  classBand: 'class',
  region: STUDY_AREA,
  scale: 30,
  seed: 42,
  geometries: true
});

print('擬似教師サンプル数:', trainingSamples.size());

// ============================================
// 11. Random Forest 分類
// ============================================

// 教師サンプルにバンド値を抽出
var trainingWithFeatures = featureStack.sampleRegions({
  collection: trainingSamples,
  properties: ['class'],
  scale: 30,
  tileScale: 8
});

// 分類器を訓練
var classifier = ee.Classifier.smileRandomForest(RF_TREES)
  .train({
    features: trainingWithFeatures,
    classProperty: 'class',
    inputProperties: featureStack.bandNames()
  });

// 分類実行
var classified = featureStack.classify(classifier);

// ============================================
// 12. 地図表示
// ============================================

Map.setCenter(143.7, 44.2, 10);

// True Color (Landsat 5: B3=Red, B2=Green, B1=Blue)
var composite = collection.median();
var visTC = {bands: ['SR_B3', 'SR_B2', 'SR_B1'], min: 0, max: 0.3};
Map.addLayer(composite.clip(STUDY_AREA), visTC, 'True Color 2000');

// NDVI（8月）
var visNDVI = {min: 0.2, max: 0.9, palette: ['red', 'yellow', 'green', 'darkgreen']};
Map.addLayer(ndviAug.clip(STUDY_AREA), visNDVI, '8月 NDVI');

// 擬似分類（閾値ベース）
var visClass = {min: 1, max: 7, palette: CLASS_PALETTE};
Map.addLayer(pseudoClass, visClass, '擬似分類（閾値ベース）', false);

// Random Forest 分類結果
Map.addLayer(classified.clip(STUDY_AREA), visClass, '土地利用分類 2000');

// 対象地域の境界
Map.addLayer(ee.Image().paint(STUDY_AREA, 0, 2), {palette: 'blue'}, '対象地域');

// ============================================
// 13. 凡例パネル
// ============================================

var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

var legendTitle = ui.Label({
  value: '土地利用分類 2000年',
  style: {fontWeight: 'bold', fontSize: '16px', margin: '0 0 4px 0'}
});
legend.add(legendTitle);

var note = ui.Label({
  value: '※ 擬似教師データ使用',
  style: {fontSize: '10px', color: 'red', margin: '0 0 8px 0'}
});
legend.add(note);

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
// 14. エクスポート
// ============================================

Export.image.toDrive({
  image: classified.toInt8(),
  description: 'Mombetsu_LandUse_2000',
  folder: 'LandUse_Mombetsu',
  region: STUDY_AREA,
  scale: 30,
  maxPixels: 1e13,
  crs: 'EPSG:32654'
});

// ============================================
// 15. 精度の限界についての注意
// ============================================

print('');
print('========================================');
print('⚠️ 精度の限界について');
print('========================================');
print('');
print('2000年の分類は擬似教師データ（スペクトル閾値）を');
print('使用しているため、以下の限界があります:');
print('');
print('1. 直接的な教師データがない');
print('2. 閾値は経験則に基づいており、最適化されていない');
print('3. 2000年当時の土地利用パターンと異なる可能性');
print('');
print('結果は「参考値」として扱い、');
print('変化検出等の相対的な比較に使用してください。');
print('========================================');
print('');
print('右側の Tasks タブからエクスポートを実行してください');
