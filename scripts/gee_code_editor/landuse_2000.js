/*
========================================================
北海道紋別市・湧別町 土地利用分類 2000年版
Google Earth Engine Code Editor (JavaScript)

Landsat 5 TM を使用（強化版雲マスキング）

修正点:
- 雲マスキングを大幅に強化
- 複数年のデータを使用（1999-2001年）
- 品質フィルタリングを追加
- 雲のないピクセルのみを使用
========================================================
*/

// ============================================
// 1. 設定パラメータ
// ============================================

var STUDY_AREA = ee.Geometry.Rectangle([143.0, 43.8, 144.5, 44.6]);

// 複数年のデータを使用（2000年中心±1年）
var TARGET_YEAR = 2000;
var USE_MULTI_YEAR = true;  // true: 1999-2001年を使用

var CLOUD_COVER_MAX = 30;  // シーン全体の雲量（マスキングで除去するので緩め）
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
// 2. 強化版雲マスキング関数
// ============================================

function maskL5CloudsStrict(image) {
  /*
  Landsat 5 Collection 2 QA_PIXEL ビット構成:
  Bit 0: Fill
  Bit 1: Dilated Cloud
  Bit 2: Cirrus (L8/L9のみ、L5では未使用)
  Bit 3: Cloud
  Bit 4: Cloud Shadow
  Bit 5: Snow
  Bit 6: Clear
  Bit 7: Water
  Bits 8-9: Cloud Confidence (0=None, 1=Low, 2=Medium, 3=High)
  Bits 10-11: Cloud Shadow Confidence
  Bits 12-13: Snow/Ice Confidence
  */

  var qa = image.select('QA_PIXEL');

  // 基本的な雲・雲影マスク
  var cloudBit = 1 << 3;
  var cloudShadowBit = 1 << 4;
  var dilatedCloudBit = 1 << 1;
  var snowBit = 1 << 5;

  // Cloud Confidence が Medium以上 (2 or 3) を除外
  // Bits 8-9: 00=None, 01=Low, 10=Medium, 11=High
  var cloudConfidenceBits = (3 << 8);  // bits 8-9
  var cloudConfidence = qa.bitwiseAnd(cloudConfidenceBits).rightShift(8);
  var lowCloudConfidence = cloudConfidence.lt(2);  // 0 or 1 のみ許可

  // Cloud Shadow Confidence
  var shadowConfidenceBits = (3 << 10);
  var shadowConfidence = qa.bitwiseAnd(shadowConfidenceBits).rightShift(10);
  var lowShadowConfidence = shadowConfidence.lt(2);

  // 複合マスク
  var mask = qa.bitwiseAnd(cloudBit).eq(0)
    .and(qa.bitwiseAnd(cloudShadowBit).eq(0))
    .and(qa.bitwiseAnd(dilatedCloudBit).eq(0))
    .and(qa.bitwiseAnd(snowBit).eq(0))
    .and(lowCloudConfidence)
    .and(lowShadowConfidence);

  return image.updateMask(mask);
}

// 追加の雲検出（スペクトル特性ベース）
function maskCloudsBySpectral(image) {
  // スケール適用後の反射率で判定
  var blue = image.select('SR_B1');
  var nir = image.select('SR_B4');
  var swir1 = image.select('SR_B5');

  // 雲は青バンドで高反射率
  var notBrightBlue = blue.lt(0.25);

  // 雲は NIR/SWIR 比が特定範囲
  // 雲: NIR ≈ SWIR (比率 ≈ 1)
  // 植生: NIR >> SWIR (比率 > 1.5)
  var nirSwirRatio = nir.divide(swir1.add(0.001));
  var notCloud = nirSwirRatio.gt(0.8);  // 雲っぽいものを除外

  // 異常に高い反射率を除外
  var notTooReflective = blue.lt(0.4).and(nir.lt(0.6));

  var spectralMask = notBrightBlue.and(notCloud).and(notTooReflective);

  return image.updateMask(spectralMask);
}

// スケールファクター適用
function applyScaleFactorsL5(image) {
  var opticalBands = image.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'])
    .multiply(0.0000275).add(-0.2);

  // 負の値をマスク（無効データ）
  var validMask = opticalBands.reduce(ee.Reducer.min()).gt(0);

  return image.addBands(opticalBands, null, true).updateMask(validMask);
}

// ============================================
// 3. データ取得
// ============================================

// 日付範囲の設定
var years = USE_MULTI_YEAR ? [1999, 2000, 2001] : [2000];
var collections = [];

years.forEach(function(year) {
  var startDate = year + '-05-01';
  var endDate = year + '-09-30';

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterBounds(STUDY_AREA)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_COVER_MAX));

  collections.push(col);
});

// コレクションを統合
var collection = ee.ImageCollection(collections[0]);
for (var i = 1; i < collections.length; i++) {
  collection = collection.merge(collections[i]);
}

// 雲マスキングとスケール適用
collection = collection
  .map(maskL5CloudsStrict)
  .map(applyScaleFactorsL5)
  .map(maskCloudsBySpectral);

print('取得画像数（雲マスキング前）:', collection.size());

// 有効ピクセル数でフィルタリング
var collectionFiltered = collection.map(function(image) {
  var validPixels = image.select('SR_B4').mask().reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: STUDY_AREA,
    scale: 1000,
    maxPixels: 1e9
  }).get('SR_B4');

  return image.set('validPixels', validPixels);
}).filter(ee.Filter.gt('validPixels', 1000));  // 最低限の有効ピクセル

print('有効画像数（フィルタリング後）:', collectionFiltered.size());

// ============================================
// 4. 植生指数の計算
// ============================================

function addIndicesL5(image) {
  var ndvi = image.normalizedDifference(['SR_B4', 'SR_B3']).rename('NDVI');

  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('SR_B4'),
      'RED': image.select('SR_B3'),
      'BLUE': image.select('SR_B1')
    }).rename('EVI');

  var ndwi = image.normalizedDifference(['SR_B2', 'SR_B4']).rename('NDWI');
  var lswi = image.normalizedDifference(['SR_B4', 'SR_B5']).rename('LSWI');
  var ndbi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDBI');

  return image.addBands([ndvi, evi, ndwi, lswi, ndbi]);
}

// ============================================
// 5. 品質重み付きコンポジット作成
// ============================================

// 各ピクセルの観測回数をカウント
var collectionWithIndices = collectionFiltered.map(addIndicesL5);
var obsCount = collectionWithIndices.select('NDVI').count().rename('obs_count');

// 中央値コンポジット（雲除去後）
var medianComposite = collectionWithIndices.median();

// 観測回数が少ない領域をマスク（信頼性低い）
var MIN_OBS = 3;  // 最低3回の観測が必要
var reliableMask = obsCount.gte(MIN_OBS);

print('観測回数マップを確認（Layers → obs_count）');

// ============================================
// 6. 月別コンポジット作成（品質フィルタ付き）
// ============================================

function createMonthlyCompositeL5(month) {
  // 複数年から該当月のデータを取得
  var monthlyImages = [];

  years.forEach(function(year) {
    var start = ee.Date.fromYMD(year, month, 1);
    var end = start.advance(1, 'month');
    var monthlyCol = collectionFiltered.filterDate(start, end);
    monthlyImages.push(monthlyCol);
  });

  // 全年の該当月データを統合
  var allMonthly = ee.ImageCollection(monthlyImages[0]);
  for (var i = 1; i < monthlyImages.length; i++) {
    allMonthly = allMonthly.merge(monthlyImages[i]);
  }

  var composite = allMonthly.median();
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
// 7. 時系列統計量
// ============================================

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
// 8. フェノロジー特徴量
// ============================================

// 複数年から月別データを取得
function getMonthlyNDVI(month) {
  var allMonthly = ee.ImageCollection([]);

  years.forEach(function(year) {
    var start = year + '-' + (month < 10 ? '0' : '') + month + '-01';
    var endMonth = month + 1;
    var end = year + '-' + (endMonth < 10 ? '0' : '') + endMonth + '-01';
    if (month === 9) {
      end = year + '-09-30';
    }

    var col = collectionFiltered.filterDate(start, end).map(addIndicesL5);
    allMonthly = allMonthly.merge(col);
  });

  return allMonthly.select('NDVI').median();
}

var ndviJuly = getMonthlyNDVI(7).rename('NDVI_july');
var ndviAug = getMonthlyNDVI(8).rename('NDVI_aug');
var ndviSept = getMonthlyNDVI(9).rename('NDVI_sept');

var ndviDiffAugJuly = ndviAug.subtract(ndviJuly).rename('NDVI_diff_aug_july');
var ndviDiffSeptAug = ndviSept.subtract(ndviAug).rename('NDVI_diff_sept_aug');
var summerMaxNdvi = ndviJuly.max(ndviAug).rename('NDVI_summer_max');

var phenology = ndviJuly.addBands(ndviAug).addBands(ndviSept)
  .addBands(ndviDiffAugJuly).addBands(ndviDiffSeptAug).addBands(summerMaxNdvi);

// ============================================
// 9. 地形データ
// ============================================

var dem = ee.Image('USGS/SRTMGL1_003');
var elevation = dem.select('elevation').rename('elevation');
var slope = ee.Terrain.slope(dem).rename('slope');
var aspect = ee.Terrain.aspect(dem).rename('aspect');
var aspectSin = aspect.multiply(Math.PI / 180).sin().rename('aspect_sin');
var aspectCos = aspect.multiply(Math.PI / 180).cos().rename('aspect_cos');

var terrain = elevation.addBands(slope).addBands(aspectSin).addBands(aspectCos);

// ============================================
// 10. 特徴量スタック（信頼性マスク適用）
// ============================================

var featureStack = monthlyStack
  .addBands(temporalStats)
  .addBands(phenology)
  .addBands(terrain)
  .updateMask(reliableMask)  // 観測回数が少ない領域をマスク
  .clip(STUDY_AREA);

print('特徴量バンド:', featureStack.bandNames());

// ============================================
// 11. 擬似教師データの作成
// ============================================

print('');
print('========================================');
print('擬似教師データ作成中...');
print('（雲除去済みデータから作成）');
print('========================================');

// 閾値ベースの擬似分類
var waterMask = ndwiMean.gt(0);
var forestMask = ndviMax.gt(0.7).and(elevation.gt(100));
var urbanMask = ndviMax.lt(0.3).and(ndwiMean.lt(0));
var croplandMask = ndviMax.gt(0.5).and(elevation.lt(100));
var cornMask = croplandMask.and(ndviAug.gt(0.75));
var grassMask = croplandMask.and(ndviAug.gte(0.5)).and(ndviAug.lte(0.75));
var otherCropMask = croplandMask.and(ndviAug.lt(0.5));
var bareMask = ndviMax.lt(0.2).and(ndwiMean.lt(0));

var pseudoClass = ee.Image(0)
  .where(bareMask, 7)
  .where(otherCropMask, 3)
  .where(grassMask, 1)
  .where(cornMask, 2)
  .where(urbanMask, 5)
  .where(forestMask, 4)
  .where(waterMask, 6)
  .rename('class')
  .updateMask(reliableMask)  // 信頼性マスク適用
  .clip(STUDY_AREA);

// ============================================
// 12. Random Forest 分類
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

var trainingWithFeatures = featureStack.sampleRegions({
  collection: trainingSamples,
  properties: ['class'],
  scale: 30,
  tileScale: 8
});

var classifier = ee.Classifier.smileRandomForest(RF_TREES)
  .train({
    features: trainingWithFeatures,
    classProperty: 'class',
    inputProperties: featureStack.bandNames()
  });

var classified = featureStack.classify(classifier);

// ============================================
// 13. 地図表示
// ============================================

Map.setCenter(143.7, 44.2, 10);

// True Color（雲除去済み）
var composite = collectionWithIndices.median();
var visTC = {bands: ['SR_B3', 'SR_B2', 'SR_B1'], min: 0, max: 0.3};
Map.addLayer(composite.clip(STUDY_AREA), visTC, 'True Color 2000（雲除去済み）');

// 観測回数マップ
var visObs = {min: 1, max: 20, palette: ['red', 'yellow', 'green']};
Map.addLayer(obsCount.clip(STUDY_AREA), visObs, '観測回数', false);

// NDVI（8月）
var visNDVI = {min: 0.2, max: 0.9, palette: ['red', 'yellow', 'green', 'darkgreen']};
Map.addLayer(ndviAug.updateMask(reliableMask).clip(STUDY_AREA), visNDVI, '8月 NDVI');

// 分類結果
var visClass = {min: 1, max: 7, palette: CLASS_PALETTE};
Map.addLayer(classified.clip(STUDY_AREA), visClass, '土地利用分類 2000');

// 信頼性マスク（観測回数が少ない領域）
var unreliableArea = reliableMask.not().selfMask();
Map.addLayer(unreliableArea.clip(STUDY_AREA), {palette: ['white']}, 'データ不足領域', false);

// 対象地域
Map.addLayer(ee.Image().paint(STUDY_AREA, 0, 2), {palette: 'blue'}, '対象地域');

// ============================================
// 14. 凡例
// ============================================

var legend = ui.Panel({
  style: {position: 'bottom-left', padding: '8px 15px'}
});

legend.add(ui.Label({
  value: '土地利用分類 2000年',
  style: {fontWeight: 'bold', fontSize: '16px', margin: '0 0 4px 0'}
}));

legend.add(ui.Label({
  value: '※ 雲除去済み・複数年データ使用',
  style: {fontSize: '10px', color: 'green', margin: '0 0 8px 0'}
}));

for (var i = 0; i < CLASS_NAMES.length; i++) {
  var row = ui.Panel({
    widgets: [
      ui.Label({style: {backgroundColor: '#' + CLASS_PALETTE[i], padding: '8px', margin: '0 0 4px 0'}}),
      ui.Label({value: (i + 1) + ': ' + CLASS_NAMES[i], style: {margin: '0 0 4px 6px'}})
    ],
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
  description: 'Mombetsu_LandUse_2000_CloudFree',
  folder: 'LandUse_Mombetsu',
  region: STUDY_AREA,
  scale: 30,
  maxPixels: 1e13,
  crs: 'EPSG:32654'
});

print('');
print('========================================');
print('✅ 雲除去強化版');
print('========================================');
print('改善点:');
print('1. QA_PIXEL の全ビットを使用した厳密な雲マスキング');
print('2. スペクトル特性による追加の雲検出');
print('3. 1999-2001年の複数年データを使用');
print('4. 観測回数が少ない領域をマスク');
print('');
print('白い領域 = データ不足（雲が多すぎて有効データなし）');
print('========================================');
