#!/usr/bin/env python3
"""
北海道紋別市・湧別町 土地利用分類 - 2000年版
Landsat 5 TM Collection 2 SR を使用した分類

【概要】
- 2000年の土地利用分類
- 2020年モデルの転用（方法A）または擬似教師データ（方法B）を使用
- Landsat 5 TM のバンド構成に対応

【精度の限界】
- 2000年は直接的な教師データがないため、精度に限界がある
- 2020年モデルを転用する場合、時代による土地利用パターンの変化の影響を受ける
- 擬似教師データを使用する場合、参照データの精度に依存する
"""

import ee
import numpy as np
import pandas as pd
import geopandas as gpd
from datetime import datetime
import json
import os

# ============================================================
# 設定パラメータ
# ============================================================

# 対象地域の設定（紋別市・湧別町の範囲）
STUDY_AREA = {
    'min_lon': 143.0,
    'max_lon': 144.5,
    'min_lat': 43.8,
    'max_lat': 44.6
}

# 分類カテゴリ
LAND_USE_CLASSES = {
    1: 'grassland',           # 牧草地
    2: 'corn_field',          # 飼料用トウモロコシ畑
    3: 'other_cropland',      # 畑作農地（その他）
    4: 'forest',              # 森林
    5: 'urban',               # 市街地
    6: 'water',               # 水域
    7: 'bare_land'            # 裸地・その他
}

# 分析対象年と期間
TARGET_YEAR = 2000
START_DATE = f'{TARGET_YEAR}-05-01'
END_DATE = f'{TARGET_YEAR}-09-30'

# 雲量閾値
CLOUD_COVER_MAX = 20

# Random Forest パラメータ
RF_PARAMS = {
    'numberOfTrees': 100,
    'minLeafPopulation': 5,
    'bagFraction': 0.7,
    'seed': 42
}


# ============================================================
# Landsat 5 TM と Landsat 8 OLI のバンド対応
# ============================================================
"""
Landsat 5 TM:           Landsat 8 OLI:
SR_B1 (Blue)     <->    SR_B2 (Blue)
SR_B2 (Green)    <->    SR_B3 (Green)
SR_B3 (Red)      <->    SR_B4 (Red)
SR_B4 (NIR)      <->    SR_B5 (NIR)
SR_B5 (SWIR1)    <->    SR_B6 (SWIR1)
SR_B7 (SWIR2)    <->    SR_B7 (SWIR2)

※ Landsat 8 には Coastal/Aerosol (B1) と Cirrus (B9) が追加
"""

# 共通バンドのマッピング
BAND_MAPPING = {
    'blue': {'l5': 'SR_B1', 'l8': 'SR_B2'},
    'green': {'l5': 'SR_B2', 'l8': 'SR_B3'},
    'red': {'l5': 'SR_B3', 'l8': 'SR_B4'},
    'nir': {'l5': 'SR_B4', 'l8': 'SR_B5'},
    'swir1': {'l5': 'SR_B5', 'l8': 'SR_B6'},
    'swir2': {'l5': 'SR_B7', 'l8': 'SR_B7'}
}


# ============================================================
# Google Earth Engine 初期化
# ============================================================

def initialize_gee():
    """Google Earth Engineを初期化"""
    try:
        ee.Initialize()
        print("GEE初期化成功")
    except Exception as e:
        print(f"GEE認証が必要です: {e}")
        ee.Authenticate()
        ee.Initialize()


# ============================================================
# 対象地域の定義
# ============================================================

def get_study_area():
    """紋別市・湧別町の対象地域を定義"""
    geometry = ee.Geometry.Rectangle([
        STUDY_AREA['min_lon'], STUDY_AREA['min_lat'],
        STUDY_AREA['max_lon'], STUDY_AREA['max_lat']
    ])
    return geometry


# ============================================================
# Landsat 5 データ取得・前処理
# ============================================================

def mask_landsat5_clouds(image):
    """
    Landsat 5 Collection 2 SR の雲マスキング
    QA_PIXELバンドを使用
    """
    qa = image.select('QA_PIXEL')

    # Bit 3: Cloud, Bit 4: Cloud Shadow
    cloud_bit_mask = 1 << 3
    cloud_shadow_bit_mask = 1 << 4

    mask = qa.bitwiseAnd(cloud_bit_mask).eq(0) \
        .And(qa.bitwiseAnd(cloud_shadow_bit_mask).eq(0))

    return image.updateMask(mask)


def apply_scale_factors_l5(image):
    """
    Landsat 5 Collection 2 SR のスケールファクター適用
    """
    optical_bands = image.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']) \
        .multiply(0.0000275).add(-0.2)

    return image.addBands(optical_bands, overwrite=True)


def get_landsat5_collection(geometry, start_date, end_date, cloud_cover_max):
    """
    Landsat 5 Collection 2 SR データを取得
    """
    collection = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2') \
        .filterBounds(geometry) \
        .filterDate(start_date, end_date) \
        .filter(ee.Filter.lt('CLOUD_COVER', cloud_cover_max)) \
        .map(mask_landsat5_clouds) \
        .map(apply_scale_factors_l5)

    print(f"取得画像数: {collection.size().getInfo()}")
    return collection


# ============================================================
# 植生指数計算（Landsat 5用）
# ============================================================

def add_vegetation_indices_l5(image):
    """
    Landsat 5 TM 用の植生指数計算

    Landsat 5 バンド対応:
    - Blue: SR_B1
    - Green: SR_B2
    - Red: SR_B3
    - NIR: SR_B4
    - SWIR1: SR_B5
    - SWIR2: SR_B7
    """
    # NDVI: (NIR - Red) / (NIR + Red)
    ndvi = image.normalizedDifference(['SR_B4', 'SR_B3']).rename('NDVI')

    # EVI: 2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1)
    evi = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
        {
            'NIR': image.select('SR_B4'),
            'RED': image.select('SR_B3'),
            'BLUE': image.select('SR_B1')
        }
    ).rename('EVI')

    # NDWI: (Green - NIR) / (Green + NIR)
    ndwi = image.normalizedDifference(['SR_B2', 'SR_B4']).rename('NDWI')

    # LSWI: (NIR - SWIR1) / (NIR + SWIR1)
    lswi = image.normalizedDifference(['SR_B4', 'SR_B5']).rename('LSWI')

    # NDBI: (SWIR1 - NIR) / (SWIR1 + NIR)
    ndbi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDBI')

    return image.addBands([ndvi, evi, ndwi, lswi, ndbi])


# ============================================================
# 月別コンポジット作成（Landsat 5用）
# ============================================================

def create_monthly_composites_l5(collection, geometry, year):
    """
    5-9月の月別中央値コンポジットを作成（Landsat 5用）
    """
    months = [5, 6, 7, 8, 9]
    composites = []

    for month in months:
        start = f'{year}-{month:02d}-01'
        if month == 9:
            end = f'{year}-{month:02d}-30'
        else:
            end = f'{year}-{month+1:02d}-01'

        monthly_col = collection.filterDate(start, end)

        # 月別中央値コンポジット
        composite = monthly_col.median().set('month', month)

        # 植生指数を追加
        composite = add_vegetation_indices_l5(composite)

        # バンド名にプレフィックスを追加
        # Landsat 5: B1(Blue), B2(Green), B3(Red), B4(NIR), B5(SWIR1), B7(SWIR2)
        bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7',
                 'NDVI', 'EVI', 'NDWI', 'LSWI', 'NDBI']
        renamed_bands = [f'M{month}_{b}' for b in bands]

        composite = composite.select(bands).rename(renamed_bands)
        composites.append(composite)

    # 全月をスタック
    stacked = composites[0]
    for comp in composites[1:]:
        stacked = stacked.addBands(comp)

    return stacked


# ============================================================
# 時系列統計量の計算（Landsat 5用）
# ============================================================

def calculate_temporal_statistics_l5(collection, geometry):
    """
    植生指数の時系列統計量を計算（Landsat 5用）
    """
    collection_with_indices = collection.map(add_vegetation_indices_l5)

    indices = ['NDVI', 'EVI', 'NDWI', 'LSWI']
    stats_image = ee.Image()

    for idx in indices:
        idx_collection = collection_with_indices.select(idx)

        max_val = idx_collection.max().rename(f'{idx}_max')
        min_val = idx_collection.min().rename(f'{idx}_min')
        mean_val = idx_collection.mean().rename(f'{idx}_mean')
        std_val = idx_collection.reduce(ee.Reducer.stdDev()).rename(f'{idx}_std')
        cv = std_val.divide(mean_val.abs().add(0.001)).rename(f'{idx}_cv')

        stats_image = stats_image.addBands([max_val, min_val, mean_val, std_val, cv])

    return stats_image


def calculate_phenology_features_l5(collection, geometry, year):
    """
    フェノロジー特徴量を計算（Landsat 5用）
    """
    july_col = collection.filterDate(f'{year}-07-01', f'{year}-07-31').map(add_vegetation_indices_l5)
    aug_col = collection.filterDate(f'{year}-08-01', f'{year}-08-31').map(add_vegetation_indices_l5)
    sept_col = collection.filterDate(f'{year}-09-01', f'{year}-09-30').map(add_vegetation_indices_l5)

    july_ndvi = july_col.select('NDVI').median().rename('NDVI_july')
    aug_ndvi = aug_col.select('NDVI').median().rename('NDVI_aug')
    sept_ndvi = sept_col.select('NDVI').median().rename('NDVI_sept')

    ndvi_diff_aug_july = aug_ndvi.subtract(july_ndvi).rename('NDVI_diff_aug_july')
    ndvi_diff_sept_aug = sept_ndvi.subtract(aug_ndvi).rename('NDVI_diff_sept_aug')
    ndvi_growth_rate = sept_ndvi.subtract(july_ndvi) \
        .divide(july_ndvi.abs().add(0.001)).rename('NDVI_growth_rate')
    summer_max_ndvi = july_col.merge(aug_col).select('NDVI').max().rename('NDVI_summer_max')

    phenology = ee.Image.cat([
        july_ndvi, aug_ndvi, sept_ndvi,
        ndvi_diff_aug_july, ndvi_diff_sept_aug,
        ndvi_growth_rate, summer_max_ndvi
    ])

    return phenology


# ============================================================
# 地形データ
# ============================================================

def add_terrain_features(geometry):
    """
    SRTM DEMから地形特徴量を計算
    """
    dem = ee.Image('USGS/SRTMGL1_003')

    elevation = dem.select('elevation').rename('elevation')
    slope = ee.Terrain.slope(dem).rename('slope')
    aspect = ee.Terrain.aspect(dem).rename('aspect')
    aspect_sin = aspect.multiply(np.pi / 180).sin().rename('aspect_sin')
    aspect_cos = aspect.multiply(np.pi / 180).cos().rename('aspect_cos')

    terrain = ee.Image.cat([elevation, slope, aspect_sin, aspect_cos])

    return terrain


# ============================================================
# 特徴量スタック作成（Landsat 5用）
# ============================================================

def create_feature_stack_l5(geometry, year=2000):
    """
    Landsat 5用の特徴量スタックを作成
    """
    print("Landsat 5 データを取得中...")
    collection = get_landsat5_collection(
        geometry,
        f'{year}-05-01',
        f'{year}-09-30',
        CLOUD_COVER_MAX
    )

    print("月別コンポジットを作成中...")
    monthly_composites = create_monthly_composites_l5(collection, geometry, year)

    print("時系列統計量を計算中...")
    temporal_stats = calculate_temporal_statistics_l5(collection, geometry)

    print("フェノロジー特徴量を計算中...")
    phenology = calculate_phenology_features_l5(collection, geometry, year)

    print("地形データを追加中...")
    terrain = add_terrain_features(geometry)

    feature_stack = ee.Image.cat([
        monthly_composites,
        temporal_stats,
        phenology,
        terrain
    ]).clip(geometry)

    print("特徴量スタック作成完了")
    return feature_stack


# ============================================================
# 方法A: 2020年モデルの転用
# ============================================================

def create_common_band_feature_stack_l5(geometry, year=2000):
    """
    2020年モデルと互換性のある共通バンドのみの特徴量スタックを作成

    ※ Landsat 5とLandsat 8で共通して計算できる特徴量のみを使用
    - 植生指数（NDVI, EVI, NDWI, LSWI）は共通
    - 地形データは共通
    - 月別コンポジットは植生指数のみ使用
    """
    print("Landsat 5 データを取得中（共通バンド版）...")
    collection = get_landsat5_collection(
        geometry,
        f'{year}-05-01',
        f'{year}-09-30',
        CLOUD_COVER_MAX
    )

    # 植生指数のみの月別コンポジット
    months = [5, 6, 7, 8, 9]
    composites = []

    for month in months:
        start = f'{year}-{month:02d}-01'
        if month == 9:
            end = f'{year}-{month:02d}-30'
        else:
            end = f'{year}-{month+1:02d}-01'

        monthly_col = collection.filterDate(start, end)
        composite = monthly_col.median().set('month', month)
        composite = add_vegetation_indices_l5(composite)

        # 植生指数のみ（共通特徴量）
        bands = ['NDVI', 'EVI', 'NDWI', 'LSWI', 'NDBI']
        renamed_bands = [f'M{month}_{b}' for b in bands]
        composite = composite.select(bands).rename(renamed_bands)
        composites.append(composite)

    stacked = composites[0]
    for comp in composites[1:]:
        stacked = stacked.addBands(comp)

    # 時系列統計
    temporal_stats = calculate_temporal_statistics_l5(collection, geometry)

    # フェノロジー
    phenology = calculate_phenology_features_l5(collection, geometry, year)

    # 地形
    terrain = add_terrain_features(geometry)

    feature_stack = ee.Image.cat([
        stacked,
        temporal_stats,
        phenology,
        terrain
    ]).clip(geometry)

    return feature_stack


def transfer_model_from_2020(classifier_2020, feature_stack_2000, common_band_names):
    """
    2020年のモデルを2000年に適用

    注意:
    - 2020年と2000年で共通のバンド名である必要がある
    - モデル転用の精度は限定的
    """
    classified = feature_stack_2000.select(common_band_names).classify(classifier_2020)
    return classified


# ============================================================
# 方法B: 擬似教師データの作成
# ============================================================

def create_pseudo_training_from_existing_data(geometry, year=2000):
    """
    既存のデータセットから擬似教師データを作成

    利用可能なデータソース:
    1. 国土数値情報 土地利用細分メッシュデータ
    2. 農林業センサス（集落別データ）
    3. 森林簿データ
    4. Copernicus Global Land Cover
    5. ESA WorldCover (2020年のみ、参考)

    ※ これらのデータは別途ダウンロード・アップロードが必要
    """
    print("擬似教師データを作成中...")

    # 方法1: Copernicus Global Land Cover (2000年に近い2015年版を使用)
    # 注: 完全に同年ではないため、精度に注意
    try:
        cgls = ee.ImageCollection('COPERNICUS/Landcover/100m/Proba-V-C3/Global') \
            .filterDate('2015-01-01', '2015-12-31') \
            .first()

        # CGLSの土地被覆クラスを本プロジェクトのクラスにマッピング
        # CGLS classes: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V-C3_Global
        remapped = cgls.select('discrete_classification').remap(
            [0, 20, 30, 40, 50, 60, 70, 80, 90, 100, 111, 112, 113, 114, 115, 116, 121, 122, 123, 124, 125, 126, 200],
            [7, 3, 3, 3, 5, 7, 7, 6, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 6],  # プロジェクトクラスにマッピング
            7  # デフォルト値
        ).rename('class')

        print("Copernicus Global Land Cover からマッピング完了")
        return remapped

    except Exception as e:
        print(f"CGLS読み込みエラー: {e}")
        return None


def create_pseudo_training_from_spectral_thresholds(feature_stack, geometry):
    """
    スペクトル閾値を使用した擬似教師データ作成

    この方法は教師なし分類に近いアプローチ:
    1. NDVI閾値で植生/非植生を分離
    2. NDWI閾値で水域を抽出
    3. NDBI閾値で市街地を抽出
    4. 8月NDVIの高低で牧草地/トウモロコシを分離
    """
    print("スペクトル閾値による擬似クラス作成中...")

    # 必要なバンドを取得
    ndvi_max = feature_stack.select('NDVI_max')
    ndwi_mean = feature_stack.select('NDWI_mean')
    ndbi_mean = feature_stack.select('NDBI_mean') if 'NDBI_mean' in feature_stack.bandNames().getInfo() else None
    ndvi_aug = feature_stack.select('NDVI_aug')
    elevation = feature_stack.select('elevation')

    # 閾値ベースの分類
    # 水域: NDWI > 0
    water = ndwi_mean.gt(0).multiply(6)

    # 森林: NDVI_max > 0.7 かつ elevation > 100m（平野部を除外）
    forest = ndvi_max.gt(0.7).And(elevation.gt(100)).multiply(4)

    # 市街地: NDBI > 0 かつ NDVI < 0.3
    if ndbi_mean is not None:
        urban = ndbi_mean.gt(0).And(ndvi_max.lt(0.3)).multiply(5)
    else:
        urban = ndvi_max.lt(0.2).multiply(5)

    # 農地（高NDVI平野部）
    cropland_mask = ndvi_max.gt(0.5).And(elevation.lt(100))

    # トウモロコシ: 8月NDVIが特に高い（> 0.8）
    corn = cropland_mask.And(ndvi_aug.gt(0.8)).multiply(2)

    # 牧草地: 8月NDVIが中程度（0.5-0.8）
    grassland = cropland_mask.And(ndvi_aug.gte(0.5)).And(ndvi_aug.lte(0.8)).multiply(1)

    # その他農地
    other_crop = cropland_mask.And(ndvi_aug.lt(0.5)).multiply(3)

    # 裸地/その他
    bare = ndvi_max.lt(0.2).And(ndwi_mean.lt(0)).multiply(7)

    # 統合（優先度順）
    pseudo_class = water \
        .where(forest.gt(0), forest) \
        .where(urban.gt(0), urban) \
        .where(corn.gt(0), corn) \
        .where(grassland.gt(0), grassland) \
        .where(other_crop.gt(0), other_crop) \
        .where(bare.gt(0), bare)

    return pseudo_class.rename('class')


def sample_pseudo_training_points(pseudo_class_image, geometry, samples_per_class=500):
    """
    擬似分類画像からサンプルポイントを抽出
    """
    # 各クラスからストラタファイドサンプリング
    samples = pseudo_class_image.stratifiedSample(
        numPoints=samples_per_class,
        classBand='class',
        region=geometry,
        scale=30,
        seed=42,
        geometries=True
    )

    return samples


# ============================================================
# Random Forest 分類
# ============================================================

def train_random_forest(feature_stack, training_data, band_names):
    """Random Forest分類器を訓練"""
    training_samples = feature_stack.sampleRegions(
        collection=training_data,
        properties=['class'],
        scale=30,
        tileScale=8
    )

    print(f"訓練サンプル数: {training_samples.size().getInfo()}")

    classifier = ee.Classifier.smileRandomForest(**RF_PARAMS)
    trained_classifier = classifier.train(
        features=training_samples,
        classProperty='class',
        inputProperties=band_names
    )

    return trained_classifier


def classify_image(feature_stack, classifier):
    """画像を分類"""
    classified = feature_stack.classify(classifier)
    return classified


# ============================================================
# 精度評価
# ============================================================

def evaluate_with_reference_data(classified, reference_data, geometry):
    """
    参照データとの比較による精度評価

    ※ 2000年の実測データがない場合は、
    別時点の土地利用データとの比較になるため参考値として扱う
    """
    if reference_data is None:
        print("警告: 参照データがないため精度評価をスキップします")
        return None

    # サンプリング
    samples = classified.addBands(reference_data).sample(
        region=geometry,
        scale=30,
        numPixels=5000,
        seed=42
    )

    # 混同行列
    confusion = samples.errorMatrix('class', 'classification')

    return {
        'confusion_matrix': confusion.getInfo(),
        'overall_accuracy': confusion.accuracy().getInfo(),
        'kappa': confusion.kappa().getInfo()
    }


# ============================================================
# 結果のエクスポート
# ============================================================

def export_classification_to_drive(classified_image, geometry, description, folder='LandUse_Mombetsu'):
    """分類結果をGoogle Driveにエクスポート"""
    task = ee.batch.Export.image.toDrive(
        image=classified_image.toInt8(),
        description=description,
        folder=folder,
        region=geometry,
        scale=30,
        maxPixels=1e13,
        crs='EPSG:32654'
    )
    task.start()
    print(f"エクスポートタスク開始: {description}")
    return task


# ============================================================
# メイン処理
# ============================================================

def main_method_a(classifier_2020=None):
    """
    方法A: 2020年モデルの転用

    2020年で訓練したモデルを2000年に適用
    ※ 事前に landuse_classification_2020.py でモデルを訓練する必要がある
    """
    print("="*60)
    print("方法A: 2020年モデルの転用による2000年分類")
    print("="*60)
    print("\n【精度の限界】")
    print("- 2020年の土地利用パターンを基準としているため、")
    print("  2000年当時の特有のパターンを捉えられない可能性がある")
    print("- 20年間の土地利用変化により、同じスペクトル特性でも")
    print("  異なる土地利用の可能性がある")

    initialize_gee()
    geometry = get_study_area()

    # 共通バンドの特徴量スタック作成
    print("\n特徴量スタック作成（共通バンド版）...")
    feature_stack = create_common_band_feature_stack_l5(geometry, TARGET_YEAR)

    band_names = feature_stack.bandNames().getInfo()
    print(f"特徴量数: {len(band_names)}")

    if classifier_2020 is None:
        print("\n警告: 2020年のモデルが提供されていません。")
        print("先に landuse_classification_2020.py を実行してモデルを取得してください。")
        return None

    # モデル転用による分類
    print("\n2020年モデルを適用中...")
    classified = feature_stack.classify(classifier_2020)

    # エクスポート
    export_classification_to_drive(
        classified, geometry,
        description='Mombetsu_LandUse_2000_MethodA'
    )

    return classified


def main_method_b():
    """
    方法B: 擬似教師データによる分類

    既存の土地被覆データやスペクトル閾値から擬似教師データを作成
    """
    print("="*60)
    print("方法B: 擬似教師データによる2000年分類")
    print("="*60)
    print("\n【精度の限界】")
    print("- 擬似教師データは実測値ではないため、本質的な精度限界がある")
    print("- 特に農地細分類（牧草地/トウモロコシ/畑作）の精度は")
    print("  スペクトル閾値の設定に大きく依存する")
    print("- 結果は「参考値」として扱い、変化検出等の")
    print("  相対的な比較に使用することを推奨")

    initialize_gee()
    geometry = get_study_area()

    # 特徴量スタック作成
    print("\n特徴量スタック作成...")
    feature_stack = create_feature_stack_l5(geometry, TARGET_YEAR)

    band_names = feature_stack.bandNames().getInfo()
    print(f"特徴量数: {len(band_names)}")

    # 擬似教師データ作成（スペクトル閾値版）
    print("\n擬似教師データ作成...")
    pseudo_class = create_pseudo_training_from_spectral_thresholds(feature_stack, geometry)

    # サンプルポイント抽出
    print("\nサンプルポイント抽出中...")
    training_samples = sample_pseudo_training_points(pseudo_class, geometry, samples_per_class=500)

    # Random Forest訓練・分類
    print("\nRandom Forest分類...")
    classifier = train_random_forest(feature_stack, training_samples, band_names)
    classified = classify_image(feature_stack, classifier)

    # エクスポート
    export_classification_to_drive(
        classified, geometry,
        description='Mombetsu_LandUse_2000_MethodB'
    )

    print("\n処理完了!")
    return classified, classifier


def main():
    """
    2000年土地利用分類のメイン処理

    方法A（モデル転用）と方法B（擬似教師データ）の両方を実行可能
    """
    print("="*60)
    print("北海道紋別市・湧別町 土地利用分類 2000")
    print("="*60)

    # 方法Bを実行（単独で実行可能）
    classified, classifier = main_method_b()

    return classified, classifier


if __name__ == '__main__':
    classified, classifier = main()
