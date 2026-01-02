#!/usr/bin/env python3
"""
北海道紋別市・湧別町 土地利用分類 - 2020年版
Landsat 8 Collection 2 SR を使用したRandom Forest分類

【概要】
- 農地（牧草地・飼料用トウモロコシ畑・畑作農地）・森林・市街地・水域等の識別
- 5-9月の時系列コンポジットを使用
- Random Forestで植生指数の時系列統計と地形データを特徴量として使用

【教師データの配置場所】
data/training_data/2020/ に以下のフォーマットでGeoJSONまたはShapefileを配置
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
# 実際の座標は調整が必要
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
TARGET_YEAR = 2020
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
    # 簡易的な矩形範囲（実際は市町村境界を使用推奨）
    geometry = ee.Geometry.Rectangle([
        STUDY_AREA['min_lon'], STUDY_AREA['min_lat'],
        STUDY_AREA['max_lon'], STUDY_AREA['max_lat']
    ])
    return geometry


def get_mombetsu_yubetsu_boundary():
    """
    紋別市・湧別町の正確な境界を取得
    注: GEE上の日本行政区界データセットを使用
    """
    # FAO GAUL または他の行政区界データを使用
    # 日本の市町村境界データがGEEにない場合は、
    # 独自にアップロードしたアセットを使用

    # 例: カスタムアセットを使用する場合
    # boundary = ee.FeatureCollection('users/YOUR_USERNAME/mombetsu_yubetsu')

    # 簡易版: 矩形範囲を使用
    return get_study_area()


# ============================================================
# Landsat 8 データ取得・前処理
# ============================================================

def mask_landsat8_clouds(image):
    """
    Landsat 8 Collection 2 SR の雲マスキング
    QA_PIXELバンドを使用
    """
    # QA_PIXELバンドのビット情報
    qa = image.select('QA_PIXEL')

    # Bit 3: Cloud, Bit 4: Cloud Shadow
    cloud_bit_mask = 1 << 3
    cloud_shadow_bit_mask = 1 << 4

    # 雲と雲影のないピクセルを抽出
    mask = qa.bitwiseAnd(cloud_bit_mask).eq(0) \
        .And(qa.bitwiseAnd(cloud_shadow_bit_mask).eq(0))

    return image.updateMask(mask)


def apply_scale_factors(image):
    """
    Landsat 8 Collection 2 SR のスケールファクター適用
    """
    # Surface Reflectance バンド
    optical_bands = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']) \
        .multiply(0.0000275).add(-0.2)

    # Surface Temperature バンド（オプション）
    thermal_band = image.select(['ST_B10']).multiply(0.00341802).add(149.0)

    return image.addBands(optical_bands, overwrite=True) \
        .addBands(thermal_band, overwrite=True)


def get_landsat8_collection(geometry, start_date, end_date, cloud_cover_max):
    """
    Landsat 8 Collection 2 SR データを取得
    """
    collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
        .filterBounds(geometry) \
        .filterDate(start_date, end_date) \
        .filter(ee.Filter.lt('CLOUD_COVER', cloud_cover_max)) \
        .map(mask_landsat8_clouds) \
        .map(apply_scale_factors)

    print(f"取得画像数: {collection.size().getInfo()}")
    return collection


# ============================================================
# 植生指数計算
# ============================================================

def add_vegetation_indices(image):
    """
    植生指数を計算して追加

    Landsat 8 バンド対応:
    - Blue: SR_B2
    - Green: SR_B3
    - Red: SR_B4
    - NIR: SR_B5
    - SWIR1: SR_B6
    - SWIR2: SR_B7
    """
    # NDVI: (NIR - Red) / (NIR + Red)
    ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI')

    # EVI: 2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1)
    evi = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
        {
            'NIR': image.select('SR_B5'),
            'RED': image.select('SR_B4'),
            'BLUE': image.select('SR_B2')
        }
    ).rename('EVI')

    # NDWI: (Green - NIR) / (Green + NIR) - 水域検出
    ndwi = image.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI')

    # LSWI: (NIR - SWIR1) / (NIR + SWIR1) - 植生水分
    lswi = image.normalizedDifference(['SR_B5', 'SR_B6']).rename('LSWI')

    # NDBI: (SWIR1 - NIR) / (SWIR1 + NIR) - 市街地検出
    ndbi = image.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI')

    return image.addBands([ndvi, evi, ndwi, lswi, ndbi])


# ============================================================
# 月別コンポジット作成
# ============================================================

def create_monthly_composites(collection, geometry, year):
    """
    5-9月の月別中央値コンポジットを作成
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
        composite = add_vegetation_indices(composite)

        # バンド名にプレフィックスを追加
        bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7',
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
# 時系列統計量の計算
# ============================================================

def calculate_temporal_statistics(collection, geometry):
    """
    植生指数の時系列統計量を計算
    - 最大値、最小値、平均値、標準偏差、変動係数
    - 7-9月のNDVI差分（トウモロコシ識別用）
    """
    # 植生指数を追加した画像コレクション
    collection_with_indices = collection.map(add_vegetation_indices)

    # 統計量計算
    indices = ['NDVI', 'EVI', 'NDWI', 'LSWI']
    stats_image = ee.Image()

    for idx in indices:
        idx_collection = collection_with_indices.select(idx)

        # 基本統計
        max_val = idx_collection.max().rename(f'{idx}_max')
        min_val = idx_collection.min().rename(f'{idx}_min')
        mean_val = idx_collection.mean().rename(f'{idx}_mean')
        std_val = idx_collection.reduce(ee.Reducer.stdDev()).rename(f'{idx}_std')

        # 変動係数 (CV = std / mean)
        cv = std_val.divide(mean_val.abs().add(0.001)).rename(f'{idx}_cv')

        stats_image = stats_image.addBands([max_val, min_val, mean_val, std_val, cv])

    return stats_image


def calculate_phenology_features(collection, geometry, year):
    """
    フェノロジー（生育段階）特徴量を計算
    牧草地 vs トウモロコシ識別に重要
    """
    # 7月、8月、9月のNDVIを取得
    july_col = collection.filterDate(f'{year}-07-01', f'{year}-07-31').map(add_vegetation_indices)
    aug_col = collection.filterDate(f'{year}-08-01', f'{year}-08-31').map(add_vegetation_indices)
    sept_col = collection.filterDate(f'{year}-09-01', f'{year}-09-30').map(add_vegetation_indices)

    july_ndvi = july_col.select('NDVI').median().rename('NDVI_july')
    aug_ndvi = aug_col.select('NDVI').median().rename('NDVI_aug')
    sept_ndvi = sept_col.select('NDVI').median().rename('NDVI_sept')

    # トウモロコシは7-8月にNDVIが急上昇（牧草地より高い）
    # 8月-7月のNDVI差分
    ndvi_diff_aug_july = aug_ndvi.subtract(july_ndvi).rename('NDVI_diff_aug_july')

    # 9月-8月のNDVI差分
    ndvi_diff_sept_aug = sept_ndvi.subtract(aug_ndvi).rename('NDVI_diff_sept_aug')

    # 7-9月のNDVI変化率
    ndvi_growth_rate = sept_ndvi.subtract(july_ndvi) \
        .divide(july_ndvi.abs().add(0.001)).rename('NDVI_growth_rate')

    # 夏季（7-8月）のNDVI最大値
    summer_max_ndvi = july_col.merge(aug_col).select('NDVI').max().rename('NDVI_summer_max')

    phenology = ee.Image.cat([
        july_ndvi, aug_ndvi, sept_ndvi,
        ndvi_diff_aug_july, ndvi_diff_sept_aug,
        ndvi_growth_rate, summer_max_ndvi
    ])

    return phenology


# ============================================================
# 地形データの追加
# ============================================================

def add_terrain_features(geometry):
    """
    SRTM DEMから地形特徴量を計算
    """
    # SRTM 30m DEM
    dem = ee.Image('USGS/SRTMGL1_003')

    # 標高
    elevation = dem.select('elevation').rename('elevation')

    # 傾斜
    slope = ee.Terrain.slope(dem).rename('slope')

    # 方位
    aspect = ee.Terrain.aspect(dem).rename('aspect')

    # 方位を sin/cos に変換（連続値として扱う）
    aspect_sin = aspect.multiply(np.pi / 180).sin().rename('aspect_sin')
    aspect_cos = aspect.multiply(np.pi / 180).cos().rename('aspect_cos')

    terrain = ee.Image.cat([elevation, slope, aspect_sin, aspect_cos])

    return terrain


# ============================================================
# 特徴量スタックの作成
# ============================================================

def create_feature_stack(geometry, year=2020):
    """
    分類用の特徴量スタックを作成
    """
    print("Landsat 8 データを取得中...")
    collection = get_landsat8_collection(
        geometry,
        f'{year}-05-01',
        f'{year}-09-30',
        CLOUD_COVER_MAX
    )

    print("月別コンポジットを作成中...")
    monthly_composites = create_monthly_composites(collection, geometry, year)

    print("時系列統計量を計算中...")
    temporal_stats = calculate_temporal_statistics(collection, geometry)

    print("フェノロジー特徴量を計算中...")
    phenology = calculate_phenology_features(collection, geometry, year)

    print("地形データを追加中...")
    terrain = add_terrain_features(geometry)

    # 全特徴量をスタック
    feature_stack = ee.Image.cat([
        monthly_composites,
        temporal_stats,
        phenology,
        terrain
    ]).clip(geometry)

    print("特徴量スタック作成完了")
    return feature_stack


# ============================================================
# 教師データの読み込み
# ============================================================

def load_training_data_from_geojson(geojson_path):
    """
    GeoJSONファイルから教師データを読み込み

    【教師データのフォーマット】
    各フィーチャーには 'class' または 'landuse' プロパティが必要:
    - 1: grassland (牧草地)
    - 2: corn_field (飼料用トウモロコシ畑)
    - 3: other_cropland (畑作農地)
    - 4: forest (森林)
    - 5: urban (市街地)
    - 6: water (水域)
    - 7: bare_land (裸地・その他)
    """
    gdf = gpd.read_file(geojson_path)

    # GeoJSONをGEE FeatureCollectionに変換
    features = []
    for idx, row in gdf.iterrows():
        geom = row.geometry.__geo_interface__
        props = {'class': int(row.get('class', row.get('landuse', 0)))}
        feature = ee.Feature(ee.Geometry(geom), props)
        features.append(feature)

    return ee.FeatureCollection(features)


def load_training_data_from_asset(asset_id):
    """
    GEE Asset から教師データを読み込み

    Asset にアップロードする場合:
    1. Google Cloud にShapefile/GeoJSONをアップロード
    2. GEE Code Editor で Asset として登録
    3. asset_id を指定して読み込み
    """
    return ee.FeatureCollection(asset_id)


def create_sample_training_data(geometry):
    """
    サンプル教師データを作成（テスト用）

    実際の分析では、以下の方法で教師データを作成:
    1. Google Earth Engine Code Editor で手動でポイントを作成
    2. QGIS/ArcGIS でポリゴンを作成しGeoJSONでエクスポート
    3. 現地調査データをGeoJSON/Shapefileに変換
    """
    # サンプルポイント（実際のデータに置き換える必要あり）
    sample_points = ee.FeatureCollection([
        # 牧草地サンプル
        ee.Feature(ee.Geometry.Point([143.5, 44.2]), {'class': 1}),
        ee.Feature(ee.Geometry.Point([143.6, 44.3]), {'class': 1}),

        # トウモロコシ畑サンプル
        ee.Feature(ee.Geometry.Point([143.7, 44.1]), {'class': 2}),
        ee.Feature(ee.Geometry.Point([143.8, 44.2]), {'class': 2}),

        # 森林サンプル
        ee.Feature(ee.Geometry.Point([143.4, 44.4]), {'class': 4}),
        ee.Feature(ee.Geometry.Point([143.5, 44.5]), {'class': 4}),

        # 水域サンプル
        ee.Feature(ee.Geometry.Point([143.9, 44.0]), {'class': 6}),
    ])

    print("警告: サンプル教師データを使用中。実際の分析では正確な教師データを使用してください。")
    return sample_points


# ============================================================
# Random Forest 分類
# ============================================================

def train_random_forest(feature_stack, training_data, band_names):
    """
    Random Forest分類器を訓練
    """
    # 教師データにバンド値を抽出
    training_samples = feature_stack.sampleRegions(
        collection=training_data,
        properties=['class'],
        scale=30,
        tileScale=8
    )

    print(f"訓練サンプル数: {training_samples.size().getInfo()}")

    # Random Forest 分類器
    classifier = ee.Classifier.smileRandomForest(**RF_PARAMS)

    # 訓練
    trained_classifier = classifier.train(
        features=training_samples,
        classProperty='class',
        inputProperties=band_names
    )

    return trained_classifier


def classify_image(feature_stack, classifier):
    """
    画像を分類
    """
    classified = feature_stack.classify(classifier)
    return classified


# ============================================================
# 精度評価
# ============================================================

def evaluate_accuracy(feature_stack, training_data, classifier, band_names, test_split=0.3):
    """
    分類精度を評価
    - 全体精度 (Overall Accuracy)
    - カッパ係数 (Kappa Coefficient)
    - 混同行列 (Confusion Matrix)
    - クラス別精度 (Producer's/User's Accuracy)
    """
    # サンプルデータをランダムに分割
    training_data_with_random = training_data.randomColumn('random')
    train_set = training_data_with_random.filter(ee.Filter.lt('random', 1 - test_split))
    test_set = training_data_with_random.filter(ee.Filter.gte('random', 1 - test_split))

    # 訓練データでサンプリング
    train_samples = feature_stack.sampleRegions(
        collection=train_set,
        properties=['class'],
        scale=30,
        tileScale=8
    )

    # テストデータでサンプリング
    test_samples = feature_stack.sampleRegions(
        collection=test_set,
        properties=['class'],
        scale=30,
        tileScale=8
    )

    # 分類器を訓練
    classifier = ee.Classifier.smileRandomForest(**RF_PARAMS).train(
        features=train_samples,
        classProperty='class',
        inputProperties=band_names
    )

    # テストデータで分類
    validated = test_samples.classify(classifier)

    # 混同行列を計算
    confusion_matrix = validated.errorMatrix('class', 'classification')

    # 精度指標
    accuracy = confusion_matrix.accuracy()
    kappa = confusion_matrix.kappa()
    producers_accuracy = confusion_matrix.producersAccuracy()
    consumers_accuracy = confusion_matrix.consumersAccuracy()

    return {
        'confusion_matrix': confusion_matrix.getInfo(),
        'overall_accuracy': accuracy.getInfo(),
        'kappa': kappa.getInfo(),
        'producers_accuracy': producers_accuracy.getInfo(),
        'consumers_accuracy': consumers_accuracy.getInfo()
    }


def print_accuracy_report(accuracy_results):
    """
    精度評価レポートを出力
    """
    print("\n" + "="*60)
    print("精度評価レポート")
    print("="*60)

    print(f"\n全体精度 (Overall Accuracy): {accuracy_results['overall_accuracy']:.4f}")
    print(f"カッパ係数 (Kappa): {accuracy_results['kappa']:.4f}")

    print("\n混同行列:")
    cm = accuracy_results['confusion_matrix']
    print(np.array(cm))

    print("\nProducer's Accuracy (クラス別):")
    for i, acc in enumerate(accuracy_results['producers_accuracy']):
        if i < len(LAND_USE_CLASSES):
            class_name = list(LAND_USE_CLASSES.values())[i]
            print(f"  {class_name}: {acc[0]:.4f}")

    print("\nUser's Accuracy (クラス別):")
    for i, acc in enumerate(accuracy_results['consumers_accuracy']):
        if i < len(LAND_USE_CLASSES):
            class_name = list(LAND_USE_CLASSES.values())[i]
            print(f"  {class_name}: {acc[0]:.4f}")


# ============================================================
# 特徴量重要度の分析
# ============================================================

def get_feature_importance(classifier, band_names):
    """
    Random Forestの特徴量重要度を取得
    """
    importance = classifier.explain()
    importance_dict = importance.getInfo()

    if 'importance' in importance_dict:
        imp_values = importance_dict['importance']
        importance_df = pd.DataFrame({
            'feature': band_names,
            'importance': [imp_values.get(b, 0) for b in band_names]
        }).sort_values('importance', ascending=False)

        return importance_df
    else:
        print("特徴量重要度を取得できませんでした")
        return None


# ============================================================
# 結果のエクスポート
# ============================================================

def export_classification_to_drive(classified_image, geometry, description, folder='LandUse_Mombetsu'):
    """
    分類結果をGoogle Driveにエクスポート
    """
    task = ee.batch.Export.image.toDrive(
        image=classified_image.toInt8(),
        description=description,
        folder=folder,
        region=geometry,
        scale=30,
        maxPixels=1e13,
        crs='EPSG:32654'  # UTM Zone 54N
    )
    task.start()
    print(f"エクスポートタスク開始: {description}")
    print(f"タスクID: {task.id}")
    return task


def export_classification_to_asset(classified_image, geometry, asset_id):
    """
    分類結果をGEE Assetにエクスポート
    """
    task = ee.batch.Export.image.toAsset(
        image=classified_image.toInt8(),
        description=asset_id.split('/')[-1],
        assetId=asset_id,
        region=geometry,
        scale=30,
        maxPixels=1e13,
        crs='EPSG:32654'
    )
    task.start()
    print(f"Assetエクスポートタスク開始: {asset_id}")
    return task


# ============================================================
# メイン処理
# ============================================================

def main():
    """
    2020年土地利用分類のメイン処理
    """
    print("="*60)
    print("北海道紋別市・湧別町 土地利用分類 2020")
    print("="*60)

    # GEE初期化
    initialize_gee()

    # 対象地域
    geometry = get_study_area()
    print(f"\n対象地域: 紋別市・湧別町")
    print(f"範囲: {STUDY_AREA}")

    # 特徴量スタック作成
    print("\n" + "-"*40)
    print("特徴量作成")
    print("-"*40)
    feature_stack = create_feature_stack(geometry, TARGET_YEAR)

    # バンド名を取得
    band_names = feature_stack.bandNames().getInfo()
    print(f"特徴量数: {len(band_names)}")

    # =========================================================
    # 教師データの読み込み
    # =========================================================
    print("\n" + "-"*40)
    print("教師データ読み込み")
    print("-"*40)

    # 教師データのパス設定
    # 以下のいずれかの方法で教師データを読み込む:

    # 方法1: ローカルのGeoJSONファイル
    training_data_path = 'data/training_data/2020/training_points.geojson'

    if os.path.exists(training_data_path):
        print(f"教師データを読み込み中: {training_data_path}")
        training_data = load_training_data_from_geojson(training_data_path)
    else:
        # 方法2: GEE Asset
        # asset_id = 'users/YOUR_USERNAME/mombetsu_training_2020'
        # training_data = load_training_data_from_asset(asset_id)

        # 方法3: サンプルデータ（テスト用）
        print("警告: 教師データが見つかりません。サンプルデータを使用します。")
        print(f"教師データを配置してください: {training_data_path}")
        training_data = create_sample_training_data(geometry)

    # =========================================================
    # Random Forest 分類
    # =========================================================
    print("\n" + "-"*40)
    print("Random Forest 分類")
    print("-"*40)

    # 分類器を訓練
    classifier = train_random_forest(feature_stack, training_data, band_names)

    # 分類実行
    classified = classify_image(feature_stack, classifier)
    print("分類完了")

    # =========================================================
    # 精度評価
    # =========================================================
    print("\n" + "-"*40)
    print("精度評価")
    print("-"*40)

    accuracy_results = evaluate_accuracy(
        feature_stack, training_data, classifier, band_names
    )
    print_accuracy_report(accuracy_results)

    # =========================================================
    # 特徴量重要度
    # =========================================================
    print("\n" + "-"*40)
    print("特徴量重要度")
    print("-"*40)

    importance_df = get_feature_importance(classifier, band_names)
    if importance_df is not None:
        print("\n上位20の重要な特徴量:")
        print(importance_df.head(20).to_string(index=False))

    # =========================================================
    # 結果のエクスポート
    # =========================================================
    print("\n" + "-"*40)
    print("結果エクスポート")
    print("-"*40)

    # Google Driveにエクスポート
    export_task = export_classification_to_drive(
        classified,
        geometry,
        description='Mombetsu_LandUse_2020'
    )

    print("\n処理完了!")
    print("エクスポートはGEEのタスクマネージャーで確認してください。")

    return classified, classifier, accuracy_results


if __name__ == '__main__':
    classified, classifier, accuracy = main()
