# 教師データ作成ガイド

## 概要

このドキュメントでは、土地利用分類のための教師データの作成方法と配置場所について説明します。

## 配置場所

```
data/training_data/
├── 2020/                          # 2020年用教師データ
│   ├── training_points.geojson    # ★ ここにGeoJSONを配置
│   ├── training_points.shp        # または Shapefile
│   └── training_points_template.geojson  # テンプレート（参考）
└── 2000/                          # 2000年用補助データ（任意）
    └── pseudo_training.geojson
```

## クラス定義

| クラスID | クラス名 | 日本語名 | 説明 |
|---------|---------|---------|------|
| 1 | grassland | 牧草地 | 牧草・採草地 |
| 2 | corn_field | 飼料用トウモロコシ畑 | デントコーン等 |
| 3 | other_cropland | 畑作農地 | 小麦・てんさい・馬鈴薯等 |
| 4 | forest | 森林 | 針葉樹・広葉樹・混交林 |
| 5 | urban | 市街地 | 建物・道路・人工構造物 |
| 6 | water | 水域 | 河川・湖沼・海 |
| 7 | bare_land | 裸地・その他 | 露岩・砂地・その他 |

## 教師データのフォーマット

### GeoJSON形式（推奨）

```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "properties": {
        "class": 1,
        "class_name": "grassland",
        "description": "牧草地サンプル",
        "confidence": "high"
      },
      "geometry": {
        "type": "Point",
        "coordinates": [143.5, 44.2]
      }
    }
  ]
}
```

### 必須プロパティ

- `class`: クラスID（1-7の整数）

### 推奨プロパティ

- `class_name`: クラス名（確認用）
- `description`: サンプルの説明
- `confidence`: 確信度（high/medium/low）
- `date_collected`: 収集日

### ジオメトリタイプ

- **Point**（推奨）: 個別のサンプルポイント
- **Polygon**: 均質な領域（内部からランダムサンプリング）

## 教師データの作成方法

### 方法1: QGIS を使用

1. QGISを起動し、衛星画像を読み込む
2. 新規レイヤー作成（GeoJSON/Shapefile）
3. 属性テーブルに `class` フィールド（整数）を追加
4. 各土地利用クラスのサンプルポイントをデジタイズ
5. GeoJSON形式でエクスポート

### 方法2: Google Earth Engine Code Editor を使用

```javascript
// GEE Code Editorでポイントを作成
var grassland = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([143.5, 44.2]), {class: 1}),
  ee.Feature(ee.Geometry.Point([143.51, 44.21]), {class: 1}),
]);

var corn = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([143.7, 44.1]), {class: 2}),
]);

// 全クラスを統合
var training = grassland.merge(corn);

// GeoJSONでエクスポート
Export.table.toDrive({
  collection: training,
  description: 'training_points_2020',
  fileFormat: 'GeoJSON'
});
```

### 方法3: Google Earth Pro を使用

1. Google Earth Proでポイントを作成
2. KML形式で保存
3. QGISまたはオンラインツールでGeoJSONに変換

## サンプル数のガイドライン

高精度な分類を達成するための推奨サンプル数:

| クラス | 最小サンプル数 | 推奨サンプル数 |
|-------|--------------|--------------|
| 牧草地 | 50 | 100-200 |
| トウモロコシ畑 | 50 | 100-200 |
| 畑作農地 | 50 | 100-200 |
| 森林 | 30 | 50-100 |
| 市街地 | 20 | 30-50 |
| 水域 | 20 | 30-50 |
| 裸地 | 20 | 30-50 |

### サンプリングの注意点

1. **空間的分散**: サンプルは対象地域全体に分散させる
2. **純粋なピクセル**: 複数クラスが混在しない均質な場所を選ぶ
3. **代表性**: 各クラスの様々なバリエーションを含める
4. **時期**: 2020年5-9月の画像で確認可能な場所を選ぶ

## 牧草地 vs トウモロコシ畑の識別ポイント

この分類で最も重要な識別:

### 牧草地の特徴
- 5-9月を通じて比較的安定したNDVI（0.5-0.7程度）
- 刈取り後にNDVIが一時的に低下することがある
- 区画が比較的不規則

### トウモロコシ畑の特徴
- 5月は低いNDVI（播種直後）
- 7-8月にNDVIが急上昇（0.8以上）
- 9月に収穫でNDVI低下
- 区画が整形で規則的

### 識別のコツ
- **8月のNDVI**が最も識別に有効
- トウモロコシは8月に牧草地より明らかに高いNDVI
- 時系列プロファイルの形状が異なる

## GEE Asset へのアップロード

大規模な教師データはGEE Assetにアップロード推奨:

1. [GEE Code Editor](https://code.earthengine.google.com/) にアクセス
2. Assets タブ → NEW → Table upload
3. GeoJSON/Shapefile をアップロード
4. Asset ID をスクリプトで指定

```python
# Pythonでの読み込み
training_data = ee.FeatureCollection('users/YOUR_USERNAME/training_2020')
```

## 精度向上のためのTips

1. **現地調査データ**: 可能であれば現地確認データを使用
2. **航空写真との照合**: Google Earth等の高解像度画像で確認
3. **既存データの活用**: 農地台帳、森林簿等の公的データ
4. **クロスバリデーション**: データを分割して精度検証
5. **複数回の検証**: 誤分類が多いクラスのサンプルを追加

## トラブルシューティング

### Q: GeoJSONが読み込めない
A: 以下を確認
- 座標系がWGS84 (EPSG:4326) か
- JSONの構文が正しいか
- `class` プロパティが整数か

### Q: サンプル数が足りないと言われる
A: 各クラス最低20サンプル以上を確保

### Q: 精度が低い
A: 以下を試す
- サンプルの品質を確認（混在ピクセルを除外）
- サンプル数を増やす
- 空間的な分散を改善
