# 北海道紋別市・湧別町 土地利用分類

Landsat衛星画像を使用した土地利用分類プロジェクト

## 概要

北海道紋別市と湧別町を対象に、2000年と2020年の土地利用を識別するプロジェクトです。
特に農地の細分類（牧草地・飼料用トウモロコシ畑・畑作農地）の高精度識別を目指します。

### 対象クラス

| ID | クラス | 説明 |
|----|-------|------|
| 1 | 牧草地 | 牧草・採草地 |
| 2 | 飼料用トウモロコシ畑 | デントコーン等 |
| 3 | 畑作農地 | 小麦・てんさい・馬鈴薯等 |
| 4 | 森林 | 針葉樹・広葉樹・混交林 |
| 5 | 市街地 | 建物・道路・人工構造物 |
| 6 | 水域 | 河川・湖沼・海 |
| 7 | 裸地・その他 | 露岩・砂地・その他 |

## プロジェクト構成

```
remotesensing_mombetsu/
├── README.md                      # このファイル
├── requirements.txt               # Python依存パッケージ
├── scripts/
│   ├── landuse_classification_2020.py  # 2020年分類スクリプト
│   └── landuse_classification_2000.py  # 2000年分類スクリプト
└── data/
    ├── training_data/
    │   ├── TRAINING_DATA_GUIDE.md      # 教師データ作成ガイド
    │   ├── 2020/
    │   │   ├── training_points.geojson     # ★教師データをここに配置★
    │   │   └── training_points_template.geojson  # テンプレート
    │   └── 2000/
    │       └── (補助データ)
    └── output/
        ├── 2020/                   # 2020年出力結果
        └── 2000/                   # 2000年出力結果
```

## クイックスタート

### 1. 環境準備

```bash
# 仮想環境作成（推奨）
python -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate   # Windows

# 依存パッケージインストール
pip install -r requirements.txt
```

### 2. Google Earth Engine 認証

```bash
earthengine authenticate
```

### 3. 教師データの準備

**重要**: 分類を実行する前に教師データを準備してください。

1. `data/training_data/TRAINING_DATA_GUIDE.md` を参照
2. GeoJSON形式で教師データを作成
3. `data/training_data/2020/training_points.geojson` として保存

### 4. 分類の実行

```bash
# 2020年の分類
python scripts/landuse_classification_2020.py

# 2000年の分類
python scripts/landuse_classification_2000.py
```

## 使用データ

### 衛星データ

| 年 | センサー | データセット |
|----|---------|------------|
| 2000 | Landsat 5 TM | LANDSAT/LT05/C02/T1_L2 |
| 2020 | Landsat 8 OLI | LANDSAT/LC08/C02/T1_L2 |

- Collection 2 Surface Reflectance (SR) を使用
- 対象期間: 5-9月（生育期）
- 雲量閾値: 20%未満

### 特徴量

1. **月別コンポジット** (5-9月)
   - 分光バンド (Blue, Green, Red, NIR, SWIR1, SWIR2)
   - 植生指数 (NDVI, EVI, NDWI, LSWI, NDBI)

2. **時系列統計量**
   - 最大値、最小値、平均値、標準偏差、変動係数

3. **フェノロジー特徴量**
   - 7月/8月/9月のNDVI
   - 月間NDVI差分
   - 成長率
   - 夏季最大NDVI

4. **地形データ** (SRTM 30m DEM)
   - 標高、傾斜、方位

## 分類手法

### 2020年（教師あり）

- Random Forest分類器
- 訓練データ: 現地調査データ等
- 精度評価: 交差検証

### 2000年（教師データなし）

2つの方法を提供:

**方法A: モデル転用**
- 2020年で訓練したモデルを適用
- 共通バンドのみ使用
- 精度限界: 20年間の変化の影響

**方法B: 擬似教師データ**
- スペクトル閾値による自動分類
- 既存土地被覆データからのサンプリング
- 精度限界: 参照データの精度に依存

## 牧草地 vs トウモロコシ畑の識別

この分類で最も重要な識別ポイント:

```
        牧草地              トウモロコシ畑
5月   |====NDVI中=====|   |__NDVI低__(播種直後)
6月   |====NDVI中=====|   |/ーNDVI上昇ー\
7月   |====NDVI中=====|   |/===NDVI高===\
8月   |====NDVI中=====|   |====NDVI最高====| ← 最も識別しやすい
9月   |====NDVI中=====|   |__NDVI低下__(収穫)
```

- **8月のNDVI**が最も識別に有効
- トウモロコシは8月にNDVI > 0.8
- 牧草地は0.5-0.7程度で安定

## 出力

- 分類結果: Google Driveにエクスポート (GeoTIFF)
- 精度評価: コンソール出力（混同行列、Overall Accuracy、Kappa）
- 特徴量重要度: 上位特徴量をリスト表示

## 対象地域

- 紋別市（北海道オホーツク総合振興局）
- 湧別町（北海道オホーツク総合振興局）
- 概略範囲: 143.0°E - 144.5°E, 43.8°N - 44.6°N

## 注意事項

### 精度の限界

- 2000年は直接的な教師データがないため、精度に限界がある
- 牧草地とトウモロコシの識別は7-8月の画像品質に依存
- 雲の多い年は精度が低下する可能性

### データ要件

- Google Earth Engine アカウントが必要
- 大規模な処理にはGEEの計算リソース制限に注意

## ライセンス

MIT License

## 参考文献

- Landsat Collection 2: https://www.usgs.gov/landsat-missions/landsat-collection-2
- Google Earth Engine: https://earthengine.google.com/
