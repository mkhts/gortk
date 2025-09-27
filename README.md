# gortk

(English is below.)

# 概要
- RTK-GNSS の勉強のために書いたコードです。RTK-GNSS のアルゴリズムを Go 言語で実装しています。拙い実装で恐縮ですが、無償で公開致します。
- 単独測位、DGPS までの処理は、書籍「[GPSのための実用プログラミング](https://www.tdupress.jp/book/b349356.html)」を参考にしました。非常にわかりやすい本で GNSS のアルゴリズムについて理解したい人には大変オススメの書籍です。著者である坂井氏に大変感謝致します。
- RTK の処理については、言わずと知れた「RTKLIB」のコードと、雑誌トランジスタ技術に付属していた「RTKコア」のコード、rtkeplorer氏の RTKLIB Demo5 のコードを大いに参考にさせて頂きました。高須氏、久保氏、rtkexplorer氏に感謝致します。LAMBDA法および対流圏遅延補正に関するコードについてはそのまま移植して利用させて頂きました。
- 行列演算については、[gonum/matライブラリ](https://pkg.go.dev/gonum.org/v1/gonum/mat)を利用しています。

# 処理について
- RINEX ファイルを読み込んで処理するオフライン型のプログラムです。
関数的にはエポック単位で処理するようになっていますが、リアルタイムにデータをストリームから読み込んで処理する実装にはなっていません。
- 処理時間を意識したコードになっていないことと、整数値アンビギュイティを解く際に Fix するまで衛星ペアを変えて何度かリトライする処理になっているので、実行時間は、RTKLIB に比べて倍程度に長くなっています。
- RTKLIB でいう AR モードは、instantaneous（エポック独立処理） と continueous（カルマンフィルタ処理）モードを実装しています。fix & hold モードに対応する処理はありません。
- カルマンフィルターの状態ベクトルとして、未知である３次元位置と、各衛星についての基準局と移動局とでの整数値アンビギュイティの一重差を未知数としています。当初は、教科書通りに整数値アンビギュイティの二重差を未知数としていたのですが、時間経過に応じて最大仰角衛星が変わるたびに二重差の組み合わせが変化するため、カルマンフィルタの誤差共分散行列（P）を都度リセットする必要があり、そのたびに測位解が飛んでしまうので、最大仰角衛星に依らない一重差方式に変更しました。RTKLIB と同じです。
- Float 解の算出後、整数値アンビギュイティを解く際に、１回目のトライで fix しなかった場合は、仰角が低い衛星を１つずつ外しながら、最大 10 回までリトライしています。それでも Fix しなかった場合は、Glonass を外して、再度 Float 解の算出からやり直しています。それでも Fix しない場合に Float になります。
- 基準局と移動局との時間差である age については、オフラインということもあって、マイナス age　すなわち未来の基準局のデータも許容しています。（最大±15秒）
- エポック独立である instantaneous モードの際も、Float 解の算出にはカルマンフィルタの計算を用いています。カルマンフィルタの観測更新を繰り返すことで（3回）、観測方程式を非線形最小二乗法の繰り返し計算で解く場合と実質的に同じ計算結果になるようにしています。
- カルマンフィルタの時間更新については、基本的には静止モデルで、ダイナミズムは入っていません。
- モードによらず、状態ベクトルの３次元位置については、DGPS の測位解を、観測更新前の初期推定値としてセットしています。前エポックの値をそのまま据え置くスタティックモードには対応していません。
- 入力として読み込む RINEX ファイルのバージョンは、3.02 もしくは 3.04 を想定しています。

# ファイル構成

主たるファイルは以下の４つです。

| ファイル | 内容 |
| ------ | ------ |
| cmd/gortk/main.go | メイン処理 |
| calcspp.go | 単独測位計算 |
| calcfloat.go | Float解の計算 |
| solveamb.go | Fix解の計算 |

その他のファイルは以下の通りです。

| ファイル | 内容 |
| ------ | ------ |
| satpos.go | 衛星の３次元位置計算 |
| rinex.go | RINEX ファイルの読み込み |
| obs.go | 観測データをしまう構造体の定義と関連処理 |
| nav.go | 航法データをしまう構造体の定義と関連処理 |
| pos.go | 位置情報をしまう構造体の定義と関連処理 |
| gtime.go | GPS時刻をしまう構造体の定義と関連処理 |
| const.go | 各種定数値の定義 |
| solvels.go | 最小二乗の行列計算（Gonumライブラリ利用）|
| lambda.go | LAMBDA 法|
| tropo.go | 対流圏遅延補正|
| misc.go | ミニ関数など |

# ビルド
```
cd cmd/gortk
go build
```

# 使い方

### 単独測位
```
gortk -p 0 rover.obs base.nav
```

### DGPS
```
gortk -p 1 -l "35.73101206 139.7396917 80.33" rover.obs base.obs base.nav
```
* -l オプションに与える基準局位置は、ダブルクオーテーションで囲む必要があります。

### RTK
```
gortk -p 2 -l "35.73101206 139.7396917 80.33" rover.obs base.obs base.nav
```

### オプション表示
```
gortk -h
```
- RTKLIB と共通するものはなるべく同じオプション名になるようにしています。
- ただし、-l、-ts、-te オプションについては、その後に続く引数をダブルクオーテーションで囲む必要があります。

# Overview
- This is code I wrote to study RTK-GNSS. It implements RTK-GNSS algorithms in Go. Although the implementation is rather modest, I am releasing it publicly and free of charge.  
- For single-point positioning and DGPS, I referred to the book ["Practical Programming for GPS"](https://www.tdupress.jp/book/b349356.html). It is an excellent and easy-to-understand resource, highly recommended for anyone who wants to gain a deeper understanding of GNSS algorithms. I am deeply grateful to the author, Mr. Sakai.  
- For RTK processing, I drew heavily on the well-known "RTKLIB" source code, the "RTK Core" code that came with the magazine *Transistor Gijutsu*, and RTKLIB Demo5 by rtkexplorer. I would like to thank Mr. Takasu, Mr. Kubo, and rtkexplorer. The code for the LAMBDA method and tropospheric delay correction was directly ported for use here.  
- For matrix operations, I use the [gonum/mat library](https://pkg.go.dev/gonum.org/v1/gonum/mat).  

# Processing
- This is an offline program that processes RINEX files. Although functions are written to handle epoch-by-epoch processing, it does not support real-time streaming input.  
- Since the code is not optimized for runtime performance, and because the ambiguity resolution step retries with different satellite pairs until a fix is achieved, execution time is about twice as long as RTKLIB.  
- The AR modes equivalent to RTKLIB’s **instantaneous** (epoch-independent) and **continuous** (Kalman filter) are implemented. There is no implementation for **fix & hold** mode.  
- The state vector of the Kalman filter includes the unknown 3D position and single-difference integer ambiguities between the base and rover for each satellite. Initially, I followed textbooks and used double-difference ambiguities as unknowns. However, since the reference satellite changes with time (as the highest-elevation satellite changes), the combinations of double-differences change as well, requiring the error covariance matrix (P) to be reset each time, which caused jumps in the positioning solution. To avoid this, I switched to the single-difference approach, as in RTKLIB.  
- After computing a float solution, if ambiguity resolution does not succeed on the first attempt, the program retries up to 10 times, removing one low-elevation satellite at a time. If it still fails, GLONASS is excluded and the process restarts from float computation. If ambiguity resolution still fails, the result remains a float solution.  
- For the **age** (time difference between base and rover observations), since this is offline processing, negative ages (future base data relative to rover) are tolerated, up to ±15 seconds.  
- Even in instantaneous mode (epoch-independent), the float solution is computed using the Kalman filter. By performing the measurement update several times (3 iterations), the result becomes equivalent to solving the nonlinear least-squares problem iteratively.  
- The time update of the Kalman filter assumes a static model without any dynamic modeling.  
- Regardless of mode, the DGPS solution is used as the initial estimate for the state vector’s 3D position before measurement update. A pure static mode (holding the previous epoch’s position) is not implemented.  
- Input RINEX files are assumed to be version 3.02 or 3.04.  

# File Structure

The main files are as follows:

| File | Description |
| ------ | ------ |
| cmd/gortk/main.go | Main process |
| calcspp.go | Single-point positioning |
| calcfloat.go | Float solution calculation |
| solveamb.go | Fixed solution calculation |

Other supporting files:

| File | Description |
| ------ | ------ |
| satpos.go | Satellite 3D position calculation |
| rinex.go | RINEX file parsing |
| obs.go | Structures and functions for observation data |
| nav.go | Structures and functions for navigation data |
| pos.go | Structures and functions for position data |
| gtime.go | Structures and functions for GPS time |
| const.go | Constant definitions |
| solvels.go | Least-squares matrix calculations (using Gonum) |
| lambda.go | LAMBDA method |
| tropo.go | Tropospheric delay correction |
| misc.go | Miscellaneous helper functions |

# Build
```sh
cd cmd/gortk
go build
```

# Usage

## Single-point positioning
```sh
gortk -p 0 rover.obs base.nav
```

## DGPS
```sh
gortk -p 1 -l "35.73101206 139.7396917 80.33" rover.obs base.obs base.nav
```
> Note: The reference station coordinates given with the `-l` option must be enclosed in double quotes.

## RTK
```sh
gortk -p 2 -l "35.73101206 139.7396917 80.33" rover.obs base.obs base.nav
```

## Show options
```sh
gortk -h
```
- Where possible, option names are consistent with RTKLIB.
- However, arguments to `-l`, `-ts`, and `-te` options must be enclosed in double quotes.
