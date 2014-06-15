FDTD法を用いたCSEMのフルウェーブインバージョン・シミュレータです。

fwi3d_cpu_posiinf  : 送受信器のポジション・計算グリッド・タイムステップ等基本となる設定ファイル
fwi3d_cpu_sigmodel : モデル内の電気伝導度の設定ファイル
fwi3d_cpu_waveform : Fictitious domainでの送信波形と実空間での送信波形を設定するファイル
fwi3d_cpu_rec2rec  : フォワード計算用のメインファイル
fwi3d_cpu_inv_fwi  : インバージョン計算用のメインファイル
fwi3d_cpu_alloli1  : フォワード計算用の変数用のヘッダー
fwi3d_cpu_alloex1  : フォワード計算用の変数用のヘッダー。

gmt_slice.rb: フォワード計算の結果を動画にはき出すために利用するconvertファイル(gmtファイルのプリ編集ファイルです）
gmt_slice.gmt : フォワード計算の結果を画像にします。最後に動画にするためにconvert等を使います
gmt_slice.sh: 上記、.rb, .gmt ファイルをまとめたものです。

Ver1.0: CPUにおけるFDTD法のシミュレーションコード作成
			Mittet(2010)によるFictitious wave domainでの高速化の組み込み
Ver1.4: Inversionの完成
Ver1.5: GPU計算の組み込み
