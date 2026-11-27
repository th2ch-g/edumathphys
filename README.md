# EduMathPhys

大学数学・物理を学部基礎から大学院まで日本語で体系的に独学するための静的サイト。

- 定義・定理・証明の番号付き環境と科目横断の相互参照（ホバープレビュー付き）
- 解答付き演習問題と、その場で採点できる確認テスト
- 学習進捗の記録（ブラウザ内保存）、科目間の前提関係グラフ、用語索引（日英）、記号表
- 全科目の章・節レベルの目次と執筆状況

## 開発

```sh
npm ci
npm run dev       # http://localhost:4321
npm run check     # content consistency check
npm run build     # static build into dist/
```

GitHub Pages のプロジェクトサイトとして配信する場合は `SITE_BASE=/<repo>/` と `SITE_URL=https://<user>.github.io` を指定してビルドする。

執筆規約と記法は [docs/STYLE.md](docs/STYLE.md) を参照。

## ライセンス

- プログラム: [MIT](LICENSE)
- 本文（`src/content/`, `curriculum/`）: [CC BY-SA 4.0](LICENSE-CONTENT)
