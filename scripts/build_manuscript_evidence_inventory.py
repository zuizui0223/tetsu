"""Inventory a fixed Git snapshot; do not infer ecological support from filenames."""
from pathlib import Path
import csv
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
BASE = "da8cc065a61cd4c155b1a6a291fcf91c498a79e8"
OUT = ROOT / "docs" / "evidence_inventory"


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args])


def family(path):
    name = Path(path).name
    versions = [int(v) for v in re.findall(r"_v(\d+)", name)]
    v = versions[0] if versions else None
    if path.startswith("research/"):
        return "F12 国外頻度・周辺比較・広域追試"
    if v is not None:
        bands = [
            (0, 16, "F01 写真・フェノロジー・管理・衛星の測定可能性"),
            (17, 73, "F02 GBIF広域分布・ネットワーク・記録過程"),
            (74, 74, "F03 広域種置換・入れ子性"),
            (75, 97, "F04 列車頻度・貨物流動と広域植物記録"),
            (98, 101, "F05 MIREN道路・内陸の組成と外来種"),
            (102, 110, "F06 鉄道微生息地・生活史・散布形質"),
            (113, 127, "F07 鉄道対道路の形質・頻度重み・研究依存性"),
            (128, 135, "F08 種子輸送・生存・野生化の段階別証拠"),
            (141, 146, "F09 種同一性・機能群・経路の識別"),
            (147, 149, "F10 既往メタ解析・適応進化"),
            (150, 163, "F09 種同一性・機能群・経路の識別"),
            (164, 171, "F11 種内形質・写真・ITVの識別"),
            (172, 173, "F13 種同一性と機能群の留保予測"),
            (174, 176, "F11 種内形質・写真・ITVの識別"),
            (180, 188, "F12 国外頻度・周辺比較・広域追試"),
            (190, 197, "F14 地域プール・参照群集選択"),
            (198, 200, "F13 種同一性と機能群の留保予測"),
            (210, 211, "F13 種同一性と機能群の留保予測"),
            (221, 222, "F11 種内形質・写真・ITVの識別"),
            (223, 233, "F15 列車頻度・線区・機能幅と地点間差"),
            (236, 244, "F14 地域プール・参照群集選択"),
            (245, 255, "F12 国外頻度・周辺比較・広域追試"),
            (256, 268, "F16 局所種プール・被度・存続と地域留保"),
            (269, 272, "F17 局所共有・種置換・在来外来"),
            (273, 284, "F18 局所結合と遠隔路線共有の対比"),
            (285, 285, "F19 草丈構造・種子重・帰無制約"),
            (286, 289, "F20 鉄道対非鉄道・外部対照の支持域"),
            (291, 298, "F19 草丈構造・種子重・帰無制約"),
            (299, 302, "F21 標高勾配・地点道路地域留保予測"),
            (303, 303, "F22 欧州外鉄道・同面積・形質定義"),
        ]
        for lo, hi, label in bands:
            if lo <= v <= hi:
                return label
    if name in {"README.md", "claim_evidence_crosswalk.md", "sources.md",
                "tetsu_manuscript_positive_evidence_synthesis_ja.md"}:
        return "F00 既存統合・出典案内（年代により更新）"
    if name == "rail_environment_traits_protocol.md":
        return "F06 鉄道微生息地・生活史・散布形質"
    return "F01 写真・フェノロジー・管理・衛星の測定可能性"


def write_csv(path, rows, columns):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True)
    paths = git("ls-tree", "-r", "--name-only", BASE).decode().splitlines()
    # Read blobs from one immutable Git tree, including binary and implementation files.
    request = "".join(f"{BASE}:{p}\n" for p in paths).encode()
    raw = subprocess.run(["git", "-C", str(ROOT), "cat-file", "--batch"],
                         input=request, stdout=subprocess.PIPE, check=True).stdout
    cursor = 0
    manifest, reports = [], []
    for path in paths:
        end = raw.index(b"\n", cursor)
        _, kind, size = raw[cursor:end].decode().split()
        size = int(size)
        blob = raw[end + 1:end + 1 + size]
        cursor = end + 2 + size
        assert kind == "blob" and len(blob) == size
        is_report = path.endswith(".md")
        role = "報告・設計・既存総説" if is_report else (
            "集計・文献転記・図等" if path.startswith(("data/", "research/")) else
            "実装・テスト・再現手順等")
        manifest.append(dict(path=path, bytes=size, sha256=hashlib.sha256(blob).hexdigest(), role=role))
        if not is_report:
            continue
        text = blob.decode("utf-8-sig")
        headings = re.findall(r"^#{1,3}\s+(.+)", text, re.M)
        # Excerpts are navigational evidence, not automated support judgements.
        candidates = [line.strip() for line in text.splitlines()
                      if line.strip() and not line.startswith(("#", "|", "```", "- "))]
        reports.append(dict(path=path, family=family(path),
                            title=headings[0] if headings else "", headings=" / ".join(headings),
                            opening_excerpt=" ".join(candidates[:2])[:700],
                            review="本文・結果表の確認対象。肯定/否定の機械分類なし",
                            treatment="同じ原研究の版・感度分析は独立証拠に加算しない"))
    write_csv(OUT / "repository_files.csv", manifest, ["path", "bytes", "sha256", "role"])
    write_csv(OUT / "all_reports.csv", reports,
              ["path", "family", "title", "headings", "opening_excerpt", "review", "treatment"])
    counts = {label: sum(r["family"] == label for r in reports)
              for label in sorted({r["family"] for r in reports})}
    summary = dict(snapshot=BASE, tracked_files=len(paths), markdown_reports=len(reports),
                   docs_reports=sum(r["path"].startswith("docs/") for r in reports),
                   research_reports=sum(r["path"].startswith("research/") for r in reports),
                   families=counts,
                   scope="Snapshot inventory, not a rerun of every analysis or a census of external data.")
    (OUT / "coverage.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    lines = ["# 全報告の索引", "", f"対象：`{BASE}`。{len(paths)}追跡ファイル、{len(reports)} Markdown。",
             "", "ファイルの存在・本文・見出しを全件取得した目録。実装・図・集計も別CSVへ全件収録した。"
             "全数値の原データ再計算を意味しない。版数・原著数・独立野外調査数を区別する。", ""]
    for label in counts:
        lines.extend([f"## {label}（{counts[label]}文書）", ""])
        for row in reports:
            if row["family"] == label:
                lines.append(f"- [{row['path']}](../../{row['path']}) — {row['title']}")
        lines.append("")
    (OUT / "all_reports_index_ja.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
