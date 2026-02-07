#!/usr/bin/env python3
"""Create a curated, deduplicated code snapshot for GitHub.

Goal:
- Reorganize project code into a GitHub-friendly structure.
- Do NOT change runtime logic or visualization styling.

What this script does:
- Copies code/notebook files into a new folder (default: github_open_source_release/).
- Removes non-English/garbled/junk header comments (hash-style comments only).
- Removes UTF-8 BOM by re-encoding to UTF-8 without BOM.
- Deduplicates identical files by content hash.
- Writes a MANIFEST.tsv mapping source -> destination.

This is a best-effort packaging tool; the original repository is not modified.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


INCLUDE_EXTS = {".R", ".py", ".sh", ".ipynb"}


def _detect_encoding(data: bytes) -> str:
    for enc in ("utf-8-sig", "utf-8", "gb18030", "latin-1"):
        try:
            data.decode(enc)
            return enc
        except UnicodeDecodeError:
            continue
    return "utf-8"


def _decode_code_bytes(data: bytes) -> str:
    """Decode source bytes for code files (best-effort).

    Strategy:
    - Prefer UTF-8 if the byte stream is *mostly* valid UTF-8.
    - Otherwise fall back to GB18030 (common on Windows for CN locales).

    Always returns a str (uses replacement on decode failures).
    """

    # Fast path: strict UTF-8 (with BOM support)
    for enc in ("utf-8-sig", "utf-8"):
        try:
            return data.decode(enc, errors="strict")
        except UnicodeDecodeError:
            pass

    # Heuristic: if most bytes form valid UTF-8 sequences, treat as UTF-8.
    try:
        ignored = data.decode("utf-8", errors="ignore")
        valid_bytes = len(ignored.encode("utf-8"))
        ratio = valid_bytes / max(1, len(data))
    except Exception:
        ratio = 0.0

    if ratio >= 0.95:
        return data.decode("utf-8", errors="replace")

    return data.decode("gb18030", errors="replace")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _is_hash_in_string(line: str, idx: int) -> bool:
    in_s = False
    in_d = False
    esc = False
    for ch in line[:idx]:
        if esc:
            esc = False
            continue
        if ch == chr(92):
            esc = True
            continue
        if ch == "'" and not in_d:
            in_s = not in_s
        elif ch == '"' and not in_s:
            in_d = not in_d
    return in_s or in_d


def _clean_hash_comments(text: str) -> str:
    non_ascii = re.compile(r"[^\x00-\x7F]")
    many_q = re.compile(r"\?{2,}")
    junk_pats = [
        re.compile(r"biowolf", re.I),
        re.compile(r"video\s*source", re.I),
        re.compile(r"wechat|weixin|wx", re.I),
        re.compile(r"qq", re.I),
        re.compile(r"foxmail", re.I),
    ]

    def drop_comment(c: str) -> bool:
        if not c:
            return True
        if non_ascii.search(c):
            return True
        if "???" in c or many_q.search(c):
            return True
        if any(p.search(c) for p in junk_pats):
            return True
        return False

    out_lines: List[str] = []
    for line in text.splitlines(keepends=False):
        raw = line
        s = raw.lstrip()

        if s.startswith("#"):
            c = s[1:].strip()
            if drop_comment(c):
                continue
            out_lines.append(raw.rstrip())
            continue

        if "#" in raw:
            idx = raw.find("#")
            if idx != -1 and not _is_hash_in_string(raw, idx):
                c = raw[idx + 1 :].strip()
                if drop_comment(c):
                    raw = raw[:idx].rstrip()

        out_lines.append(raw.rstrip())

    # Collapse blank line runs.
    collapsed: List[str] = []
    blank_run = 0
    for l in out_lines:
        if l.strip() == "":
            blank_run += 1
            if blank_run <= 2:
                collapsed.append("")
            continue
        blank_run = 0
        collapsed.append(l)

    return "\n".join(collapsed).rstrip() + "\n"


def _fix_windows_paths_in_r_strings(text: str) -> str:
    """Fix single-backslash Windows paths inside R string literals.

    Some scripts contain paths like "D:\\...". If a file is decoded through a
    lossy codec path, a double-backslash sequence may collapse into a single
    backslash, which makes R parsing fail (e.g., "\\p" is an invalid escape).

    This function only targets string literals that look like Windows drive
    paths (e.g., "D:\\...") and only doubles *single* backslashes.
    """

    out: list[str] = []
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]
        if ch != '"':
            out.append(ch)
            i += 1
            continue

        # Enter a string literal.
        out.append(ch)
        i += 1
        start = i
        escaped = False
        buf: list[str] = []
        while i < n:
            c = text[i]
            if escaped:
                buf.append(c)
                escaped = False
                i += 1
                continue
            if c == chr(92):
                buf.append(c)
                escaped = True
                i += 1
                continue
            if c == '"':
                break
            buf.append(c)
            i += 1

        s = "".join(buf)
        if re.match(r"^[A-Za-z]:\\", s):
            # Replace single backslashes with double backslashes.
            s = re.sub(r"(?<!\\)\\(?!\\)", r"\\\\", s)
        out.append(s)

        # Closing quote (if any)
        if i < n and text[i] == '"':
            out.append('"')
            i += 1

    return "".join(out)


def _categorize(rel: Path) -> Path:
    s = rel.as_posix()

    if rel.parts[:1] == ("scripts",):
        return Path("scripts") / rel.name

    # Numbered project steps
    if s.startswith("2."):
        return Path("r/single_cell/02_standard_pipeline") / rel.name
    if s.startswith("3."):
        return Path("r/single_cell/03_annotation") / rel.name
    if s.startswith("4."):
        return Path("r/single_cell/04_lactylation_score") / rel.name
    if s.startswith("5."):
        return Path("r/single_cell/05_hdwgcna") / rel.name
    if s.startswith("6.2."):
        return Path("r/bulk/06_degs_and_cor") / rel.name
    if s.startswith("6."):
        return Path("r/bulk/06_consensus_clustering") / rel.name
    if s.startswith("7."):
        return Path("r/bulk/07_download") / rel.name

    if s.startswith("8."):
        if "14.MLdata" in s:
            return Path("r/ml/benchmark/14_mldata") / rel.name
        if "15.ML" in s:
            return Path("r/ml/benchmark/15_ml") / rel.name
        if "模型ROC" in s or "ROC" in s:
            return Path("r/ml/benchmark/roc") / rel.name
        if "混淆矩阵" in s:
            return Path("r/ml/benchmark/confusion_matrix") / rel.name
        if "随机森林" in s:
            return Path("r/ml/benchmark/random_forest") / rel.name
        return Path("r/ml/benchmark/misc") / rel.name

    if s.startswith("9."):
        return Path("r/ml/gmm_lr") / rel.name
    if s.startswith("10."):
        return Path("r/ml/mbds_scoring") / rel.name
    if s.startswith("11.深度学习"):
        return Path("r/ml/deep_learning") / rel.name
    if s.startswith("11.ROC"):
        return Path("r/ml/roc") / rel.name
    if s.startswith("12."):
        return Path("r/figures/boxplots") / rel.name
    if s.startswith("13."):
        return Path("r/pathway_quantification") / rel.name
    if s.startswith("14."):
        return Path("r/single_cell/cell_communication") / rel.name
    if s.startswith("16."):
        return Path("r/single_cell/trajectory") / rel.name
    if s.startswith("17."):
        return Path("r/genes/rcc2_function") / rel.name

    if "上皮细胞拟时序" in s:
        return Path("r/single_cell/epithelial_pseudotime") / rel.name
    if "数据补充测试" in s:
        return Path("r/extras") / rel.name

    if rel.suffix == ".R":
        return Path("r/misc") / rel.name
    if rel.suffix == ".ipynb":
        return Path("notebooks") / rel.name
    return rel


def _rename_numeric_r(dest_rel: Path) -> Path:
    """Rename numeric R scripts (e.g., 1.R/2.R) to descriptive names.

    Only applies to a small, curated set of known destinations.
    """
    if dest_rel.suffix != ".R":
        return dest_rel

    name = dest_rel.name
    parent = dest_rel.parent.as_posix()

    mapping: Dict[Tuple[str, str], str] = {
        ("r/ml/benchmark/14_mldata", "1.R"): "build_train_test_mldata_combat.R",
        ("r/ml/benchmark/14_mldata", "2.R"): "build_train_test_mldata_no_combat_quiet_labels.R",
        ("r/ml/benchmark/15_ml", "1.R"): "run_ml_benchmark_113_combinations.R",
        ("r/ml/benchmark/random_forest", "1.R"): "random_forest_feature_importance.R",
        ("r/ml/mbds_scoring", "1.R"): "mbds_scoring_analysis.R",
        ("r/single_cell/trajectory", "1.R"): "trajectory_aec_monocle2_analysis.R",
        ("r/single_cell/cell_communication", "1.R"): "cellchat_retinoic_acid_high_low_analysis.R",
        ("r/single_cell/epithelial_pseudotime", "1.R"): "epithelial_pseudotime_monocle2_analysis.R",
        ("r/pathway_quantification", "1.R"): "quantify_metabolic_pathway_scores.R",
        ("r/figures/boxplots", "1.R"): "plot_gene_expression_violin_boxplots.R",
        ("r/extras", "1.R"): "supplementary_geo_processing.R",
        ("r/genes/rcc2_function", "1.R"): "rcc2_gsea_and_correlation.R",
        ("r/genes/rcc2_function", "2.R"): "rcc2_expression_group_comparison.R",
    }

    new_name = mapping.get((parent, name))
    if not new_name:
        return dest_rel
    return dest_rel.with_name(new_name)


def _should_include(
    rel: Path,
    *,
    exclude_segments: Iterable[str],
    exclude_if_contains: Iterable[str],
    include_exts: Iterable[str],
) -> bool:
    if rel.suffix not in set(include_exts):
        return False
    if any(seg in set(exclude_segments) for seg in rel.parts):
        return False
    rel_s = str(rel)
    for token in exclude_if_contains:
        if token in rel_s:
            return False
    return True


def run(*, out_dir: Path, root: Path) -> None:
    exclude_segments = {
        "github_open_source",
        "github_open_source_manuscript_order",
        "github_open_source_release",
        str(out_dir.name),
        "ml_summary_tables",
        ".git",
        "__pycache__",
        ".ipynb_checkpoints",
    }
    exclude_if_contains = [
        "废弃",
        "ANCA-GN_transcriptomics-main",
    ]

    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    seen_hash_to_dest: Dict[str, str] = {}
    manifest: List[str] = [
        "\t".join(
            [
                "source_relpath",
                "dest_relpath",
                "sha256",
                "status",
                "canonical_dest_relpath",
            ]
        )
    ]

    included = 0
    copied = 0
    deduped = 0

    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(root)
        if not _should_include(
            rel,
            exclude_segments=exclude_segments,
            exclude_if_contains=exclude_if_contains,
            include_exts=INCLUDE_EXTS,
        ):
            continue

        included += 1
        dest_rel = _categorize(rel)
        dest_rel = _rename_numeric_r(dest_rel)
        dest = out_dir / dest_rel

        data = path.read_bytes()
        if path.suffix == ".ipynb":
            cleaned = data
        else:
            txt = _decode_code_bytes(data)
            txt = txt.replace("\r\n", "\n").replace("\r", "\n")
            txt = _clean_hash_comments(txt)
            if path.suffix == ".R":
                txt = _fix_windows_paths_in_r_strings(txt)
            cleaned = txt.encode("utf-8")

        h = _sha256(cleaned)
        if h in seen_hash_to_dest:
            deduped += 1
            canonical = seen_hash_to_dest[h]
            manifest.append(
                "\t".join([str(rel), str(dest_rel), h, "dedup_skipped", canonical])
            )
            continue

        if dest.exists():
            # Collision with different content: append short hash.
            dest = dest.with_name(f"{dest.stem}__{h[:8]}{dest.suffix}")
            dest_rel = dest.relative_to(out_dir)

        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(cleaned)
        seen_hash_to_dest[h] = str(dest_rel)
        copied += 1
        manifest.append("\t".join([str(rel), str(dest_rel), h, "copied", str(dest_rel)]))

    (out_dir / "MANIFEST.tsv").write_text("\n".join(manifest) + "\n", encoding="utf-8")

    print(f"Created: {out_dir}")
    print(f"Included candidates: {included}")
    print(f"Copied unique files: {copied}")
    print(f"Dedup skipped files: {deduped}")


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Create a curated, deduplicated code snapshot for GitHub.",
    )
    p.add_argument(
        "--out_dir",
        type=Path,
        default=Path("github_open_source_release"),
        help="Output directory to create.",
    )
    args = p.parse_args(argv)
    run(out_dir=args.out_dir, root=Path("."))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
