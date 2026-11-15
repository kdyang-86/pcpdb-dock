#!/usr/bin/env python3
"""下载 sc-PDB 中的一个 case。

sc-PDB（http://bioinfo-pharma.u-strasbg.fr/scPDB/）是从 PDB 整理出的
可成药结合位点数据库。一个 case 对应一个蛋白–配体结合位点，编号形如
``1oir_1``（4 位 PDB ID + 配体序号）。

每个 case 会保存到单独目录，包含：

    protein.mol2    蛋白（含配体附近的离子和辅因子）
    ligand.mol2     配体
    ligand.sdf      配体（SDF）
    site.mol2       结合位点
    cavity6.mol2    空腔（VolSite）
    ints_C.mol2     相互作用示意图（以配体为中心）
    ints_M.mol2     相互作用示意图（另一种表示）
    IFP.txt         相互作用指纹

用法：

    python download_case.py 1oir_1
    python download_case.py 1oir -o data
"""

from __future__ import annotations

import argparse
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

BASE_URL = "http://bioinfo-pharma.u-strasbg.fr/scPDB"
# 站点当前发布目录。条目页里的路径是 ressources/2016/entries/<case_id>/。
RELEASE = "2016"

FILES = (
    "protein.mol2",
    "ligand.mol2",
    "ligand.sdf",
    "site.mol2",
    "cavity6.mol2",
    "ints_C.mol2",
    "ints_M.mol2",
    "IFP.txt",
)

CASE_ID_RE = re.compile(r"^[0-9][a-z0-9]{3}_[1-9][0-9]*$")
PDB_ID_RE = re.compile(r"^[0-9][a-z0-9]{3}$")

USER_AGENT = "scpdb-download-case/1.0"
TIMEOUT = 60
RETRIES = 3


class DownloadError(RuntimeError):
    pass


def case_file_url(case_id: str, filename: str) -> str:
    return f"{BASE_URL}/ressources/{RELEASE}/entries/{case_id}/{filename}"


def request(url: str, method: str = "GET") -> urllib.request.Request:
    return urllib.request.Request(url, method=method, headers={"User-Agent": USER_AGENT})


def fetch(url: str, method: str = "GET") -> tuple[int, bytes]:
    last_error: Exception | None = None
    for attempt in range(1, RETRIES + 1):
        try:
            with urllib.request.urlopen(request(url, method), timeout=TIMEOUT) as response:
                status = getattr(response, "status", 200)
                body = b"" if method == "HEAD" else response.read()
                return status, body
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return 404, b""
            last_error = exc
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc
        if attempt < RETRIES:
            time.sleep(attempt)
    raise DownloadError(f"请求失败: {url} ({last_error})")


def case_exists(case_id: str) -> bool:
    status, _ = fetch(case_file_url(case_id, "protein.mol2"), method="HEAD")
    return status == 200


def resolve_cases(query: str) -> list[str]:
    """把 ``1oir_1`` 或 ``1oir`` 解析成 case 编号列表。"""
    query = query.strip().lower()
    if CASE_ID_RE.fullmatch(query):
        if not case_exists(query):
            raise DownloadError(f"sc-PDB 中没有 case {query}")
        return [query]
    if not PDB_ID_RE.fullmatch(query):
        raise DownloadError(
            f"编号格式不对: {query}。请使用 PDB ID（如 1oir）或 case ID（如 1oir_1）。"
        )

    found: list[str] = []
    index = 1
    while index <= 50:
        case_id = f"{query}_{index}"
        if not case_exists(case_id):
            break
        found.append(case_id)
        index += 1
    if not found:
        raise DownloadError(f"sc-PDB 中没有 PDB {query} 的 case")
    return found


def download_case(case_id: str, out_dir: Path) -> Path:
    dest = out_dir / case_id
    dest.mkdir(parents=True, exist_ok=True)
    print(f"下载 {case_id} -> {dest}")
    for filename in FILES:
        url = case_file_url(case_id, filename)
        status, body = fetch(url)
        if status != 200 or not body:
            raise DownloadError(f"{case_id}/{filename} 下载失败 (HTTP {status})")
        path = dest / filename
        path.write_bytes(body)
        print(f"  {filename:16} {len(body):8} bytes")
    return dest


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="下载 sc-PDB 的一个（或同一 PDB 下的全部）结合位点 case。"
    )
    parser.add_argument(
        "case",
        help="case 编号（1oir_1）或 PDB ID（1oir）。只给 PDB ID 时会下载该结构的全部 case。",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        default=Path("scpdb_cases"),
        help="输出目录，默认 scpdb_cases/",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        case_ids = resolve_cases(args.case)
        for case_id in case_ids:
            download_case(case_id, args.out)
    except DownloadError as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1
    print(f"完成，共 {len(case_ids)} 个 case。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
