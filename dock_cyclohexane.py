#!/usr/bin/env python3
"""扫描 case 目录，把每个齐全的 case 与给定小分子对接。

用法：

    python dock_cyclohexane.py
    python dock_cyclohexane.py scpdb_cases

默认读取 ``scpdb_cases`` 下的每个子文件夹。某个 case 缺少标准文件时，
把缺失项写入 ``docking/missing_files.txt``。``protein.mol2`` 或
``ligand.mol2`` 缺失时跳过该 case 的对接。
"""

from __future__ import annotations

import csv
import os
import re
from datetime import datetime
import subprocess
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

CASES_DIR = Path("scpdb_cases")
DOCKING_DIR = Path("docking")
MISSING_LOG = DOCKING_DIR / "missing_files.txt"
BINDING_TABLE = DOCKING_DIR / "binding_energies.csv"
BINDING_FIELDS = ("配体", "靶点", "蛋白", "case", "结合能_kcal_mol", "完成时间", "结果目录")
# 靶点用常用简称，蛋白用全称。未知 PDB 只保留编号。
TARGETS = {
    "1oir": ("CDK2", "细胞周期蛋白依赖性激酶 2"),
    "1hsg": ("HIV-1 蛋白酶", "HIV-1 蛋白酶"),
    "1iep": ("ABL1", "酪氨酸蛋白激酶 ABL1"),
    "2pmd": ("翻译起始因子 IF2γ", "翻译起始因子 IF2 γ 亚基"),
}
VINA = Path.home() / "AppData/Local/AutoDock-Vina/vina.exe"
SMILES = "COC1=CC(=CC(=C1O)OC)[C@@H]2[C@H]3CO[C@@H]([C@H]3CO2)C4=CC(=C(C=C4)O)OC"
LIGAND_NAME = "medioresinol"
SEED = 1
CASE_FILES = (
    "protein.mol2",
    "ligand.mol2",
    "ligand.sdf",
    "site.mol2",
    "cavity6.mol2",
    "ints_C.mol2",
    "ints_M.mol2",
    "IFP.txt",
)
REQUIRED_FOR_DOCKING = ("protein.mol2", "ligand.mol2")
STANDARD = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "KCX",
}


def parse_substructures(lines: list[str]) -> dict[int, tuple[str, str, int]]:
    residues: dict[int, tuple[str, str, int]] = {}
    section = None
    for raw in lines:
        if raw.startswith("@<TRIPOS>"):
            section = raw.split(">", 1)[1]
            continue
        if section != "SUBSTRUCTURE":
            continue
        parts = raw.split()
        if len(parts) < 7 or not parts[0].isdigit() or parts[3] not in {"RESIDUE", "GROUP"}:
            continue
        if len(parts) >= 12 and parts[11].isdigit():
            resseq = int(parts[11])
        else:
            tail = parts[-1]
            match = re.fullmatch(r"[A-Za-z](\d+)", tail)
            resseq = int(match.group(1) if match else tail)
        residues[int(parts[0])] = (parts[6].upper(), parts[5], resseq)
    return residues


def parse_mol2(path: Path) -> list[dict]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    residues = parse_substructures(lines)
    atoms: list[dict] = []
    section = None
    for raw in lines:
        if raw.startswith("@<TRIPOS>"):
            section = raw.split(">", 1)[1]
            continue
        if section != "ATOM":
            continue
        parts = raw.split()
        if len(parts) < 8:
            continue
        sybyl = parts[5]
        element = "H" if sybyl == "H" else sybyl.split(".", 1)[0]
        if element == "H":
            continue
        resname, chain, resseq = residues[int(parts[6])]
        if resname not in STANDARD:
            continue
        atoms.append(
            {
                "name": parts[1],
                "x": float(parts[2]),
                "y": float(parts[3]),
                "z": float(parts[4]),
                "element": element.capitalize() if len(element) > 1 else element,
                "resname": resname,
                "resseq": resseq,
                "chain": chain,
            }
        )
    return atoms


def ligand_box(path: Path) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    points: list[tuple[float, float, float]] = []
    section = None
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith("@<TRIPOS>"):
            section = raw.split(">", 1)[1]
            continue
        if section != "ATOM":
            continue
        parts = raw.split()
        if len(parts) < 6 or parts[5] == "H":
            continue
        points.append((float(parts[2]), float(parts[3]), float(parts[4])))
    if not points:
        raise RuntimeError(f"{path} 里没有配体重原子")
    axes = list(zip(*points))
    count = len(points)
    center = tuple(sum(values) / count for values in axes)
    size = tuple(round(max(values) - min(values) + 10.0, 1) for values in axes)
    return center, size


def collapse_unresolved_sidechains(atoms: list[dict]) -> list[str]:
    """把只解析到主链或 CB 的残基改成 GLY/ALA，避免残缺侧链无法匹配模板。"""
    backbone = {"N", "CA", "C", "O", "OXT"}
    grouped: dict[tuple[str, int, str], list[dict]] = {}
    for atom in atoms:
        grouped.setdefault((atom["chain"], atom["resseq"], atom["resname"]), []).append(atom)
    changed = []
    for (chain, resseq, resname), group in grouped.items():
        names = {atom["name"] for atom in group}
        if names - backbone - {"CB"}:
            continue
        new_name = "ALA" if "CB" in names else "GLY"
        if new_name == resname:
            continue
        for atom in group:
            atom["resname"] = new_name
        changed.append(f"{resname}{resseq}->{new_name}")
    return changed


def pdb_atom_name(name: str) -> str:
    if len(name) >= 4:
        return name[:4]
    return f" {name:<3}"


def write_pdb(atoms: list[dict], path: Path) -> None:
    lines = []
    for serial, atom in enumerate(atoms, start=1):
        record = "ATOM" if atom["resname"] not in {"HOH", "WAT"} else "HETATM"
        name = pdb_atom_name(atom["name"])
        element = atom["element"][:2].rjust(2)
        lines.append(
            f"{record:<6}{serial:5d} {name} {atom['resname']:>3} {atom['chain']}"
            f"{atom['resseq']:4d}    "
            f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}  1.00  0.00          {element}"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def prepare_ligand(path: Path, out_dir: Path) -> None:
    sys.path.insert(0, str(out_dir / "_pydeps"))
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    if mol is None:
        raise RuntimeError(f"无法解析 SMILES: {SMILES}")
    embedded = False
    for offset in range(8):
        params = AllChem.ETKDGv3()
        params.randomSeed = SEED + offset
        params.enforceChirality = True
        if AllChem.EmbedMolecule(mol, params) == 0:
            embedded = True
            break
    if not embedded:
        raise RuntimeError("配体三维构象生成失败")
    status = AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    if status == -1:
        raise RuntimeError("配体 MMFF 优化失败")
    Chem.MolToMolFile(mol, str(path.with_suffix(".mol")))
    setups = MoleculePreparation().prepare(mol)
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(error)
    path.write_text(pdbqt, encoding="utf-8")


def list_case_dirs(cases_dir: Path) -> list[Path]:
    if not cases_dir.is_dir():
        raise RuntimeError(f"case 目录不存在: {cases_dir}")
    return sorted(path for path in cases_dir.iterdir() if path.is_dir())


def missing_files(case_dir: Path) -> list[str]:
    missing = []
    for name in CASE_FILES:
        path = case_dir / name
        if not path.is_file():
            missing.append(name)
        elif path.stat().st_size == 0:
            missing.append(f"{name}（空文件）")
    return missing


def can_dock(missing: list[str]) -> bool:
    lacking = {item.split("（", 1)[0] for item in missing}
    return not any(name in lacking for name in REQUIRED_FOR_DOCKING)


def write_missing_report(cases_dir: Path, rows: list[tuple[str, list[str]]], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"检查目录: {cases_dir.resolve()}", ""]
    incomplete = [(case_id, missing) for case_id, missing in rows if missing]
    if not incomplete:
        lines.append("没有缺少文件的 case。")
    else:
        lines.append("缺少文件的 case：")
        for case_id, missing in incomplete:
            status = "已跳过对接" if not can_dock(missing) else "仍对接（蛋白和配体文件可用）"
            lines.append(f"{case_id}\t{status}\t缺少: {', '.join(missing)}")
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def dock_case(case_dir: Path) -> None:
    case_id = case_dir.name
    out_dir = DOCKING_DIR / f"{case_id}_{LIGAND_NAME}"
    out_dir.mkdir(parents=True, exist_ok=True)
    center, size = ligand_box(case_dir / "ligand.mol2")
    print(
        f"\n== {case_id} ==\n"
        f"盒子中心 {' '.join(f'{value:.3f}' for value in center)}  "
        f"尺寸 {' '.join(f'{value:.1f}' for value in size)}"
    )
    atoms = parse_mol2(case_dir / "protein.mol2")
    truncated = collapse_unresolved_sidechains(atoms)
    print("截断侧链改为 ALA/GLY:", ", ".join(truncated) if truncated else "无")
    unusual = sorted(
        {f"{atom['chain']}:{atom['resname']}{atom['resseq']}" for atom in atoms if atom["resname"] not in STANDARD - {"KCX"}}
    )
    print("非标准残基:", ", ".join(unusual) if unusual else "无")
    protein_pdb = out_dir / "protein.pdb"
    write_pdb(atoms, protein_pdb)
    print(f"重原子数: {len(atoms)} -> {protein_pdb}")

    kcx_template = out_dir / "kcx_template.json"
    kcx_template.write_text(
        """{
  "ambiguous": {"KCX": ["KCX"]},
  "residue_templates": {
    "KCX": {
      "smiles": "C(=O)C([H])(N[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=O)O[H]",
      "atom_name": ["C", "O", "CA", "HA", "N", "H", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "CE", "HE2", "HE3", "NZ", "HZ", "CX", "OQ2", "OQ1", "HQ1"],
      "link_labels": {"0": "C-term", "4": "N-term"}
    }
  }
}
""",
        encoding="utf-8",
    )
    gemmi_stub = out_dir / "_pydeps"
    gemmi_stub.mkdir(exist_ok=True)
    (gemmi_stub / "gemmi.py").write_text(
        "class _Cif:\n"
        "    @staticmethod\n"
        "    def read(source):\n"
        "        raise RuntimeError('gemmi is not installed; KCX should use the local template')\n"
        "cif = _Cif()\n",
        encoding="utf-8",
    )
    meeko_script = Path(sys.prefix) / "Scripts" / "mk_prepare_receptor.exe"
    if not meeko_script.exists():
        meeko_script = Path(sys.executable).with_name("mk_prepare_receptor.exe")
    receptor_cmd = [
        str(meeko_script),
        "--read_pdb",
        str(protein_pdb),
        "-o",
        str(out_dir / "receptor"),
        "-p",
        "--box_center",
        *(f"{value:.3f}" for value in center),
        "--box_size",
        *(f"{value:.1f}" for value in size),
        "--delete_bad_res_from_box_radius",
        "2",
        "--add_templates",
        str(kcx_template),
    ]
    print("准备受体...")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(gemmi_stub) + ((";" + env["PYTHONPATH"]) if env.get("PYTHONPATH") else "")
    subprocess.run(receptor_cmd, check=True, env=env)

    ligand_pdbqt = out_dir / f"{LIGAND_NAME}.pdbqt"
    print("准备配体...")
    prepare_ligand(ligand_pdbqt, out_dir)

    receptor_pdbqt = out_dir / "receptor.pdbqt"
    if not receptor_pdbqt.exists():
        candidates = sorted(out_dir.glob("receptor*.pdbqt"))
        if len(candidates) != 1:
            raise RuntimeError(f"找不到唯一的受体 PDBQT: {candidates}")
        receptor_pdbqt = candidates[0]

    log_path = out_dir / "vina.log"
    docked_path = out_dir / "docked.pdbqt"
    vina_cmd = [
        str(VINA),
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        f"{center[0]:.3f}",
        "--center_y",
        f"{center[1]:.3f}",
        "--center_z",
        f"{center[2]:.3f}",
        "--size_x",
        f"{size[0]:.1f}",
        "--size_y",
        f"{size[1]:.1f}",
        "--size_z",
        f"{size[2]:.1f}",
        "--exhaustiveness",
        "16",
        "--seed",
        str(SEED),
        "--out",
        str(docked_path),
    ]
    print("开始对接...")
    result = subprocess.run(vina_cmd, check=True, capture_output=True, text=True)
    log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    print(result.stdout, flush=True)
    affinity = best_affinity(result.stdout)
    record_binding(case_id, LIGAND_NAME, affinity, out_dir)
    target, protein = target_names(case_id)
    print(f"已记录: {target} / {protein} / {LIGAND_NAME} = {affinity:.3f} kcal/mol")


def target_names(case_id: str) -> tuple[str, str]:
    pdb_id = case_id.split("_", 1)[0]
    return TARGETS.get(pdb_id, (pdb_id.upper(), pdb_id.upper()))


def best_affinity(vina_text: str) -> float:
    for line in vina_text.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                return float(parts[1])
            except ValueError:
                continue
    raise RuntimeError("Vina 输出里没有最佳结合能")


def load_binding_rows() -> list[dict[str, str]]:
    if not BINDING_TABLE.is_file():
        return []
    with BINDING_TABLE.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_binding_rows(rows: list[dict[str, str]]) -> None:
    BINDING_TABLE.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda row: float(row["结合能_kcal_mol"]))
    with BINDING_TABLE.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=BINDING_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def record_binding(case_id: str, ligand: str, affinity: float, out_dir: Path, finished_at: str | None = None) -> None:
    target, protein = target_names(case_id)
    row = {
        "配体": ligand,
        "靶点": target,
        "蛋白": protein,
        "case": case_id,
        "结合能_kcal_mol": f"{affinity:.3f}",
        "完成时间": finished_at or datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "结果目录": str(out_dir).replace("\\", "/"),
    }
    rows = [item for item in load_binding_rows() if not (item["配体"] == ligand and item["case"] == case_id)]
    rows.append(row)
    write_binding_rows(rows)


def backfill_binding_table() -> int:
    """把已有 vina.log 补进表格。同一配体和 case 已有记录时不覆盖。"""
    existing = {(row["配体"], row["case"]) for row in load_binding_rows()}
    added = 0
    if not DOCKING_DIR.is_dir():
        return 0
    for log_path in sorted(DOCKING_DIR.glob("*/vina.log")):
        folder = log_path.parent.name
        ligand = None
        case_id = None
        for name in ("medioresinol", "cyclohexane"):
            suffix = f"_{name}"
            if folder.endswith(suffix):
                ligand = name
                case_id = folder[: -len(suffix)]
                break
        if ligand is None or case_id is None:
            continue
        if case_id == "1oir":
            case_id = "1oir_1"
        if (ligand, case_id) in existing:
            continue
        affinity = best_affinity(log_path.read_text(encoding="utf-8", errors="replace"))
        finished = datetime.fromtimestamp(log_path.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
        record_binding(case_id, ligand, affinity, log_path.parent, finished)
        existing.add((ligand, case_id))
        added += 1
    return added


def main() -> int:
    added = backfill_binding_table()
    if added:
        print(f"已从已有结果补入 {added} 行到 {BINDING_TABLE}")
    cases_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else CASES_DIR
    case_dirs = list_case_dirs(cases_dir)
    rows = [(case_dir.name, missing_files(case_dir)) for case_dir in case_dirs]
    write_missing_report(cases_dir, rows, MISSING_LOG)
    print(f"文件检查已写入 {MISSING_LOG}")
    for case_id, missing in rows:
        if missing:
            print(f"{case_id} 缺少: {', '.join(missing)}")

    failed = []
    for case_dir, (_case_id, missing) in zip(case_dirs, rows):
        if not can_dock(missing):
            print(f"跳过 {case_dir.name}：缺少对接所需的 protein.mol2 或 ligand.mol2")
            continue
        try:
            dock_case(case_dir)
        except Exception as exc:
            failed.append(case_dir.name)
            print(f"{case_dir.name} 失败: {exc}", file=sys.stderr)
    if failed:
        print("未完成:", ", ".join(failed), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
