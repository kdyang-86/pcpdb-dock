# scpdb

从 sc-PDB 下载蛋白结合位点，并用 AutoDock Vina 对接小分子。

## 环境

- Python 3
- RDKit、Meeko、NumPy、SciPy
- AutoDock Vina 1.2.7，命令名为 ina

## 下载一个 case

`powershell
python download_case.py 1oir_1
`

文件保存在 scpdb_cases/<case编号>/。

## 对接

`powershell
python dock_cyclohexane.py
`

脚本会扫描 scpdb_cases 下的每个 case。缺少文件时写入 docking/missing_files.txt。每个 case 对接完成后，最佳结合能写入 docking/binding_energies.csv。

当前配体是中脂素（medioresinol）。已包含 1oir_1、1hsg_1、1iep_1、2pmd_1 的下载数据和对接结果。
