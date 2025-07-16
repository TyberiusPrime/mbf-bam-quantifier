from pathlib import Path
import shutil


for fn in 'counts.tsv.stats.tsv', 'matrix.mtx.stats.tsv':

    for p in Path('test_cases').glob(f'**/{fn}'):
        if 'actual' in str(p):
            print(p)
            target = p.parent.parent.parent / f"{p.parent.parent.name}/{fn}"
            if target.exists():
                shutil.copy(p, target)
