from pathlib import Path
import shutil

path = Path('test_cases')
if Path('.').absolute().name == 'test_cases':
    path = Path('.')

for fn in 'counts.tsv.stats.tsv', 'matrix.mtx.stats.tsv':

    for p in path.glob(f'**/{fn}'):
        if 'actual' in str(p):
            print(p)
            target = p.parent.parent.parent / f"{p.parent.name}/{fn}"
            if target.exists():
                shutil.copy(p, target)
