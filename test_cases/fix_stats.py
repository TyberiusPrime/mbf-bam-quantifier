from pathlib import Path
import shutil

for p in Path('test_cases').glob('**/counts.tsv.stats.tsv'):
    if 'actual' in str(p):
        print(p)
        target = p.parent.parent.parent / "output/counts.tsv.stats.tsv"
        print(target)
        shutil.copy(p, target)
