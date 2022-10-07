import glob
import pandas as pd

for output_table in glob.glob('out/processed/*'):
    read_table = pd.read_csv(output_table)
    display(read_table[0])