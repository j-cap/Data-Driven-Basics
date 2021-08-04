
"""

Creates a labeled dataset from the Bloodcell dataset in 
     C:\Users\WeberJ\Documents\Data\dataset2-master\dataset2-master


Possible labels:
    - Eosinophil
    - Lymphocyte
    - Monocyte
    - Neutrophil

"""

from pathlib import Path
import os

EOSINOPHIL = 0
LYMPHOCYTE = 1
MONOCYTE = 2
NEUTROPHIL = 3

TYPES = ["EOSINOPHIL", "LYMPHOCYTE", "MONOCYTE", "NEUTROPHIL"]

DATA_PATH = Path('C:\\Users\\WeberJ\\Documents\\Data\\dataset2-master\\dataset2-master')

for folder in os.listdir(DATA_PATH / 'images/'):
    for type_ in TYPES:
        for img in os.listdir(DATA_PATH / 'images/' / folder / type_):
            if str(folder) is "TEST":
                print("Test")
            elif str(folder) is "TRAIN":
                print("Train")
            else:
                pass
            print(img)
