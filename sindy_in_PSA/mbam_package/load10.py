import numpy as np, pandas as pd
B='/mnt/user-data/attachments/670b9ebf-6058-52c2-944f-f5b88ac16043/'
GF={1:'098952f1-G1_Fag_10s.csv',2:'f37729d0-G2_Fag_10s.csv',3:'723a8281-G3_Fag_10s.csv'}
BF={2:'4699d727-IVPQ_B2_H.csv',7:'63d82e4d-IVPQ_B7_H.csv',8:'61d5b1d3-IVPQ_B8_H.csv'}
G={k:pd.read_csv(B+v) for k,v in GF.items()}
BU={k:pd.read_csv(B+v) for k,v in BF.items()}
