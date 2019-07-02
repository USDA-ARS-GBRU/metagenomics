import json
import pandas as pd
import os

wdir = os.getcwd()

list = []
for lfile in os.listdir("rqc_data"):
    file = os.path.join(wdir, "rqc_data", lfile)
    if file.endswith('.json'):
        bfile = file
        with open(file) as js:
            js_dict = json.load(js)
            tr = js_dict['scaffoldStats1.txt']['desc']['TotalReads']
            tb = js_dict['scaffoldStats1.txt']['desc']['TotalBases']
            contam = js_dict['scaffoldStats1.txt']['desc']['ReadsMatched']
            pctcontam = js_dict['scaffoldStats1.txt']['desc']['PctReadsMatched']
            list.append((os.path.basename(bfile), tr, tb, contam, pctcontam))

cols = ['SampleID', 'TotalReads', 'TotalBases', 'Contaminants', "Percent_Contaminants"]
result = pd.DataFrame(list, columns=cols)
result.to_csv("rqc_data/parse_json.csv")
