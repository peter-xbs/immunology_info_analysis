# _*_ coding:utf-8 _*_

import requests

# data = requests.get('https://www.frontiersin.org/articles/file/downloadfile/317773_supplementary-materials_presentations_1_zip/octet-stream/Presentation%201.ZIP/1/317773')
#
# with open('tmp.zip', 'wb') as fo:
#     fo.write(data.content)

import pickle

with open('gene2go.pickle', 'rb') as f:
    dic = pickle.load(f)
    for key in dic.items():
        print(key)