# _*_ coding:utf-8 _*_

import requests

# data = requests.get('https://www.frontiersin.org/articles/file/downloadfile/317773_supplementary-materials_presentations_1_zip/octet-stream/Presentation%201.ZIP/1/317773')
#
# with open('tmp.zip', 'wb') as fo:
#     fo.write(data.content)

import pickle

#
# with open('illum_mus_ref8_v2.anno') as f:
#     with open('illum_mus_ref8_v2.pickle', 'wb') as fo:
#         dic = {}
#         for line in f:
#             line_list = line.strip().split('\t')
#             if len(line_list) < 2:
#                 continue
#             ID, gene = line_list[0], line_list[1]
#             dic[ID] = gene
#         pickle.dump(dic, fo)

with open('config/gene2go.pickle', 'rb') as f:
    dic = pickle.load(f)
    for key in dic:
        print(len(dic[key]))

