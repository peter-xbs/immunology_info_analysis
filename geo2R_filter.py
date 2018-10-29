# _*_ coding:utf-8 _*_

import pickle
import requests

with open('gene2go.pickle', 'rb') as f:
    gene2go_dic = pickle.load(f)
    print(len(gene2go_dic))

# url = 'https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1539095023&job_key=JSID_01_592643_130.14.22.10_9000_geo2r'
#
# res = requests.get(url).content.decode()
# with open('GSE60290.diff', 'w') as fo:
#     fo.write(res)

key_words = ['mitochondrial']
scaned_set = set()
with open('GSE60290_filter.diff', 'w') as fo:
    with open('GSE60290.diff') as f:
        i = 0
        header_line = f.readline()
        fo.write(header_line)
        # header_dic = {}
        # for idx, item in enumerate(header_line.split('\t')):
        #     header_dic[item] = idx
        for line in f:
            flag = False
            line = line.replace('"', '')
            line_list = line.split('\t')[1:-1]
            gene_symbol = line_list[-1]
            if not gene_symbol.strip():
                continue
            if not gene_symbol in gene2go_dic:
                i += 1
                print(gene_symbol)
                continue

            go_infos = gene2go_dic[gene_symbol]
            for go in go_infos:
                for key in key_words:
                    if key in go:
                        flag = True
                        break
            if not flag:
                continue

            if gene_symbol in scaned_set:
                continue

            scaned_set.add(gene_symbol)

            new_line = '\t'.join(line_list)+'\n'
            fo.write(new_line)
    print(i)






