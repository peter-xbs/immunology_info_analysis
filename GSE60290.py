# _*_ coding:utf-8 _*_

import requests

# url = 'https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1538533969&job_key=JSID_01_588503_130.14.22.10_9000_geo2r'
#
# data = requests.get(url).content.decode()
#
# with open('gse60290.txt', 'w') as fo:
#     fo.write(data)
#
with open('GSE53986_series_matrix_format_filt.txt') as f:
    gene_set = set()
    for line in f:
        line_list = line.split('\t')
        gene = line_list[0]
        gene_set.add(gene)

with open('gse60290.txt') as f:
    for line in f:
        print(line)