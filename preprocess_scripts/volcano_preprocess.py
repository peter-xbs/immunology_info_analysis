# _*_ coding:utf-8 _*_

import pickle
import requests
import os


# url = 'https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1539095023&job_key=JSID_01_592643_130.14.22.10_9000_geo2r'
#
# url = 'https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1541686008&job_key=JSID_01_614458_130.14.22.10_9000_geo2r'
# res = requests.get(url).content.decode()
# with open('GSE43075.diff', 'w') as fo:
#     fo.write(res)

key1 = "glycolytic process" #糖酵解
key2 = "gluconeogenesis" # 糖异生
key3 = "lipid metabolic process" # 脂代谢
key4 = "cellular amino acid metabolic process" # 氨基酸代谢过程
key5 = "carbohydrate metabolic process" # 碳酸代谢过程
key6 = "glutamate metabolic process " # 谷氨酸代谢
key7 = "glucose metabolic process" # 葡萄糖代谢过程
key8 = "reactive oxygen species metabolic process" # ROS 代谢
key9 = "one-carbon metabolic process" # 单碳循环
key10 = "glutamine metabolic process" #
key11 = "NAD metabolic process"
key_words = [key1, key2, key3, key4, key5, key6, key7, key8, key9, key10, key11]

home_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_dir = os.path.join(home_dir, 'VOLCANO_HEATMAP_PREPROCESS')
input_name = 'GSE43075'
input_diff = os.path.join(home_dir, input_name+'.diff')
output_diff = os.path.join(output_dir, input_name+'.process.diff')

go_config = os.path.join(home_dir, 'config/gene2go.pickle')
with open(go_config, 'rb') as f:
    gene2go_dic = pickle.load(f)


# annotation
with open(output_diff, 'w') as fo:
    with open(input_diff) as f:
        i = 0
        scaned_set = set()
        # header_dic = {}
        # for idx, item in enumerate(header_line.split('\t')):
        #     header_dic[item] = idx
        for idx, line in enumerate(f):
            supp_list = [] # 追加各个基因是否属于某个信号通路
            flag = False
            line = line.replace('"', '')
            line_list = line.split('\t')[1:-1]

            if idx == 0:
                supp_list = ['key'+str(i) for i in range(1, len(key_words)+1)]

            else:
                gene_symbol = line_list[-1]
                if not gene_symbol.strip():
                    continue
                if not gene_symbol in gene2go_dic:
                    go_infos = ('', '', '')
                else:
                    go_infos = gene2go_dic[gene_symbol]

                for key in key_words:
                    flag = 'N'
                    for go in go_infos:
                        if key in go:
                            flag = 'Y'
                            break
                    supp_list.append(flag)


                if gene_symbol in scaned_set:
                    continue

                scaned_set.add(gene_symbol)

            new_line = '\t'.join(line_list+supp_list)+'\n'
            fo.write(new_line)
    print(i)